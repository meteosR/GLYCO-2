#!/usr/bin/env python
# coding: utf-8

# Version 3:
#    - Added cython support for worker function
#    - count points argument construction is now done inside the multiprocessing worker
#        - Reduced memory usage for explicit argument construction
#        - Reduced time and CPU usage, now nproc is proportional to CPU usage

# TODO:
#     - Add arg parse support
#     - Add freesasa support !!
#     - Add output folder, do all calculations there
#     - Add epitope definition of protein residues
#     - Actually calculate aggregate metrics like in glyco after collecting results
#     - Benchmark KD-tree
#     - Benchmark RAM usage:
#         - Save intermediate results for large bundles if memory consumption is too large
#     - Handle multiprocessing failures
#         - Failure to spawn process
#         - Failure inside a child: report and start again ?
#         - Failure to spawn process
#     - Add general progress bars, keep track of children progress
#     - PDB view support (GUI)
#     - fixed radius query for d?


# For pdb extraction
from itertools import *
from Bio import PDB, SeqIO, SeqUtils

# For multiprocessing
import multiprocessing
from multiprocessing.pool import Pool

from natsort import natsort_keygen

# Data structures and math
from sklearn.neighbors import KDTree
import numpy as np
import pandas as pd
from collections import defaultdict

# For profiling
import time
import os
import re

# Cython worker functions to speed up critical parts
from bundle_cylinder_calc import c_worker

import warnings

warnings.filterwarnings("ignore")


def log(file_in, msg):
    with open(file_in, "a") as f:
        f.write(str(msg) + "\n")


def renumber_pdb_glycans(file_name, wd, prefix, chain_names=None, type_names=[], by="chain"):
    with open(file_name, "r") as f:
        lines = f.readlines()

    out_name = "{}/result/{}_renumbered.pdb".format(wd, prefix)

    chains_num_old = defaultdict(lambda: None)
    chains_num_new = defaultdict(lambda: 1)

    with open(out_name, "w") as f:
        for line in lines:
            if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
                resname = line[17:21].strip()
                chain = line[21].strip()
                resid = line[22:28].strip()

                if chains_num_old[chain] is None:
                    chains_num_old[chain] = resid

                if by == "chain":

                    if chain in chain_names:

                        if chains_num_old[chain] != resid:
                            chains_num_new[chain] += 1
                            chains_num_old[chain] = resid

                        new_num = str(chains_num_new[chain]).center(6, " ")
                        new_line = line[:22] + new_num + line[28:]
                    else:
                        new_line = line

                if by == "type":

                    if resname in type_names:

                        if chains_num_old[chain] != resid:
                            chains_num_new[chain] += 1
                            chains_num_old[chain] = resid

                        new_num = str(chains_num_new[chain]).center(6, " ")
                        new_line = line[:22] + new_num + line[28:]
                    else:
                        new_line = line
                f.write(new_line)
            else:
                f.write(line)

    return out_name


def get_struct(file_name):
    parser = PDB.PDBParser()
    struct = parser.get_structure('test_1', file_name)
    return struct


def get_atoms_data(file_name, type_names, surface_set=None):
    atoms_data = []

    # Glycan and protein sets will be selected by residue name, line by line

    with open(file_name, "r") as f:
        lines = f.readlines()

    for line in lines:

        resname = line[17:21].strip()

        if resname in type_names:

            chain = line[21].strip()
            resid = line[22:28].strip()
            atomtype = line[12:16].strip()
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except:
                continue

            key = chain, resname, resid

            if (surface_set is None) or (surface_set is not None and key in surface_set):

                if "H" not in atomtype:
                    # Get atom coordinates
                    coords = np.array([x, y, z], dtype=np.float32).reshape(3, )  # .astype(dtype=float32)

                    # Create id for atom
                    atom_id = "{}_{}_{}_{}".format(chain, resname, resid, atomtype)

                    # Append data
                    atoms_data.append((coords, atom_id))

    return np.array(atoms_data, dtype=object)


def get_surface_keys(file_name, wd, prefix, fresasa_path, surface_threshold, protein_alphabet, probe_radius):
    ABS_data = defaultdict(float)

    # Run freesasa
    print("Running freesasa", flush=True)
    outrsa = "{}/result/{}.rsa".format(wd, prefix)
    os.system("{} {} -w -p {} --rsa > {}".format(fresasa_path, file_name, probe_radius, outrsa))

    print("{} {} -p 1.4 --rsa > {}".format(fresasa_path, file_name, outrsa), flush=True)

    # Read rsa output file
    regex = re.compile("RES\s+([A-Z]{3})\s+([A-Z])\s*([a-zA-Z]*\s*\d+\w*)\s*(\d*[.,]?\d*)\s*")

    surface_set = []
    with open(outrsa, "r") as f:
        for line in f.readlines():
            try:
                m = regex.match(line)
                if m is not None:
                    residue, chain, pos, value_ABS = m.groups()
                    value_ABS = float(value_ABS)

                    if value_ABS >= surface_threshold and residue in protein_alphabet:
                        surface_set.append((chain, residue, pos))
                        key = "_".join([chain, residue, pos])
                        ABS_data[key] = value_ABS
            except:
                pass

    if os.path.isfile(outrsa):
        os.remove(outrsa)

    return set(surface_set), ABS_data


def get_subset(residue_list_file, protein_alphabet):
    surface_set = []

    with open(residue_list_file, "r") as f:
        lines = f.readlines()

    for line in lines:
        line = line.split()
        if len(line) == 3:
            # LYS A  121
            resn = line[0]
            chain = line[1]
            resid = line[2]
            if resn in protein_alphabet:
                surface_set.append((chain, resn, resid))
    return set(surface_set)


def apply_bfactor(df_res_counts, in_name, wd, out_file):
    # Make bfactor pdb
    bfactor_mapper = df_res_counts
    bfactor_mapper["protein_res_id"] = bfactor_mapper["Protein Chain"] + "_" + bfactor_mapper["Protein residue"] + "_" + \
                                       bfactor_mapper["Protein residue position"]
    bfactor_mapper = bfactor_mapper.set_index("protein_res_id")["Glycan_density"].to_dict()
    print("Applying bfactor")

    with open(in_name, "r") as f:
        lines = f.readlines()

    out_name = "{}/result/{}".format(wd, out_file + "_bfactor.pdb")

    print(out_name)
    with open(out_name, "w") as f:
        for line in lines:
            if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
                # atomnum = line[6:11].strip()
                # atomtype = line[12:16].strip()
                resname = line[17:21].strip()
                chain = line[21].strip()
                resid = line[22:28].strip()

                atom_id = "{}_{}_{}".format(chain, resname, resid)

                if atom_id in bfactor_mapper:
                    b = bfactor_mapper[atom_id]
                    val = str('{:>6.2f}'.format(b))
                    new_line = line[:60] + val + line[66:]
                else:
                    val = str('{:>6.2f}'.format(0))
                    new_line = line[:60] + val + line[66:]

                f.write(new_line)
            else:
                f.write(line)


def merge_dict_lists(added_dict=None, inplace_dict=None, val_unique=True):
    for k, v in added_dict.items():
        if k in inplace_dict.keys():
            if val_unique:
                inplace_dict[k] = list(np.unique(inplace_dict[k] + v))
            else:
                inplace_dict[k] += v
        else:
            inplace_dict[k] = v


def process_dataframes(all_results, ABS_data):
    # Retrieve results from calls
    print("\nCollecting results...", flush=True)
    all_results = {k: v.get() for k, v in all_results.items()}

    print("\nMerging temp dicts...", flush=True)

    final_dict = {}
    for final_dict_i in all_results.values():
        merge_dict_lists(added_dict=final_dict_i, inplace_dict=final_dict, val_unique=True)

    print("Making dataframe...", flush=True)

    final_df = defaultdict(list)
    all_v = []
    for k, v in final_dict.items():
        final_df["Protein_ID"].append(k)
        final_df["Glycans_atoms"].append(v)
        all_v = all_v + v

    duplicate_sum = len(all_v)
    non_duplicate_sum = len(list(np.unique(all_v)))

    final_df = pd.DataFrame(final_df)

    final_df["Glycan_density"] = final_df["Glycans_atoms"].apply(len)
    final_df['Protein Chain'], final_df['Protein residue'], final_df['Protein residue position'] = zip(
        *final_df['Protein_ID'].apply(lambda key: key.split("_")))
    final_df["Protein residue position"] = final_df["Protein residue position"].astype(dtype=str)

    final_df["SASA ABS"] = final_df["Protein_ID"].apply(lambda key: ABS_data[key] if key in ABS_data else np.nan)

    final_df = final_df.drop("Glycans_atoms", axis=1)[["Protein Chain", "Protein residue position", "Protein residue", "Glycan_density", "SASA ABS"]]

    # final_df = final_df.sort_values(["Protein Chain", "Protein residue position"])

    final_df = final_df.sort_values(
        by=["Protein Chain", "Protein residue position"],
        key=natsort_keygen()
    )

    print("\nFinished df processing...", flush=True)

    return final_df, duplicate_sum, non_duplicate_sum


############################################
# ################ Main ################## #
############################################


def main_glyco(file_name, wd, fresasa_path, protein_types=None, glycan_types=None, r=1.0, distance_cutoff=26.0,
               surface_threshold=30.0, probe_radius=1.4, nproc=32, residue_list_file=None, module_type="all_atom"):
    log_file = wd + "/log.txt"

    print("Starting\n", flush=True)
    log(log_file, "Started main_glyco")

    e_msg = ""

    start_time = time.perf_counter()

    ############################################
    # ############### Main code ############## #
    ############################################

    out_file = os.path.basename(file_name).replace(".pdb", "")
    print(out_file, flush=True)

    print("Using a pool with {} cpus.".format(nproc), flush=True)
    bundle_num = nproc * 2

    # Renumber glycan seqNum from 1 to n on each chain
    file_name = renumber_pdb_glycans(file_name, wd, out_file, chain_names=None, type_names=glycan_types, by="type")

    # Parse pdb
    struct = get_struct(file_name)

    # Get surface residues from RSA
    if module_type == "all_atom":
        print("Calculating surface residues.", flush=True)
        log(log_file, "Calculating surface residues.")
        surface_set, ABS_data = get_surface_keys(file_name, wd, out_file, fresasa_path, surface_threshold,
                                                 protein_types, probe_radius)

        if len(surface_set) == 0:
            return None, None, "Surface set is empty"

    elif module_type == "subset":
        log(log_file, "Getting subset")
        print("Module sub")
        path_sub = residue_list_file
        print(path_sub)
        surface_set = get_subset(path_sub, protein_types)
        ABS_data = {}

        if len(surface_set) == 0:
            return None, None, "Surface set is empty"

    else:
        surface_set = None
        ABS_data = {}

    # # Get atomic data for glycan, protein # #

    print("\nExtracting glycan atoms.", flush=True)
    log(log_file, "Extracting glycan atoms.")
    print("glycan_types", glycan_types)
    glycan_data_list = get_atoms_data(file_name, type_names=glycan_types, surface_set=None)

    print("\nExtracting protein atoms.", flush=True)
    log(log_file, "Extracting protein atoms.")
    protein_data_list = get_atoms_data(file_name, type_names=protein_types, surface_set=surface_set)

    print("\nFound {} glycan heavy atoms.".format(len(glycan_data_list)), flush=True)
    print("Found {} protein heavy atoms with RSA > {}.".format(len(protein_data_list), surface_threshold), flush=True)

    log(log_file, "Found {} glycan heavy atoms.".format(len(glycan_data_list)))
    log(log_file, "Found {} protein heavy atoms.".format(len(protein_data_list)))

    # Get protein heavy atoms that can be counted on cylinder collisions
    log(log_file, "Extracting collision atoms.")
    protein_collision_data_list = get_atoms_data(file_name, type_names=protein_types, surface_set=None)

    print("\nExtracted collision atoms.", flush=True)

    if len(glycan_data_list) == 0:
        return None, None, "No protein atoms were found."
    if len(protein_data_list) == 0:
        return None, None, "No glycan atoms were found."
    if len(protein_collision_data_list) == 0:
        return None, None, "No protein atoms were found."

    log(log_file, "Building kd tree.")
    all_coords = np.stack(protein_collision_data_list[:, 0]).astype(dtype=np.float64)

    print("\nBuilt coordinate matrix.", flush=True)

    # Add coordinates to KD-Tree for fast radius queries in 3D space
    tree = KDTree(all_coords, leaf_size=3)

    print("\nBuilt kd-tree.", flush=True)

    # # Preprocess data to send to multiple CPUs # #

    # Set the longest of {protein, glycan} lists to be the outter loop when generating pairs, and partition outter
    # into sublists to send to multipel CPUs with full inner list
    if len(glycan_data_list) > len(protein_data_list):
        outter_data = glycan_data_list
        inner_data = protein_data_list
        glycan_out = True
    else:
        outter_data = protein_data_list
        inner_data = glycan_data_list
        glycan_out = False

    # Split calculations into multiple bundles
    bundle_outter_data = np.array_split(outter_data, bundle_num, axis=0)

    set_ = "Glycan" if glycan_out else "Protein"
    print("Each bundle has {} atoms from set={} to calculate. There are {} bundles.".format(
        bundle_outter_data[0].shape[0], set_, bundle_num), flush=True)

    # Store results
    all_results = {}

    # Define a pool: nproc CPUs will compute simultaneously at all times, and taking care of all jobs in a queue
    log(log_file, "Getting pool.")
    tp = Pool(nproc)

    # Loop through all our jobs and add them to the queue
    for job_id, outter_data in enumerate(bundle_outter_data):

        if glycan_out:
            glycan_mini = outter_data
            protein_mini = inner_data
        else:
            glycan_mini = inner_data
            protein_mini = outter_data

        # Add worker to queue
        res = tp.apply_async(c_worker, (job_id, glycan_mini, protein_mini, r, distance_cutoff, all_coords, tree))

        # Save AsyncResult object
        all_results[job_id] = res

        # Close and join pool, waits for all jobs to be added to queue on main thread
    print("Added all jobs to queue.\n", flush=True)
    log(log_file, "Added all jobs to queue")
    tp.close()
    log(log_file, "Closed pool.")
    tp.join()
    log(log_file, "Joined pool")

    # Process all results from processes, make dataframes and final result files
    log(log_file, "Processing dataframes.")
    df_res_counts, duplicate_sum, non_duplicate_sum = process_dataframes(all_results, ABS_data)

    # Cleanup
    if os.path.isfile(file_name):
        os.remove(file_name)

    end_time = time.perf_counter()
    log(log_file, "Done calculating. Time={} min.".format(round((end_time - start_time) / 60, 2)))
    print("Done. Execution time=", (end_time - start_time) / 3600, "hours.", flush=True)

    return df_res_counts, duplicate_sum, non_duplicate_sum, struct, e_msg


"""
conda create -n GLYCO python=3.8
conda activate GLYCO
pip install -r requirements.txt
python3 -m pip install Cython && python3 setup.py build_ext --inplace
"""
