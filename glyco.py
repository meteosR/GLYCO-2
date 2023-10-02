import pandas as pd
import os
import sys
import pathlib
import argparse
import glob
from core_glyco import main_glyco, apply_bfactor, log
from utils import data_validator
import multiprocessing
from concurrent.futures import ProcessPoolExecutor as Pool_out
import time
import datetime

#######################################
# ########### USER INPUTS ########### #
#######################################

# how to run:
# python3 glyco.py -pdb glyco_examples/1_ebola_man5_frame_10__BGLN_BMAN_AMAN__23.pdb -cutoff 23 -glycans BGLN,BMAN,AMAN -out_folder 1_ebola_man5_frame_10 -ncpu 12 -module all_atom
# python3 glyco.py -pdb glyco_examples/2_ha_man5_frame_10__BGL_BMA_AMA__23.pdb -cutoff 23 -glycans BGL,BMA,AMA -out_folder 2_ha_man5_frame_10__BGL_BMA_AMA__23 -ncpu 12 -module all_atom
# python3 glyco.py -in_folder /afs/csail.mit.edu/u/m/mreveiz/data_mnt_rgb/glyco/GLYCO-2/glyco_examples/small_md -cutoff 23 -sur_cutoff 70  -glycans BGL,BMA,AMA -out_folder debug -ncpu 16 -nproc_out 2 -module all_atom -average


parser = argparse.ArgumentParser()
parser.add_argument('-pdb', type=str, help='Enter your name of pdb. ex) input.pdb')
parser.add_argument('-cutoff', type=float, help='Enter your distance cutoff for glycans in Angstrom. ex) 30',
                    default=23)
parser.add_argument('-cyl_radius', type=float, help='Enter the cylinder radius', default=1.4)
parser.add_argument('-probe_radius', type=float, help='BSA probe radius', default=1.4)
parser.add_argument('-glycans', type=str,
                    help='Enter your name of glycans with comma "," for separators. ex) BMA,AMA,BGLN', required=True)
parser.add_argument('-module', type=str, help='Enter your module name, either all_atom or subset. ex) all_atom',
                    required=True)
parser.add_argument('-residue', type=str, help='Enter your protein residue list. ex) residuelist.txt', default=None)
parser.add_argument('-freesasa', type=str, help='Enter your path of Freesasa executable, ex) dependencies/freesasa',
                    default="dependencies/freesasa")
parser.add_argument('-ncpu', type=int, help='Enter the number of workers to use for the computation ex) 9',
                    default=min(multiprocessing.cpu_count(), 8))
parser.add_argument('-nproc_out', type=int, help='Enter the number of parallel PDBs to process) 1',
                    default=1)
parser.add_argument('-sur_cutoff', type=float, help='Enter surface area cutoff of FreeSASA in Angstrom^2. ex) 40',
                    default=30)
parser.add_argument('-average', action='store_true',
                    help='Add if using input folder with multiple pdbs of the same protein.', default=False)
parser.add_argument('-in_folder', type=str, help='Input folder where input pdb files are located. ex) input',
                    default=None)
parser.add_argument('-out_folder', type=str, help='Output folder where results will be saved. ex) output',
                    required=True)
parser.add_argument('-verbose', action='store_true', help='How much logging info to be displayed', default=False)
args = parser.parse_args()

submission_data = {}

if args.in_folder is None:
    submission_data["pdb_files"] = [args.pdb]
else:
    in_folder = args.in_folder
    submission_data["pdb_files"] = list(glob.glob(os.path.join(in_folder, "*.pdb")))

protein_alphabet = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HSD', 'HIS', 'HIE', 'HID', 'HSE',
                    'HSP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

submission_data["glycan_names"] = args.glycans
submission_data["module_type"] = args.module
submission_data["residue_list_file"] = args.residue
submission_data["distance_cutoff"] = float(args.cutoff)
submission_data["cylinder_radius"] = float(args.cyl_radius)
submission_data["surface_threshold"] = float(args.sur_cutoff)
submission_data["protein_alphabet"] = protein_alphabet
submission_data["probe_radius"] = float(args.probe_radius)
submission_data["nproc_out"] = int(args.nproc_out)
submission_data["verbose"] = verbose = args.verbose
fresasa_path = args.freesasa
nproc = args.ncpu
out_path = args.out_folder
pdb_files = submission_data["pdb_files"]
average_pdbs = args.average

#######################################
# ########### Validation ############ #
#######################################

valid_submission, error_message_validation = data_validator(submission_data)

submission_data["glycan_names"] = [x.strip() for x in submission_data["glycan_names"]]

if not valid_submission:
    print("Error: ", error_message_validation)
    sys.exit(1)

if not pathlib.Path(fresasa_path).exists():
    print("Freesasa path not found.")
    sys.exit(1)

if len(pdb_files) <= 1 and average_pdbs:
    print("Cannot average pdbs when only 1 PDB is provided.")
    sys.exit(1)

if pathlib.Path(out_path).exists():
    print("Output folder '{}' already exists. Please prodive another folder.".format(out_path))
    sys.exit(1)
else:
    os.mkdir(out_path)
    os.mkdir(out_path + "/result")
    wd = out_path
    print("Output folder '{}' was created.".format(out_path))

##############################################################

# Save input params
param_file = wd + "/params_in.txt"
for k, v in submission_data.items():
    if "folder" not in k:
        log(param_file, str(k) + ": " + str(v))


##############################################################


def run_glyco(file_name, wd, fresasa_path, nproc, submission_data):
    if verbose: print("About to start analysis for " + str(file_name), flush=True)

    res = main_glyco(file_name, wd, fresasa_path, protein_alphabet, submission_data["glycan_names"],
                     submission_data["cylinder_radius"], submission_data["distance_cutoff"],
                     submission_data["surface_threshold"], submission_data["probe_radius"], nproc,
                     submission_data["residue_list_file"], submission_data["module_type"], verbose)

    df_res_counts, duplicate_sum, non_duplicate_sum, struct, e_msg = res

    # Write pdb with bfactors
    out_file = os.path.basename(file_name).replace(".pdb", "")
    df_res_counts.to_csv(out_path + "/result/final_df_{}.csv".format(out_file), sep=",", index=False)  # new
    apply_bfactor(df_res_counts, file_name, wd, out_file, verbose)
    if verbose: print("Finished apply_bfactor", flush=True)

    # Write sums
    with open(out_path + "/result/glysums_{}.txt".format(out_file), "w") as f:
        f.write(str(duplicate_sum) + "\n" + str(non_duplicate_sum) + "\n")


def main():
    if len(pdb_files) == 1:
        # Single pdb to calculate
        run_glyco(pdb_files[0], wd, fresasa_path, nproc, submission_data)
    else:
        # Multiple pdbs

        # Create pool
        nproc_out = submission_data["nproc_out"]
        OUTER_POOL = Pool_out(nproc_out)

        # Submit all jobs
        for input_file in pdb_files:
            OUTER_POOL.submit(run_glyco, input_file, wd, fresasa_path, nproc, submission_data)

        # Shutdown and wait for children to finish
        OUTER_POOL.shutdown(wait=True)

        # Average frames if requested
        if average_pdbs:
            if verbose: print("Averaging frames...")
            all_results_df = pd.DataFrame()

            output_files = list(glob.glob(os.path.join(wd + "/result/", "final_df_*.csv")))
            for idx, o_file in enumerate(output_files):
                df = pd.read_csv(o_file, sep=",")
                df["frame"] = idx
                all_results_df = pd.concat((all_results_df, df), ignore_index=True)

                all_results_df_averaged = \
                all_results_df.groupby(["Protein Chain", "Protein residue", "Protein residue position"])[
                    ["Glycan_density", "SASA ABS"]].mean().reset_index()
                all_results_df_averaged.to_csv(out_path + "/result/final_df_averaged.csv", sep=",", index=False)

    print("Done processing!", flush=True)


start = time.time()
print("Starting glyco!")
main()
end = time.time()
print('TIME SPENT: ', str(round(datetime.timedelta(seconds=end - start).total_seconds() / 60, 2))), " minutes."
