
import numpy as np
cimport numpy as np
import pandas as pd
from collections import defaultdict

def c_count_points( np.ndarray[np.float32_t, ndim=1] pt1, 
                    np.ndarray[np.float32_t, ndim=1] pt2, 
                    np.float32_t r, 
                    np.float32_t vec_norm, 
                    np.ndarray[np.float64_t, ndim=2] all_coords,  
                    tree):
                           
    # Type defs  
    cdef np.ndarray[np.float32_t, ndim=1] vec
    cdef np.ndarray[np.float32_t, ndim=1] sphere_origin 
    cdef np.float32_t sphere_R 
    cdef np.ndarray[np.float64_t, ndim=2] query_points
    cdef np.float32_t const 
    cdef np.ndarray[double, ndim=2] q_pt1    
    

    # We only need to query for points within the sphere of radius R and origin (pt2 + pt1) / 2 that encapsulates
    # the cylinder defined by pt1, pt2 and radius r.
    # R = sqrt(r^2 + l^2 / 4)

    vec = pt2 - pt1

    sphere_origin = (pt1 + pt2) / 2
    sphere_R = np.sqrt(r ** 2 + vec_norm ** 2 / 4)


    # Get the points that might reside in our cylinder
    query_points = all_coords[tree.query_radius(sphere_origin.reshape(1, -1), r=sphere_R, count_only=False)[0]]

    # Actually check points in cylinder
    const = r * np.linalg.norm(vec)
    q_pt1 = query_points - pt1
    

    num =  np.count_nonzero((np.dot(q_pt1, vec) > 0) & (np.dot(query_points - pt2, vec) < 0) & (np.linalg.norm(np.cross(q_pt1, vec), axis=1) <= const)) 
    return num


def c_worker( np.uint32_t job_id, 
              glycan_data_list, 
              protein_data_list, 
              np.float32_t r, 
              np.float32_t distance_cutoff, 
              np.ndarray[np.float64_t, ndim=2] all_coords, 
              tree):

    cdef np.float32_t vec_norm

    # print("Job id={} starting.".format(job_id))

    res_dict = defaultdict(list)

    for idx, glycan_data in enumerate(glycan_data_list):
        for protein_data in protein_data_list:

            coords_glycan, glycan_id = glycan_data
            coords_protein, protein_id = protein_data
            
            chain, res, pos, atom = protein_id.split("_")
            key = "_".join([chain, res, pos])

            # Calculate distance between pair
            vec_norm = np.linalg.norm(coords_glycan - coords_protein)

            # Only add pair if distance < cutoff
            if vec_norm < distance_cutoff:
                num_points = c_count_points(coords_glycan, coords_protein, r, vec_norm, all_coords, tree)
                #res["glycan_id"].append(glycan_id)
                #res["protein_id"].append(protein_id)
                #res["num_points"].append(num_points)
                
                if num_points == 0:
                    res_dict[key].append(glycan_id)

    """
    mini_df = pd.DataFrame(res_dict)
    
    mini_df = mini_df.loc[mini_df["num_points"] == 0].copy()

    mini_df ['chain_name_pro'], mini_df ['residue_name_pro'], mini_df ['residue_pos_pro'], mini_df ['atom_name_pro'] = zip(*mini_df ['protein_id'].apply(lambda key: key.split("_"))) 
    mini_df ["residue_pos_pro"] = mini_df ["residue_pos_pro"].astype(dtype=str)
    
    mini_df = mini_df.groupby(["chain_name_pro", "residue_name_pro", "residue_pos_pro"])["glycan_id"].agg(list).reset_index()
    mini_df["glycan_id"] = mini_df["glycan_id"].apply(lambda l: np.unique(np.array(l)))
    """
    
    for k,v in res_dict.items():
        res_dict[k] = list(np.unique(v))

    # print("Job id={} finished.".format(job_id))

    return res_dict


