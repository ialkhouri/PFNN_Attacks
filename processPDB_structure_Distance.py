import glob
import re

import numpy as np

from PFNN_attacks_Functions import PDB_file_to_structure_coordinates, exact_Distane_geometry_rep
import os

#################


##################################################################
#################################
#################################
## Here, we test the fucntions and calculates the output structure distance - PDB results
#################################
##################################################################
#################################
#
#
#
# arr = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/Jha_sequences_with_Adv_T40/out4/")
#
# seq_name = '6OEU'
#
# PDB_file_clean = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_used_in_Jha_paper/out1/rcsb_pdb_"+seq_name+"/ranked_0.pdb"
# struc_clean = PDB_file_to_structure_coordinates(PDB_file_clean)
#
# list_toSave = []
#
# for folder_name in arr:
#     if seq_name in folder_name:
#         arr_folder = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/Jha_sequences_with_Adv_T40/out4/"+folder_name+'/')
#
#         for file in arr_folder:
#             if file == "ranked_0.pdb":
#                 PDB_file_Adver = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/Jha_sequences_with_Adv_T40/out4/"+folder_name+'/'+file
#
#
#                 struc_adver = PDB_file_to_structure_coordinates(PDB_file_Adver)
#
#                 _, RMSD = exact_Distane_geometry_rep(struc_clean, struc_adver)
#
#                 print("RMSD = ", RMSD, 'adv seq: ', [folder_name] )
#
#                 list_toSave.append(RMSD)
#
# print('Avg = ', np.mean(list_toSave))









#D_exact_struc, RMSD = exact_Distane_geometry_rep(struc_clean, struc_adver)


# ##################################################################
# #################################
# #################################
# ## Here, we test the fucntions and calculates the output structure distance - yeast results
# #################################
# ##################################################################
# #################################
#
#
#
# arr = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/Jha_sequences_with_Adv_T20/T20all/")
#
# seq_name = 'SFP1'
#
# YEAST_file_clean = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_used_in_Jha_paper/out1/S288C_YLR403W_SFP1_genomic/ranked_0.pdb"
# struc_clean = PDB_file_to_structure_coordinates(YEAST_file_clean)
#
# list_toSave = []
#
# for folder_name in arr:
#     if seq_name in folder_name:
#         arr_folder = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/Jha_sequences_with_Adv_T20/T20all/"+folder_name+'/')
#
#         for file in arr_folder:
#             if file == "ranked_0.pdb":
#                 PDB_file_Adver = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/Jha_sequences_with_Adv_T20/T20all/"+folder_name+'/'+file
#
#
#                 struc_adver = PDB_file_to_structure_coordinates(PDB_file_Adver)
#
#                 _, RMSD = exact_Distane_geometry_rep(struc_clean, struc_adver)
#
#                 print("RMSD = ", RMSD, 'adv seq: ', [folder_name] )
#
#                 list_toSave.append(RMSD)
#
# print('Avg = ', np.mean(list_toSave))
#
# #PDB_file_Adver = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/pdb_output_file_samples/2ROP_adv1_dbbc0_relaxed_rank_1_model_2.pdb"
#
# print("")

##################################################################
#################################
#################################
## Here, we test the fucntions and calculates the output structure distance - COVID results
#################################
##################################################################
#################################

COVID_file_clean = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/P59632/ranked_0.pdb"
struc_clean = PDB_file_to_structure_coordinates(COVID_file_clean)
PDB_file_Adver = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/COVID19_sequences_with_Adv_T40/outT40/COVID19_P59632_adv_1/ranked_0.pdb"
struc_adver = PDB_file_to_structure_coordinates(PDB_file_Adver)



arr = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/COVID19_sequences_with_Adv_T40/outT40/")

seq_name = 'P07711'

COVID_file_clean = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/P07711/ranked_0.pdb"
struc_clean = PDB_file_to_structure_coordinates(COVID_file_clean)

list_toSave = []

for folder_name in arr:
    if seq_name in folder_name:
        arr_folder = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/COVID19_sequences_with_Adv_T40/outT40/"+folder_name+'/')

        for file in arr_folder:
            if file == "ranked_0.pdb":
                PDB_file_Adver = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/COVID19_sequences_with_Adv_T40/outT40/"+folder_name+'/'+file


                struc_adver = PDB_file_to_structure_coordinates(PDB_file_Adver)

                _, RMSD = exact_Distane_geometry_rep(struc_clean, struc_adver)

                print("RMSD = ", RMSD, 'adv seq: ', [folder_name] )

                list_toSave.append(RMSD)

print('Avg = ', np.mean(list_toSave))

#PDB_file_Adver = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/pdb_output_file_samples/2ROP_adv1_dbbc0_relaxed_rank_1_model_2.pdb"

print("")