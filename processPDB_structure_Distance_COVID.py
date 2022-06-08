import csv
import glob
import re
from pandas import *
import numpy as np
import time

import pickle

import scipy.io


import subprocess

from PFNN_attacks_Functions import PDB_file_to_structure_coordinates, exact_Distane_geometry_rep, RMSD_between_two_structures, PyMOL_align_Ismail, GDT_score_fucntion
import os


################################################################################################################################################
#################################### Calucultaing the RMSD after alignmet using pymol ####################################
################################################################################################################################################

#### Logic is as follow:
##### create a fucntion, given two .pdb files, rename one of them, create pml file, input pymol commands, run the .pml file, get results...


file_org_path = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/sequence_COVID19_Q92985/org_sequence_COVID19_Q92985/ranked_0.pdb"

file_adv_path = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/sequence_COVID19_Q92985/adv_18_seq_COVID19_Q92985/ranked_0.pdb"

output_RMSD = RMSD_between_two_structures(file_org_path, file_adv_path)

print("output RMSD = ", output_RMSD)


# def RMSD_between_two_structures(file_org_path, file_adv_path):
#     """
#
#     :param file_org_path: path of the org structure pdb file
#     :param file_adv_path: path of the adv structure pdb file
#     :return: Given two .pdb files, rename one of them, create pml file, input pymol commands, run the .pml file, get results...
#     """
#     # renaming
#
#     os.system("cp " + file_adv_path + " " + file_adv_path[0:-12] + "ranked_0.pdb")
#
#     file_adv_path_renamed = file_adv_path[0:-12] + "ranked_0.pdb"
#
#     # creating .pml file and writeing the commands
#     pml_file = open("temp_pml_file.pml", "w")
#     pml_file.write("load " + file_org_path + "\n")
#     pml_file.write("load " + file_adv_path_renamed + "\n")
#     pml_file.write("align " + "ranked_0, " + "ranked_0_adv, " + "cycles=0")
#     pml_file.close()
#
#     # executing the .pml file
#     out_str = os.popen('pymol -c temp_pml_file.pml').read()
#
#     # extract the data
#     for item in out_str.split("\n"):
#         if "Executive: RMSD =" in item:
#             output_RMSD_str = item[21:]
#
#     output_RMSD = float(output_RMSD_str.split("(")[0])
#
#     return output_RMSD


print("__^__")


##################################################################
#################################
#################################
## Here, we test the fucntions and calculates the output structure distance - COVID results
#################################
##################################################################
#################################

############# Logic: for every out folder, get the seq name, and the original seq output.
############# Then, loop over every other folder that has "adv". For every one, get the ranked_0
############# log D_str, and at the end print the average...



dir_arr = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/")

number_of_seqs_soFar_cntr = 0

min_GDT_TS_perProtein_list = []
avg_GDT_TS_perProtein_list = []

min_GDT_HA_perProtein_list = []
avg_GDT_HA_perProtein_list = []

max_RMSD_perProtein_list = []
avg_RMSD_perProtein_list = []

max_Dstr_perProtein_list = []
avg_Dstr_perProtein_list = []

for folder_name in dir_arr:
    #### if the folder name has sequence
    if "sequence_" in folder_name:

        number_of_seqs_soFar_cntr = number_of_seqs_soFar_cntr + 1

        #print("*********** We are handelling", folder_name, "************")

        seq_name = folder_name[17:]

        original_seq_pdb =  "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/"+"sequence_COVID19_"+seq_name+"/org_sequence_COVID19_"+seq_name+"/ranked_0.pdb"
        struc_clean = PDB_file_to_structure_coordinates(original_seq_pdb)



        in_dir_arr = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/"+"sequence_COVID19_"+seq_name+"/")

        #################### from in_dir_arr, open commandsfixed.csv to get the sum of run times of running the org seq and 20 adv sequences...
        csv_file = read_csv("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/"+"sequence_COVID19_"+seq_name+"/"+"commandsfixed.csv")

        ### remove strings and convert to int

        runTimeCol_to_list = csv_file["runtimeseconds"].tolist()

        item_to_remove = []
        for index in range(len(runTimeCol_to_list)):
            # if item is a string
            if not isinstance(runTimeCol_to_list[index],int):

                # if item is a string and can be converted to integer, then convert and replace
                if runTimeCol_to_list[index].isdigit():
                    runTimeCol_to_list[index] = int(runTimeCol_to_list[index])

                # if item is a string and not a convertable to int, then remove it
                else:
                    item_to_remove.append(runTimeCol_to_list[index])
                    #runTimeCol_to_list.remove(runTimeCol_to_list[index])

        for item in item_to_remove:
            if item in runTimeCol_to_list:
                runTimeCol_to_list.remove(item)




        runTime_in_seconds = np.sum(runTimeCol_to_list)
        runTime_in_days    = runTime_in_seconds / (3600*24)


        max_D_str = 0
        max_adv = ""
        list_Dstr_to_save = []

        list_RMSD_to_save = []
        max_RMSD = 0
        max_RMSD_adv = ""

        list_GDT_TS_to_save = []
        min_GDT_TS = 100
        min_GDT_TS_adv = ""

        list_GDT_HA_to_save = []
        min_GDT_HA = 100
        min_GDT_HA_adv = ""


        for in_folder_name in in_dir_arr:
            ### if folder is adv
            if "adv" in in_folder_name:
                adv_seq_pdb = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/"+"sequence_COVID19_"+seq_name+"/"+in_folder_name+"/ranked_0.pdb"
                struc_adver = PDB_file_to_structure_coordinates(adv_seq_pdb)

                # This is to calcalate D_str
                _, D_str = exact_Distane_geometry_rep(struc_clean, struc_adver)

                RMSD = RMSD_between_two_structures(original_seq_pdb , adv_seq_pdb)

                # calculate GDT:

                adv_seq_pdb_ALIGNED = PyMOL_align_Ismail(original_seq_pdb, adv_seq_pdb)

                # maybe, we need a delay here to allow the file to be created... maybe for 1 second...
                time.sleep(1.5)

                GDT_TS = GDT_score_fucntion(original_seq_pdb, adv_seq_pdb_ALIGNED, GDT_type="TS")

                GDT_HA = GDT_score_fucntion(original_seq_pdb, adv_seq_pdb_ALIGNED, GDT_type="HA")




                list_Dstr_to_save.append(D_str)

                list_RMSD_to_save.append(RMSD)

                list_GDT_TS_to_save.append(GDT_TS)

                list_GDT_HA_to_save.append(GDT_HA)

                if D_str > max_D_str:
                    max_D_str = D_str
                    max_adv = in_folder_name

                if RMSD > max_RMSD:
                    max_RMSD = RMSD
                    max_RMSD_adv = in_folder_name

                if GDT_TS < min_GDT_TS:
                    min_GDT_TS = GDT_TS
                    min_GDT_TS_adv = in_folder_name

                if GDT_HA < min_GDT_HA:
                    min_GDT_HA = GDT_HA
                    min_GDT_HA_adv = in_folder_name

                #print("D_str = ", D_str, 'adv seq: ', [in_folder_name], "MAX so far = ", max_D_str , "[RMSD = ", [RMSD], "GDT_TS = ", [GDT_TS])
        #print("[seq, seq length] = ",[seq_name],[len(struc_clean)], 'Avg = ', np.mean(list_Dstr_to_save), "MAX = ", [max_D_str,max_adv], "runTime in Days = ", [runTime_in_days])

        # This logger to be copy paste to overleaf text... [Add RMSD and its average]
        print(seq_name, " & ", len(struc_clean), " & ",  round((len(struc_clean)-5)*100/len(struc_clean),4)   ,
              " & ", round(max_RMSD,4) , " & ", round(np.mean(list_RMSD_to_save),4) ,
              " & ", round(min_GDT_TS, 4), " & ", round(np.mean(list_GDT_TS_to_save), 4),
              " & ", round(min_GDT_HA, 4), " & ", round(np.mean(list_GDT_HA_to_save), 4),
              " & ", round(runTime_in_days, 4) , " \\\ ")
        print("\\hline")

        # # detailed logger
        # print("[seq, seq length] = ", [seq_name], [len(struc_clean)], '[Avg(Dstr), Avg(RMSD), Avg(GDT_TS)]  = ', [np.mean(list_Dstr_to_save), np.mean(list_RMSD_to_save), np.mean(list_GDT_TS_to_save)],
        #       "[Max Dstr, Max RMSD, Min GDT_TS] = ", [max_D_str,max_RMSD, min_GDT_TS], "===== Best  adv index = ", [max_adv, max_RMSD_adv, min_GDT_TS_adv],
        #       "=====runTime in Days = ", [runTime_in_days])
        # print('*********************************************************************************')

        min_GDT_TS_perProtein_list.append(min_GDT_TS)
        avg_GDT_TS_perProtein_list.append(np.mean(list_GDT_TS_to_save))

        min_GDT_HA_perProtein_list.append(min_GDT_HA)
        avg_GDT_HA_perProtein_list.append(np.mean(list_GDT_HA_to_save))

        max_RMSD_perProtein_list.append(max_RMSD)
        avg_RMSD_perProtein_list.append(np.mean(list_RMSD_to_save))

        max_Dstr_perProtein_list.append(max_D_str)
        avg_Dstr_perProtein_list.append(np.mean(list_Dstr_to_save))




print("The number of processed COVID19 proteins  =  ", number_of_seqs_soFar_cntr, "Avg of the minimums GDT_TS = ", [np.mean(min_GDT_TS_perProtein_list)])

np.save("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/output_measures/GDT_TS.npy", [min_GDT_TS_perProtein_list, avg_GDT_TS_perProtein_list])
np.save("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/output_measures/GDT_HA.npy", [min_GDT_HA_perProtein_list, avg_GDT_HA_perProtein_list])
np.save("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/output_measures/RMSD.npy", [max_RMSD_perProtein_list, avg_RMSD_perProtein_list])
#np.save("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/output_measures/Dstr.npy", [max_Dstr_perProtein_list, avg_Dstr_perProtein_list])

print("")


# Pi_best_class_specific_cifar10_matricies_PGDT = np.load("/home/ismail/pycharmProjects/SSLTL_project/Coarse_Robust/Pi_best_class_specific_cifar10_matricies_PGDT.npy",allow_pickle=True)
# scipy.io.savemat('/home/ismail/pycharmProjects/SSLTL_project/Coarse_Robust/Pi_best_class_specific_cifar10_matricies_PGDT.mat',    dict(Pi_best_class_specific_cifar10_matricies_PGDT=Pi_best_class_specific_cifar10_matricies_PGDT))
