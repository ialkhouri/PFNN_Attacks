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


##### load the oveall data from the dictionary...
adv_winners = np.load("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/output_measures/adv_winners.npy",allow_pickle=True)
scipy.io.savemat('/home/ismail/pycharmProjects/SSLTL_project/Coarse_Robust/adv_winners.mat',    dict(adv_winners=adv_winners))

adv_winners = adv_winners.item()

dir_arr = os.listdir("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/")



FIGUREs_file = open("Figure_file.txt", "w")
for folder_name in dir_arr:
    #### if the folder name has sequence
    if "sequence_" in folder_name:

        seq_name = folder_name[17:]

        # this is the original .pdb file
        original_seq_pdb = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/" + "sequence_COVID19_" + seq_name + "/org_sequence_COVID19_" + seq_name + "/ranked_0.pdb"


        # LET THE ADV WINNDER BE 1 ... [HERE WE CODE THE ADV...]

        adv = int(adv_winners[seq_name][0][4])

        adv_winner_seq_pdb = "/home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/out/" + "sequence_COVID19_" + seq_name + "/adv_"+str(adv)+"_seq_COVID19_" + seq_name + "/ranked_0_adv.pdb"

        ## PyMOL commands:
        pml_file = open("temp_pml_file_for_PNG_saving.pml", "w")
        pml_file.write("load " + original_seq_pdb + "\n")
        pml_file.write("load " + adv_winner_seq_pdb + "\n")

        pml_file.write("align " + "ranked_0, " + "ranked_0_adv, " + "cycles=0" + "\n")

        # color the org as blue and adv as red:
        pml_file.write("color blue, ranked_0" + "\n")
        pml_file.write("color red, ranked_0_adv" + "\n")

        # save image command:

        pml_file.write("png /home/ismail/anaconda3/envs/fold_attakcs_toy_env/fasta_sequences_COVID19/output_figures/"+seq_name+".png, " +  "width = 5cm, height = 3cm, dpi = 600, ray=1")

        pml_file.close()

        # [this is only used tfor debugging] reading the output of the pymol file
        out_str = os.popen('pymol -c temp_pml_file_for_PNG_saving.pml').read()


        ########################## PRINTING OUT ...


        overleaf_string = "\\begin{figure*}[t]\\centering\\includegraphics[width=16cm]{img/supp_org_vs_adv_Imgs/" + \
                          seq_name + ".png}" + \
                          "\\vspace{-0.8cm}\\caption[]{{\\small{Sequence: " + seq_name + ", " + \
                          "$n=$" + str(adv_winners[seq_name][1]) + \
                          ", Similarity = " + str(adv_winners[seq_name][2]) + " \%" \
                          ", RMSD = " + str(adv_winners[seq_name][3]) + " \\r{A}. }}} \\end{figure*}"




        FIGUREs_file.write(overleaf_string + "\n")



        """
        \begin{figure*}[t]\centering\includegraphics[width=16cm]{img/supp_org_vs_adv_Imgs/O95992.png}\vspace{-0.8cm}\caption[]{{\small{Sequence: XXXXXXX, $n=000$, RMSD = 0000\r{A}, Sim = 0000\%}}}\end{figure*}
        """







        print("We are done with seq: ", [seq_name])

FIGUREs_file.close()