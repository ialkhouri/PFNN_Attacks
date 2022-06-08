import pickle
import numpy as np
import re
import random
import os


# The below pickle file must be placed in the same directory
pickle_file = open("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/blosum62_fixed_matrix.pickle", 'rb')

matrix62_similarities = pickle.load(pickle_file)
# this list matches the the label of each row/column to the matrix62_similarities
mapper = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

residue_abbre_dect = {"ALA": "A",
                      "ARG": "R",
                      "ASN": "N",
                      "ASP": "D",
                      "CYS": "C",
                      "GLN": "Q",
                      "GLU": "E",
                      "GLY": "G",
                      "HIS": "H",
                      "ILE": "I",
                      "LEU": "L",
                      "LYS": "K",
                      "MET": "M",
                      "PHE": "F",
                      "PRO": "P",
                      "SER": "S",
                      "THR": "T",
                      "TRP": "W",
                      "TYR": "Y",
                      "VAL": "V"}

print("")

# makre sure you have the "blosum62_fixed_matrix.pickle" file in the same location as this fucntion

def BLOSUM62_distance_score(seq_1,seq_2):
    """

    :param seq_1: string 1
    :param seq_2: string 2
    :return: BLOSUM sciore based on BLOSUM62 fixed matrix. String 1 and string 2 must be of the same length!!!. Strings MUST be upper case
    """



    distance_per_residue = []

    for res_index in range(len(seq_1)):

        row_res_term1 = seq_1[res_index]
        row_ind_term1 = mapper.index(row_res_term1)

        col_res_term1 = seq_1[res_index]
        col_ind_term1 = mapper.index(col_res_term1)

        row_res_term2 = seq_1[res_index]
        row_ind_term2 = mapper.index(row_res_term2)

        col_res_term2 = seq_2[res_index]
        col_ind_term2 = mapper.index(col_res_term2)

        distance_per_residue.append(matrix62_similarities[row_ind_term1][col_ind_term1] -
                                    matrix62_similarities[row_ind_term2][col_ind_term2])

    return sum(distance_per_residue)


def exact_Distane_geometry_rep(Aseq1, Aseq2):
    """

    :param Aseq_1: list of tuples. each tuple represent the x,y,z in the paper
    :param Aseq_2: list of tuples. each tuple represent the x,y,z in the paper
    :return: the exact distance, and the RMSD... The exact distance reflects distance and the RMSD is similarity.. High RMSD, low D ==> adv example is no good and PFNN is robust... Low RMSD, high D ==> adv example is good and PFNN is not robust
    """
    d_Aseq1 = np.zeros(shape=(len(Aseq1), len(Aseq1)))
    d_Aseq2 = np.zeros(shape=(len(Aseq2), len(Aseq2)))

    diff_mat_Aseq = np.zeros(shape=(len(Aseq2), len(Aseq2)))

    for i in range(len(Aseq1)):
        for j in range(len(Aseq1)):
            if i != j:
                d_Aseq1[i][j] = (Aseq1[i][0] - Aseq1[j][0]) ** 2 + \
                                (Aseq1[i][1] - Aseq1[j][1]) ** 2 + \
                                (Aseq1[i][2] - Aseq1[j][2]) ** 2

                d_Aseq2[i][j] = (Aseq2[i][0] - Aseq2[j][0]) ** 2 + \
                                (Aseq2[i][1] - Aseq2[j][1]) ** 2 + \
                                (Aseq2[i][2] - Aseq2[j][2]) ** 2

                diff_mat_Aseq[i][j] = (d_Aseq1[i][j] - d_Aseq2[i][j]) ** 2

    ## get the exact  Distane_geometry_ representation == \sum_{i,j\in [n]} ( d_Aseq1(i,j) - d_Aseq2(i,j) )^2
    # Distane_geometry_rep =
    exact_Distane_geometry_rep = np.sum(diff_mat_Aseq)

    #RMSD = (1/exact_Distane_geometry_rep)*(10**(10))

    D_str = exact_Distane_geometry_rep*(10**(-10))


    return exact_Distane_geometry_rep, D_str

# This fucntion is still under construction
def Find_farther_residues(seq_1):
    """

    :param seq_1: input sequence
    :return: farthest residues res_a and res_b
    """
    current_max = 0
    for i in range(len(seq_1)):
    #for res in list(seq_1):
        res = seq_1[i]
        # this is to obtain the value from the blosum matrix
        res_ind_term_1 = mapper.index(res)
        # this is to get the index in the sequence

        #for other_res in list(seq_1):
        for j in range(len(seq_1)):
            other_res = seq_1[j]
            other_ind_term_1 = mapper.index(other_res)

            if res != other_res:

                # calculate distance is residues are not the same
                distance_temp = matrix62_similarities[res_ind_term_1][other_ind_term_1]
                # DOUBLE CHECK HERE!
                distance_temp = np.abs(distance_temp)

                # if distance is max, replace, else don't
                if distance_temp >= current_max:
                    current_max = distance_temp
                    # save res_a and res_b
                    res_a = res
                    res_b = other_res
                    index_a = i
                    index_b = j

                # logger
                #print("distance between", [res, other_res], " ; Indicies = ", [j, j], " = ", distance_temp,      " ; current max = ", current_max)

    return res_a, res_b, index_a, index_b

# This fucntion is still under construction
def Approximate_Distance_geometry(seq_1, seq_2, Aseq1, Aseq2):

    res_a, res_b, index_a, index_b = Find_farther_residues(seq_1)

    res_a_hat, res_b_hat, index_a_h, index_b_h = Find_farther_residues(seq_2)

    ### from Aseq1, s_a and s_b, obtain the value of d(s_a,s_b)
    # tuple from Aseq that represent res_a
    A_s_a = Aseq1[index_a]
    # tuple from Aseq that represent res_b
    A_s_b = Aseq1[index_b]

    d_Aseq1_ab = (A_s_a[0] - A_s_b[0]) ** 2 + \
                 (A_s_a[1] - A_s_b[1]) ** 2 + \
                 (A_s_a[2] - A_s_b[2]) ** 2

    ### from Aseq2, s_a_hat and s_b_hat, obtain the value of d(s_a_hat,s_b_hat)
    # tuple from Aseq that represent res_a_hat
    A_s_a_h = Aseq2[index_a_h]
    # tuple from Aseq that represent res_b_hat
    A_s_b_h = Aseq2[index_b_h]

    d_Aseq1_ab_h = (A_s_a_h[0] - A_s_b_h[0]) ** 2 + \
                   (A_s_a_h[1] - A_s_b_h[1]) ** 2 + \
                   (A_s_a_h[2] - A_s_b_h[2]) ** 2

    D_max_structure = (d_Aseq1_ab - d_Aseq1_ab_h) ** 2

    return D_max_structure


def PDB_file_to_structure_coordinates(file_path):
    """
    :param file_path: string of the pdb file path
    :return: tuple for CA 3D co-ordinates
    This is a fucntion that takes a pdb file and output a tuple of of co-ordinates with length equal to the sequence length. We read CA
    """

    pdb_file = open(file_path, 'r')

    Lines = pdb_file.readlines()

    output_list_of_tuples = []

    ## we need the number of lines
    for line in Lines:
        line_string = line
        ### we need to only consider the lines starting with ATOM and contains CA
        if 'ATOM' in line_string and "CA" in line_string:
            #print(line)
            ### look for co-ordinates in that line and add it as a tuple in the output lise
            p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
            floats = [float(i) for i in p.findall(line_string)]
            # Convert strings to float and append it to the list
            # print(floats)
            output_list_of_tuples.append((floats[0], floats[1], floats[2]))

    return output_list_of_tuples


def Sequence_from_PDB_file(file_path):
    """
    :param file_path: string of the pdb file path
    :return: sequence
    This is a fucntion that takes a pdb file and output a sequence
    """

    pdb_file = open(file_path, 'r')

    Lines = pdb_file.readlines()

    output_seq_string = ''

    ## we need the number of lines
    for line in Lines:
        line_string = line
        ### we need to only consider the lines starting with ATOM and contains CA
        if 'ATOM' in line_string and "CA" in line_string:

            res = line_string.split("CA")[1][2:5]
            res_abb = residue_abbre_dect[res]

            output_seq_string = output_seq_string + res_abb



    return output_seq_string



def fasta_to_seq_string(input_fasta_file_path):
    """

    :param input_fasta_file_path:
    :return: seq str
    """
    fasta_file_current = open(input_fasta_file_path, 'r')

    Lines = fasta_file_current.readlines()

    seq_oneLiner_list = []

    for line in Lines:
        if ">" not in line:
            # remove \n if it has one:
            if "\n" in line:
                line = line[0:-1]
            seq_oneLiner_list.append(line)

    ## convert the list of characters to one string
    seq_string = ""
    for item in seq_oneLiner_list:
        seq_string = seq_string + item

    return seq_string


def adv_seq_gen(indicies_to_switch, clean_sequence):
    """

    :param indicies_to_switch: the indices to change randomly from the clean sequence
    :param clean_sequence:
    :return: adv_sequence and BLOSUm_distance
    """
    mapper = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    def lists_Diff(li1, li2):
        return list(set(li1) - set(li2)) + list(set(li2) - set(li1))

    residues_list_to_switch = list(clean_sequence[indicies_to_switch[0]:indicies_to_switch[-1]+1])
    indices_to_switch = indicies_to_switch


    # for each residue in the list_to_change, get a list of all possibilities
    ### get a sequence by replacing each residue in the residues_list_to_switch by a randomn choice residue from the possible choices
    Adv_seq_list = []
    cntr = 0

    for index in range(len(clean_sequence)):
        #random.seed(seed_id+index)
        if index in indices_to_switch:
            Adv_seq_list.append(random.choice(lists_Diff(mapper, [residues_list_to_switch[cntr]])))
            cntr = cntr + 1
        else:
            Adv_seq_list.append(clean_sequence[index])

    BLOSUm_distance = BLOSUM62_distance_score(clean_sequence, Adv_seq_list)

    return Adv_seq_list, BLOSUm_distance


def RMSD_between_two_structures(file_org_path, file_adv_path):
    """

    :param file_org_path: path of the org structure pdb file
    :param file_adv_path: path of the adv structure pdb file
    :return: Given two .pdb files, rename one of them, create pml file, input pymol commands, run the .pml file, get results...
    """
    # renaming

    os.system("cp " + file_adv_path + " " + file_adv_path[0:-12] + "ranked_0_adv.pdb")

    file_adv_path_renamed = file_adv_path[0:-12] + "ranked_0_adv.pdb"

    # creating .pml file and writeing the commands... Go to the PyMOL documentation for this
    pml_file = open("temp_pml_file.pml", "w")
    pml_file.write("load " + file_org_path + "\n")
    pml_file.write("load " + file_adv_path_renamed + "\n")
    pml_file.write("align " + "ranked_0, " + "ranked_0_adv, " + "cycles=0")
    #pml_file.write("align " + "ranked_0, " + "ranked_0_adv, ")
    #pml_file.write("alignto ")
    pml_file.close()

    # executing the .pml file
    out_str = os.popen('pymol -c temp_pml_file.pml').read()

    # extract the data
    for item in out_str.split("\n"):
        if "Executive: RMSD =" in item:
            output_RMSD_str = item[21:]

    output_RMSD = float(output_RMSD_str.split("(")[0])

    return output_RMSD

def pdb_file_name_from_path(file_path):
    """

    :param file_path:
    :return: .pdb file name
    We need this so we can use align in PyMOL without
    """
    return file_path.split("/")[-1]


def RMSD_between_two_structures_raw(file_org_path, file_adv_path):
    """

    :param file_org_path: path of the org structure pdb file
    :param file_adv_path: path of the adv structure pdb file
    :return: Given two .pdb files, without rename one of them, create pml file, input pymol commands, run the .pml file, get results...
    """

    # creating .pml file and writeing the commands... Go to the PyMOL documentation for this
    pml_file = open("temp_pml_file.pml", "w")
    pml_file.write("load " + file_org_path + "\n")
    pml_file.write("load " + file_adv_path + "\n")

    file_org_pdb_only = pdb_file_name_from_path(file_org_path)
    file_adv_pdb_only = pdb_file_name_from_path(file_adv_path)

    pml_file.write("align " +file_org_pdb_only[:-4]+ ", " +file_adv_pdb_only[:-4]+", " + "cycles=0")

    #pml_file.write("align " + "ranked_0, " + "ranked_0_adv, " + "cycles=0")

    #pml_file.write("align " + "ranked_0, " + "ranked_0_adv, ")
    #pml_file.write("alignto ")
    pml_file.close()

    # executing the .pml file
    out_str = os.popen('pymol -c temp_pml_file.pml').read()

    # extract the data
    for item in out_str.split("\n"):
        if "Executive: RMSD =" in item:
            output_RMSD_str = item[21:]

    output_RMSD = float(output_RMSD_str.split("(")[0])

    return output_RMSD




def PyMOL_align_Ismail(org_pdb_file, adv_pdb_file):
    """

    :param org_pdb_file: original pdb file to be the target for alignment
    :param adv_pdb_file: adversarial pdb file to be aligned w.r.t. the org
    :return: aligned adv pdb file path. NOTE THAT THE ORDER MATTERS...
    """
    # make a new directory to save
    os.system("mkdir " + adv_pdb_file[0:-12] + "aligned" + " 2> /dev/null")

    # renaming

    os.system("cp " + adv_pdb_file + " " + adv_pdb_file[0:-12] + "ranked_0_adv_tobe_aligned.pdb")

    renamed_adv_pdb_file_to_be_aligned = adv_pdb_file[0:-12] + "ranked_0_adv_tobe_aligned.pdb"

    # creating .pml file and writeing the commands... Go to the PyMOL documentation for this
    pml_file = open("temp_pml_file_for_Alignning_and_resaving.pml", "w")
    pml_file.write("load " + renamed_adv_pdb_file_to_be_aligned + "\n")
    pml_file.write("load " + org_pdb_file + "\n")

    # pml_file.write("align " + "ranked_0, " + "ranked_0_adv_tobe_aligned, " + "cycles=0" + "\n")
    pml_file.write("align " + "ranked_0_adv_tobe_aligned, " + "ranked_0, " + "cycles=0" + "\n")
    # pml_file.write("align " + "ranked_0, " + "ranked_0_adv, ")
    # pml_file.write("alignto ")

    # remove original
    # pml_file.write("remove " + org_pdb_file + "\n")
    # pml_file.write("remove " + org_pdb_file + "\n")

    # save the adv as
    pml_file.write("save " + adv_pdb_file[0:-12] + "aligned/ranked_0_adv__aligned.pdb")

    pml_file.close()

    # [this is only used tfor debugging] reading the output of the pymol file
    out_str = os.popen('pymol -c temp_pml_file_for_Alignning_and_resaving.pml').read()

    adv_pdb_aligned_file_path = adv_pdb_file[0:-12] + "aligned/ranked_0_adv__aligned.pdb"

    return adv_pdb_aligned_file_path


def GDT_score_fucntion(original_pdb_file_path, adv_aligned_pdb_file_path, GDT_type = "TS"):
    """

    :param original_pdb_file_path:
    :param adv_aligned_pdb_file_path:
    :return: GDT_score
    """
    # input: original PDB file, adv aligned PDB file
    # output GDT-HA and GDT-TS GDT_HA and GDT_TS similarity fucntion: [0-100]\%

    """
    GDT score = 100 * (C1 + C2 + C3 + C4) / 4N

    // Where:
    //   C1   = Count of number of residues superposed below (threshold/4)
    //   C2   = Count of number of residues superposed below (threshold/2)
    //   C3   = Count of number of residues superposed below (threshold
    //   C4   = Count of number of residues superposed below (2*threshold)
    //   N    = Total number of residues

    // GDT_TS (Total Score)   : threshold = 4
    // GDT_HA (High Accuracy) : threshold = 2

    source: http://www.sbg.bio.ic.ac.uk/~maxcluster/index.html#GDT
    """
    org_coords = PDB_file_to_structure_coordinates(original_pdb_file_path)

    adv_coords_after_align = PDB_file_to_structure_coordinates(adv_aligned_pdb_file_path)

    length_of_seq = len(org_coords)

    #### get differences and based on the thresholds (cutoffs), get the counters C_1, C_2, C_3, and C_4

    #GDT_type = "HA"

    if GDT_type == "TS":
        threshold = 4
    if GDT_type == "HA":
        threshold = 2

    C_1_counter = 0
    C_2_counter = 0
    C_3_counter = 0
    C_4_counter = 0

    difference_list_in_Angs = []

    for index in range(length_of_seq):

        d = np.sqrt((org_coords[index][0] - adv_coords_after_align[index][0]) ** 2 +
                    (org_coords[index][1] - adv_coords_after_align[index][1]) ** 2 +
                    (org_coords[index][2] - adv_coords_after_align[index][2]) ** 2)

        difference_list_in_Angs.append(d)

        if d < threshold / 4:
            C_1_counter = C_1_counter + 1
        if d < threshold / 2:
            C_2_counter = C_2_counter + 1
        if d < threshold:
            C_3_counter = C_3_counter + 1
        if d < 2 * threshold:
            C_4_counter = C_4_counter + 1

    GDT_score = (100 / (4 * length_of_seq)) * (C_4_counter + C_3_counter + C_2_counter + C_1_counter)

    return GDT_score
