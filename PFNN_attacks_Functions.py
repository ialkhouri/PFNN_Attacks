import pickle
import numpy as np
import re

# The below pickle file must be placed in the same directory
pickle_file = open("/home/ismail/anaconda3/envs/fold_attakcs_toy_env/blosum62_fixed_matrix.pickle", 'rb')

matrix62_similarities = pickle.load(pickle_file)
# this list matches the the label of each row/column to the matrix62_similarities
mapper = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

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

    RMSD = exact_Distane_geometry_rep*(10**(-10))


    return exact_Distane_geometry_rep, RMSD


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
