# Exploring The Predictive Limit of AlphaFold Again BLOSUM Based Adversarial Sequences

- This repo is for the paper submitted to ICML2023

- We provide the original and adversarial sequences pushed to AlphFold. We use AlphaFold V2 located at "https://github.com/deepmind/alphafold". Samples of the .fasta files used are given inside the Sequences directory. 

- We generate sequences using "Find_input_search_space_ofseqs_examples.py". The comments inside the code serve as the guide to generate the adversarial sequences. 

- We process the output PDB files using "processPDB_structure_Distance_COVID.py". Two samples are given inside the "out" directory. Due, to space constraints, we will provide the complete list of the PDB files (for the reduced and full database AlphaFold configuration) upon acceptance. Make sure to unzip the folder inside "out". 

- The complete original and adversarial structures for the AlphaFold full configuration are given inside the Aliged_Protein_Structures folder. 

- The requirements to run the scripts are numpy, panda, and PyMOL. 





