import random
import time
import helper_functions as hf
import numpy as np
import pypdb as pyp
import importlib
hf = importlib.reload(hf)

# Here we define a simple class for structure objects, in which we store all information about
# a given protein
class Structure:
    def __init__(self, name):
        self.name = name
        self.ligands = set()
        self.subunits = set()
        self.functions = set()
        self.processes = set()

    def Add_Ligand(self, ligand_list):
        self.ligands = set(ligand_list)

    def Add_Subunits(self, subunits):
        self.subunits = set(subunits)

    def Add_Functions(self, function):
        self.functions = set(function)

    def Add_Processes(self, process):
        self.processes = set(process)

#Create dictonaries to keep all valid proteins and all proteins we have seen begore
Seen_Dict = {}
Structure_Dict = {}
file_save_name = "ProteinDict_ten_thousand"
seen_save_name = "Seen_Dict"
Seen_Dict = hf.readDict(seen_save_name, Seen_Dict)
Structure_Dict = hf.readDict(file_save_name, Structure_Dict)
#Create a safe save before we start modifying the protein dictonary, as some times it can get messed up
hf.saveDict(Structure_Dict, "ProteinDict_ten_thousand_Safe_Save")

# %% generate list of random structures
# set seed value
seed_val = 5
random.seed(seed_val)
# set sample size
sample_size = 10000
# get random sample of structures from database
num_structs = len(pyp.get_all())
sample_indxs = np.array(random.sample(range(num_structs), sample_size))
rand_sructs = np.array(pyp.get_all())[sample_indxs]

#Here we actually mine for new proteins
start = time.time()
for index, struct in enumerate(rand_sructs):
    #If we have already seen a protein, just skip to the next one
    if struct in Seen_Dict:
        continue
    #Make a note that we have seen this protein now
    Seen_Dict[struct] = 1
    print("Number:", index, "Struct:", struct)
    #Save the seen structure dictonary
    if index % 1 == 0:
        hf.saveDict(Seen_Dict, seen_save_name)
    #If for whatever reason the protein is already in the structure dictonary, but not in the seen set, still skip
    if struct in Structure_Dict:
        continue
    #Get information on the current protein from the databasses
    res, ligs, functions, processes = hf.get_ligands_and_uniprot(struct)
    #If the protein contains no ligands, residues, functions, or processes, then we are not interested in it
    if ligs == None or res == None or len(functions) == 0 or len(processes) == 0:
        continue
    #Create a new structure object for the given protein and add it to the proein dictornary
    new_Struct = Structure(struct)
    new_Struct.Add_Subunits(res)
    new_Struct.Add_Ligand(ligs)
    new_Struct.Add_Functions(functions)
    new_Struct.Add_Processes(processes)
    Structure_Dict[struct] = new_Struct

    #Save the protein dictonary and
    if index % 1 == 0:
        hf.saveDict(Structure_Dict, file_save_name)
        hf.saveDict(Seen_Dict, seen_save_name)
end = time.time()

print('Completitoin TIME', end - start)
hf.saveDict(Structure_Dict, file_save_name)
