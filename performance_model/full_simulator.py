import importlib

from simulator_imports import *
from os.path import exists
import pickle

## Get list of spectra, if its already present just load into object,
## otherwise re-read it from ms2 file

# Check if we have already read the spectra
if exists('spectra_list_file'):
    print("File exists")
    obj_file = open('spectra_list_file', 'rb')
    list_spectra = pickle.load(obj_file)
    obj_file.close()
else:
    print("File does not exist")
    # Open file to store the object
    obj_file = open('spectra_list_file', 'wb')
    list_spectra = []
    create_spec_index(list_spectra, 'combined-small.ms2')
    pickle.dump(list_spectra, obj_file)
    obj_file.close()

## Get all the peptides from human_proteome database

# Check if we have already read the peptides from the file
if exists('peptide_list_file'):
    print("Peptide File exists")
    obj_file = open('peptide_list_file', 'rb')
    list_peptides = pickle.load(obj_file)
    obj_file.close()
else:
    print("File does not exist")
    # Open file to store the object
    obj_file = open('peptide_list_file', 'wb')
    list_peptides = []
    for line in pep_file:
        values = line.split()
        list_peptides.append(peptide(values[0], values[1]))
    list_peptides.sort(key=lambda x: x.mass)
    pickle.dump(list_peptides, obj_file)
    obj_file.close()

print(len(list_peptides))



## Run the performance model code to estimate the latency of the system
precursor_tolerance = 0.1
total_spectra = 64
num_pcores = 2
num_pes = 2
spectrum_batch_size = 8
peptide_batch_size = 32
performance_model(list_spectra, list_peptides, precursor_tolerance, total_spectra, num_pcores, num_pes, spectrum_batch_size, peptide_batch_size)
