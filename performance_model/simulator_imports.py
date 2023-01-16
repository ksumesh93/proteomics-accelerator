# package imports
import numpy as np
import struct
import math
import sys
import re
import pickle
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from random import random
from os.path import exists
import os, os.path
import bisect


# defining important functions and classes

class spectrum:
    def __init__(self, sid, mass):
        self.id = sid
        self.mass = mass

#Define peptide class, that is going to represent the peptide information
class peptide:
    def __init__(self, sequence, mass):
        self.sequence = sequence
        self.mass = float(mass)
    def __lt__(self, other):
        return self.mass < other.mass
    def __repr__(self):
        return self.sequence + ' ' + str(self.mass)
    
# Putting the imports and definitions here
# pep_file = open("human-peptides.txt", "r")
# pep_mod_file = open("modified_peptides", "a")

modify = {'S' : '*', 'T' : '+', 'Y' : '-'} 
        
def create_spec_index(list_spectra, file_name):
    # Open the ms2 file
    exp = open(file_name, 'r')

    count = 0
    list_spectra = []

    # Read all the lines
    for line in exp:

        if line[0] == 'S':
            words = line.split()
            spec_id = int(words[1])
            count += 1
            if count % 10000 == 1:
                print(count)
        elif line[0] == 'Z':
            words = line.split()
            mass = int(float(words[2])*1000)
            list_spectra.append(spectrum(spec_id, mass))
        
      # Sort the spectral indices by their mass
    list_spectra = sorted(list_spectra, key=lambda x: x.mass, reverse = False)
    exp.close()
    return


    
def phosphorylation(pept, place):
    mod = []
    
    mass = pept.mass
    queue = [(list(pept.sequence), mass)]
    visited = set()
    visited.add(pept.sequence)
    modcount = 0
    while queue:
        seq = queue.pop(0)
        chars = seq[0]
        mass = seq[1]
        for index,amino in enumerate(chars):
            if amino == 'S' or amino == 'T' or amino == 'Y':
                mod = chars.copy()
                mod[index] = modify[amino]
                conv = ''.join(mod)
                mass += 79.966
                #tup = (conv,mass)
                if conv not in visited:
                    queue.append((mod,mass))
                    p1 = peptide(conv, mass)
                    place.append(p1)
                    visited.add(conv)
                    modcount += 1
    return modcount


def mergeKLists(lists):
        amount = len(lists)
        interval = 1
        while interval < amount:
            for i in range(0, amount - interval, interval * 2):
                lists[i] = merge2Lists(lists[i], lists[i + interval])
                
            interval *= 2

        return lists[0] if amount > 0 else None

def merge2Lists(l1, l2):
    m = len(l1)
    n = len(l2)
    
    print("merging list of size: {}, {}".format(m,n))
    
    x = m+n
    
    l_merge = [None]*x
    i = j = k = 0
    while i < m and j < n:
        if l1[i].mass <= l2[j].mass:
        #if l1[i] <= l2[j]:
            l_merge[k] = l1[i]
            k += 1
            i += 1
        else:
            l_merge[k] = l2[j]
            k += 1
            j += 1
    
    if j == n:
        while i < m:
            l_merge[k] = l1[i]
            k += 1
            i += 1
    elif i == m:
        while j < n:
            l_merge[k] = l2[j]
            k += 1
            j += 1
    
    return l_merge

### A more accurate performance model code

# Writing a class for processing element
class proc_element:
    def __init__(self, pe_id, const_latency_factor):
        self.const_latency_factor = const_latency_factor
        self.pe_id = pe_id
        #self.pe_state = pe_state
        # Each PE will add up total computation time
        self.C_computation_time = 0
        # And total communication time
        self.C_communication_time = 0
        
        # Synchronization time
        self.T_PSync = 0
        
        # Cummulative time spent on peptide communication and spectrum computation windows
        self.T_PComm = 0
        self.T_SComm = 0
        self.T_PComp = 0
        self.T_SComp = 0

        # Time spent on previous peptide 
        self.T_Pcomp_W = 0
        self.T_Pcomm_W = 0
        # Time spent on previous spectrum
        self.T_Scomp_W = 0
        self.T_Scomm_W = 0
        
        # Detail of assigned spectrum
        self.spectrum_id = 0
        self.spectrum_mass = 0
        # Peptides assigned to this PE1
        self.peptide_start = 0
        self.peptide_end = 0
        self.peptide_length = 0
        self.current_index = 0
        self.total_peptides = 0
        self.num_peptide_batches = 0
        
    def assign_spectrum(self, sid, mass):
        self.T_PComm = 0
        self.T_PComp = 0
        self.spectrum_id = sid
        self.spectrum_mass = mass
        
    def assign_candidates(self, peptide_start, peptide_end):
        self.peptide_start = peptide_start
        self.peptide_end = peptide_end
        self.current_index = peptide_start
        
    def spectrum_communicate(self, num_pes):
        # Each spectrum has 128 mz/intensity pairs = 32 packets of 256-bit long (1024 bytes) = 32 cycles 
        # We also add the pipelining delay, the PE closest to the pipeline will get it quickest
        self.T_Scomm_W = 32 * num_pes
        
        # Accumulate this time
        self.T_SComm += self.T_Scomm_W
    
    def peptide_communicate(self, list_peptides, num_pes):
        # Check if we have reached the end of our batch
        if self.current_index <= self.peptide_end:
            self.peptide_length = len(list_peptides[self.current_index].sequence)
            self.current_index += 1
        else:
            self.T_Pcomm_W = 0
            return
        
        # Delay of sending peptide to PE (each peptide in one 256 bit packet) 
        self.T_Pcomm_W = num_pes
        
        # Accumulating the peptide computation time 
        self.T_PComm += self.T_Pcomm_W
        
        self.T_Scomp_W += self.T_Pcomp_W + self.T_Pcomm_W
        self.T_SComp = self.T_Scomp_W
    
    def peptide_compute(self):
        if self.current_index > self.peptide_end:
            self.T_Pcomp_W = 0
            return
        #The latency of our data-path is as many cycles as the length of the number of fragment-ions
        # We use only b/y ions so it will be twice the length of peptide
        self.T_Pcomp_W = 2 * self.peptide_length
        
        # Accumulating the peptide computation time 
        self.T_PComp += self.T_Pcomp_W 
    
    def spectrum_compute(self, list_peptides, num_pes):
        for j in range(self.peptide_start, self.peptide_end, 1):
            self.peptide_communicate(len(list_peptides[j].sequence), num_pes)
            self.peptide_compute()
        self.T_Scomp_W = self.T_PComp + self.T_PComm
        self.T_SComp += self.T_Scomp_W
        

# Writing a class for pcore
class pcore:
    def __init__(self, pcore_id, sbuffer_size, pbuffer_size, pe_num, precursor_tolerance):
        # This attribute tells how many spectra can be held not how many bytes
        # We will do the scaling later
        self.spectrum_batch_size  = sbuffer_size
        # This attribute tells how many peptides/scores can we hold
        self.peptide_batch_size = pbuffer_size
        # num of PEs
        self.pe_num = pe_num

        #Precursor tolerance
        self.precursor_tolerance = precursor_tolerance
        
        self.pcore_id = pcore_id
        #self.pe_state = pe_state
        # Each Pcore will add up total computation time
        self.T_computation_time = 0
        # And total communication time
        self.T_communication_time = 0
        
        # Cummulative time spent on batch communication and batch computation windows
        self.T_PBComm = 0
        self.T_SBComm = 0
        self.T_PBComp = 0
        self.T_SBComp = 0

        # Time spent on previous peptide 
        self.T_PBcomp_W = 0
        self.T_PBcomm_W = 0
        self.T_PBComm_ov = 0
        # Time spent on previous spectrum
        self.T_SBcomp_W = 0
        self.T_SBcomm_W = 0
        self.T_SBComm_ov = 0
        
        # Cummulative time spent on peptide communication and spectrum computation windows
        self.T_PComm = 0
        self.T_SComm = 0
        self.T_SComm_ov = 0
        self.T_PComm_ov = 0
        # Time spent on previous peptide 
        self.T_Pcomp_W = 0
        self.T_Pcomm_W = 0
        self.T_PComp = 0
        self.T_SComp = 0
        # Time spent on previous spectrum
        self.T_SComp_W = 0
        self.T_SComm_W = 0
        
        # Total number of spectra processed
        self.total_spectra = 0
        
        # Total dot-scores computed
        self.total_scores = 0
        
        # Create instances of the PE array
        self.pe_list = []
        for i in range(self.pe_num):
            self.pe_list.append(proc_element(i,1))
            
    def spectrum_batch_communicate(self):
        # There are 32 cycles in a spectrum, assuming it takes one cycle to copy a packet
        self.T_SBcomm_W = self.spectrum_batch_size * 32
        
        overhead = 0
        window_time = 0
        if self.T_SBcomm_W <= self.T_SBcomp_W:
            overhead = 0
        # Bad overlap
        elif self.T_SBcomm_W > self.T_SBcomp_W:
            overhead = self.T_SBcomm_W - self.T_SBcomp_W
        
        window_time = max([self.T_SBcomp_W, self.T_SBcomm_W])
        # Useful print messages
        print('T_SBcomp_W: {}, T_SBcomm_W: {}, SBComm_overhead: {}, Total window time: {}'.
              format(self.T_SBcomp_W, self.T_SBcomm_W, overhead, window_time))
        
        # Accumulate this time
        self.T_SBComm += self.T_SBcomm_W
        self.T_SBComp += window_time
        self.T_SBComm_ov += overhead
        return window_time
    
    def peptide_batch_communicate(self, batch_num):
        self.T_PBcomm_W = self.peptide_batch_size 
        
        # Can we overlap it with spectrum batch computation window
        
        # # The case of perfect over-lap
        # if self.T_PBcomm_W <= self.T_PBcomp_W:
        #     self.T_PBcomm_W = 0
        # # Bad overlap
        # elif self.T_PBcomm_W > self.T_PBcomp_W:
        #     self.T_PBcomm_W = self.T_PBcomm_W - self.T_PBcomp_W
        
        # Accumulate this time
        self.T_PBComm += self.T_PBcomm_W
    
    def peptide_batch_compute(self, list_peptides, batch_num):
        batch_size = int(self.peptide_batch_size/self.pe_num)
        batch_size = self.peptide_batch_size
        x = 0
        self.T_PBcomp_W = 0
        
        while x < batch_size:
            list_pcomputes = []
            list_pcommuns = []
            # Transfer a peptide to all the PEs
            for i in range(self.pe_num):
                self.pe_list[i].peptide_communicate(list_peptides, self.pe_num)
                list_pcomputes.append(self.pe_list[i].T_Pcomp_W)
                list_pcommuns.append(self.pe_list[i].T_Pcomm_W)
                x += 1
            
            self.T_PComp_W = max(list_pcomputes)
            self.T_PComm_W = max(list_pcommuns)
        
            if self.pe_num > self.T_PComp_W:
                overhead = self.pe_num - self.T_PComp_W
            else:
                overhead = 0
                
            window_time = max([self.T_PComp_W, self.T_PComm_W])
            
            # Useful print messages
            print('T_Pcomp_W: {}, T_Pcomm_W: {}, P_Comm overhead: {}, Total window time: {}'.
                  format(self.T_PComp_W, self.T_PComm_W, overhead, window_time))
            
            self.T_PComp += window_time
            self.T_PComm_ov += overhead
            
            # All PEs compute their respective peptide
            for i in range(self.pe_num):
                self.pe_list[i].peptide_compute()
            
            # Find the synchronization time
            for i in range(self.pe_num):
                self.pe_list[i].T_PSync += self.T_PComp_W  - self.pe_list[i].T_Pcomp_W
                
            # Accumulate this time
            self.T_PBcomp_W += window_time
            
        
        
    
    def spectrum_batch_compute(self, spectrum_start, list_spectra, list_peptides):
        # This is the max of the time taken by all the PEs 
        x = 0
        self.T_SBcomp_W = -1*self.T_SComp_W
        
        while x < self.spectrum_batch_size:
            list_spec_computes = []
            # Assigning all spectra to the PEs
            for i in range(self.pe_num):
                mass = list_spectra[spectrum_start + x +i].mass/1000
                sid = list_spectra[spectrum_start + x +i].id
                self.pe_list[i].assign_spectrum(sid, mass)
                
                # Assign every PE their candidate peptides
                start = bisect.bisect_left(list_peptides, peptide("AAAAA", mass - self.precursor_tolerance))
                end = bisect.bisect_left(list_peptides, peptide("AAAAA", mass + self.precursor_tolerance))
                
                k = end- start
                print("Transfer spectrum: {} (k = {}) from pcore: {} spectrum buffer to PE : {}".
                      format(spectrum_start + x +i, k,self.pcore_id, i))
                self.pe_list[i].spectrum_communicate(self.pe_num)
                
                list_spec_computes.append(self.pe_list[i].T_Scomp_W)
                
                
                self.pe_list[i].assign_candidates(start, end)
                
                peptide_batch_size = int(self.peptide_batch_size/self.pe_num)
                self.pe_list[i].num_peptide_batches = int((end-start)/peptide_batch_size)
            
            
            #self.T_SComp_W = max(list_spec_computes)
            self.T_SComm_W = 32*self.pe_num
            if self.T_SComm_W > self.T_SComp_W:
                overhead = self.T_SComm_W - self.T_SComp_W
            else:
                overhead = 0
            window_time = max([self.T_SComp_W, self.T_SComm_W])
            
            self.T_SBcomp_W += window_time
            self.T_SComp += window_time
            self.T_SComm_ov += overhead
            print('T_SComp_W: {}, T_SComm_W: {},  SComm overhead: {}, Total window time: {}'.
                  format(self.T_SComp_W, self.T_SComm_W, overhead, window_time))
            
            # Find the maximum number of batches
            max_batch = max([pe.num_peptide_batches for pe in self.pe_list])
            self.T_SComp_W = 0
            
            
            for i in range(max_batch):
                print("Transfer peptide batch: {} from CPU to peptide buffer in pcore: {}".format(i, self.pcore_id))
                self.peptide_batch_communicate(i)
                # Useful print messages
                
                if self.T_PBcomm_W > self.T_PBcomp_W:
                    overhead = self.T_PBcomm_W - self.T_PBcomp_W
                else:
                    overhead = 0
                window_time = max([self.T_PBcomp_W, self.T_PBcomm_W])
                self.T_PBComp += window_time
                
                print('T_PBcomp_W: {}, T_PBcomm_W: {},  PBComm overhead: {}, Total window time: {}'.
                      format(self.T_PBcomp_W, self.T_PBcomm_W, overhead, window_time))
                
                self.T_PBComm_ov += overhead
                print("Computing peptide batch: {} on pcore: {}".format(i, self.pcore_id))
                self.peptide_batch_compute(list_peptides, i)
                
                self.T_SComp_W += window_time
                
            self.T_PBComp += self.T_PBcomp_W    
            self.T_SComp_W += self.T_PBcomp_W
            
            x += self.pe_num
            
        self.T_SBcomp_W += self.T_SComp_W   
        
        return self.T_SBcomp_W


# Simulate the model based on data and hardware parameters
def performance_model(list_spectra, list_peptides, tolerance, num_spectra, m, n, S_SB, S_PB):
    precursor_tolerance = tolerance
    num_pes = n
    num_pcores = m

    pcore_list = []
    x = 0


    spectrum_batch_size = S_SB
    peptide_batch_size = S_PB

    T_ACW = [0]*num_pcores
    T_AComp = [0]*num_pcores

    while x < num_spectra:
        for i in range(num_pcores):
            # Create new object for PCORE if it doesn't exist
            if len(pcore_list) <= i:
                pcore_list.append(pcore(i, spectrum_batch_size, peptide_batch_size, num_pes, precursor_tolerance))

            T_ACW[i] = pcore_list[i].spectrum_batch_communicate()
            T_AComp[i] = pcore_list[i].spectrum_batch_compute(x, list_spectra, list_peptides)

            x += spectrum_batch_size

            if x >= num_spectra:
                for i in range(num_pcores):
                    pcore_list[i].T_SBComp += T_AComp[i]
        # for i in range(num_pcores):
        #     T_ACW[i] += T_AComp[i]


    for i in range(num_pcores):
        print("Stats for pcore: {}, Total time: {}, Computation time: {}, Communication time: {}".format(
                                                    i, pcore_list[i].T_SBComp, 
                                                pcore_list[i].T_SBComp - pcore_list[i].T_SBComm_ov, pcore_list[i].T_SBComm_ov))
        print("Peptide batch time: {}, Peptide batch compute time: {}, Peptide batch communicate time: {}".format(
                                                    pcore_list[i].T_PBComp, 
                                                pcore_list[i].T_PBComp - pcore_list[i].T_PBComm_ov, pcore_list[i].T_PBComm_ov))
        print("Spectrum processing time: {}, Spectrum compute time: {}, Spectrum communicate time: {}".format(
                                                    pcore_list[i].T_SComp, 
                                                pcore_list[i].T_SComp - pcore_list[i].T_SComm_ov, pcore_list[i].T_SComm_ov))
        print("Peptide processing time: {}, Peptide compute time: {}, Peptide communicate time: {}".format(
                                                 pcore_list[i].T_PComp, 
                                                pcore_list[i].T_PComp - pcore_list[i].T_PComm_ov, pcore_list[i].T_PComm_ov))
    

        