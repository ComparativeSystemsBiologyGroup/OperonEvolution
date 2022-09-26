import os
import math
import pathlib
import numpy as np

## Dictionaries and directories
working_dir = "D:\\Dottorato\\Operon_Evolution" ## You would need to change this accordingly!
os.chdir(working_dir)
modules_present = {"identifier":[], "description":[], "type":[]}
identifiers = []
module_types = []
with open("modules_present.tsv", "r") as f:
    for line in f:
        identifiers.append(line.split("\t")[0].rstrip())
        module_types.append(line.split("\t")[2].rstrip())

def wmean(array):
    mean = np.mean(array)
    weights = []
    for i in array:
        val = mean - i
        weights.append(val)
    return np.average(weights)

#=============================================================================
## The proximity score function
def calcProximity4(starts, ends, strands, size):
    quad = size/2
    starts.sort()
    ends.sort()
    distances = []
    overflow = 0
    for i in range(len(starts)-1):
        dist = starts[i+1]-ends[i]
        if dist > quad:
            overflow += 1
            continue
        distances.append(dist)
    #print("OVERFLOW:", overflow) ## Always one tops
    if overflow:
        distances.append(starts[0] + (size - ends[-1]))
    for i in range(len(distances)):
        if distances[i] < 1:
            distances[i] = 1
    return distances
#=============================================================================

total_distances = []
## The big LOOP
for species_path in pathlib.Path("D:\\Dottorato\\Operon_Evolution\\mono_species").glob("*"): ## You would need to change this accordingly!
    os.chdir(species_path)
    
    species = os.path.basename(species_path)
    kegg_genomes = os.listdir()
    num_genomes = len(kegg_genomes)
    nmodules = 0
    modules_average = []
    ## Now we must iterate in every folder containing module-level kegg genome information
    for kegg_path in pathlib.Path(species_path).glob("*"):
        os.chdir(kegg_path)
        ## Here we have three files: data.tsv, genome.jp, nucleotides.txt
        ## data.tsv has the following columns:
            ## Species, identifier, locus tag, gene name, description, module id, chromosome, start, end, strand
        ## We first retrieve the modules present
        modules = set()
        with open("data.tsv", "r") as f:
            for line in f:
                tmp = line.split("\t")[5].split("_")[1].rstrip()
                modules.add(tmp)
        modules = list(modules)
        modules.sort() # sort for reproducibility
        nmodules = len(modules)
        ## Now compute information at module-level
        m = [""] * len(modules)
        with open("data.tsv", "r") as f:
            for mod in modules:
                starts = []
                ends = []
                strands = []
                f.seek(0) # rewind file pointer
                
                if module_types[identifiers.index(mod)]!="Pathway module":
                    continue
                m[modules.index(mod)] = mod
                for line in f:
                    start = -1
                    end = -1
                    strand = '.'
                    if line.find(mod) > -1:
                        ## Unfortunately the information isnt' always reported in the API.
                        ## We therefore need to pass through some try except statements
                        try:
                            start = int(line.split("\t")[7].replace(">", "").replace("<", ""))
                        except:
                            continue
                        try:
                            end = int(line.split("\t")[8].replace(">", "").replace("<", ""))
                        except:
                            continue                
                        try:
                            strand = line.split("\t")[9].replace(">", "").replace("<", "").rstrip()
                        except:
                            continue
                        if start != -1 and end != -1 and strand != '.':
                            starts.append(start)
                            ends.append(end)
                            strands.append(strand)
                ## check the number of genes were found for each module
                ## we want at least 3 genes
                
                if len(starts) < 3:
                    continue
                
                ## if this point is reached, it means that there were at least 3 genes found
                ## now we need to take the chromosome size into account
                chr_size = 0
                with open("nucleotides.txt", "r") as n:
                    chr_size = int(n.readline().rstrip())
                ## Now invoke a function that gathers all the information 
                proximity = calcProximity4(starts, ends, strands, chr_size)
                for element in proximity:
                    ## You need to comment/uncomment depending on the desired normalization
                    #modules_average.append(element)
                    modules_average.append(math.log2(element))
                
        os.chdir(species_path)
    if len(modules_average) < 1:
        continue

    print(f'{species.replace(" ", "_")},{num_genomes},{np.percentile(modules_average, 25)},{np.percentile(modules_average, 75)}, {np.mean(modules_average)}, {np.median(modules_average)}')
    os.chdir(working_dir)