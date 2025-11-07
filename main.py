import argparse
import os
import glob
from multiprocessing import Process, Manager, Lock
import pandas as pd
import numpy as np
import tskit
from best_subtree_finder import best_subtree_finder
from node_recalculate import node_recalculate
import re

def atoi(text):
    return int(text) if text.isdigit() else text
    
def natural_keys(text): #human sorting by increasing number in text
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]
    
def main_run(lock, ts_path, chromosome, phenotypes, nIndiv, allpositions, allr2, allslopes, mapfile, pedfile):
    ts = tskit.load(ts_path)
    nbtree = ts.get_num_trees()
    best_r2 = 0 
    for tree_number in range(nbtree):
        tree = ts.at_index(tree_number)
        if tree_number == 0:
            nodelist = list(range(nIndiv * 2))
        else:
            tree_prev = ts.at_index(tree_number - 1)
            nodelist = node_recalculate(tree_prev, tree)
        
        result = best_subtree_finder(tree, nodelist, nIndiv * 2, nIndiv * 2 * 0.05, phenotypes)
        if result[2]:
            allpositions[chromosome-1].append(int(tree.interval.mid))
            allr2[chromosome-1].append(result[0])
            allslopes[chromosome-1].append(result[1])
            with lock:
                # pseudo-map file
                with open(mapfile, "a") as file1:  
                    file1.write(str(chromosome)+" ARG"+str(tree_number)+"_"+str(chromosome)+" 0 "+str(int(tree.interval.mid))+"\n")
                # pseudo-ped file
                with open(pedfile, "a") as file1:  
                    file1.write("ARG"+str(tree_number)+"_"+str(chromosome)+" "+result[2]+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ARGWAS",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--output", help="output directory for the results")
    parser.add_argument("-p", "--ped", help="ped file")
    parser.add_argument("-t", "--trees", help="trees files folder")
    args = parser.parse_args()
    
    output_dir = args.output
    ped_file = args.ped
    trees_dir = args.trees
    
    os.makedirs(output_dir, exist_ok=True)
    lock = Lock()

    ped_table = pd.read_csv(ped_file, header=None, sep="\t")
    treefiles = glob.glob(f"{trees_dir}/*.trees")
    treefiles.sort(key=natural_keys) # sort by increasing number 
    nchr = len(treefiles)

    phenotypes = ped_table[5] - np.mean(ped_table[5])
    nIndiv = len(ped_table)

    with Manager() as manager:
        allpositions = manager.list([manager.list() for _ in range(nchr)])
        allr2 = manager.list([manager.list() for _ in range(nchr)])
        allslopes = manager.list([manager.list() for _ in range(nchr)])
        
        processes = []
        for chromosome in range(1, nchr + 1):
            mapfile=f"{output_dir}/results.map"
            pedfile=f"{output_dir}/results.ped"
            ts_path = treefiles[chromosome-1]
            p = Process(target=main_run, args=(lock, ts_path, chromosome, phenotypes, nIndiv, allpositions, allr2, allslopes,mapfile,pedfile))
            processes.append(p)
            p.start()
        
        for p in processes:
            p.join()

        # Save results
        data = []
        for chrom_idx, (positions, r2s, slopes) in enumerate(zip(allpositions, allr2, allslopes), start=1):
            for pos, r2, slo in zip(positions, r2s, slopes):
                data.append((chrom_idx, pos, r2, slo))
        results = pd.DataFrame(data, columns=["chromosome", "position", "r2", "slope"])

        results.to_csv(f"{output_dir}/results.csv", index=False)
        
        # Transform pseudo-map and pseudo-ped into real map and ped files
        # Putative QTLs ancestral state is "A" and new state is "T"
        Mapfile = pd.read_csv(mapfile, sep=" ", header=None)
        Pedfile = pd.read_csv(pedfile, sep=" ", header=None)
        
        # Sort Pedfile and Mapfile by chromosome then position
        Mapfile = Mapfile.sort_values([0, 3]).reset_index(drop=True)
        Pedfile = Pedfile.iloc[Mapfile.index].reset_index(drop=True)

        # First 6 columns of ped_table
        newPed = ped_table.iloc[:, :6].copy()

        # Transpose Pedfile and replace 0->A, 1->T
        addPed = Pedfile.iloc[:, 1:].T.replace({0: "A", 1: "T"}).to_numpy()

        # Reshape trick: split even/odd rows without a Python loop
        odd_rows = addPed[0::2, :]   # rows 0,2,4... in numpy
        even_rows = addPed[1::2, :]  # rows 1,3,5...

        # Interleave odd and even rows into new columns
        geno_matrix = np.empty((odd_rows.shape[0], odd_rows.shape[1] * 2), dtype=object)
        geno_matrix[:, 0::2] = odd_rows
        geno_matrix[:, 1::2] = even_rows

        # Convert to DataFrame
        geno_df = pd.DataFrame(geno_matrix)

        # Concatenate once
        newPed = pd.concat([newPed.reset_index(drop=True), geno_df], axis=1)
        
        # Save files
        Mapfile.to_csv(mapfile, sep="\t", header=False, index=False)
        newPed.to_csv(pedfile, sep="\t", header=False, index=False)
