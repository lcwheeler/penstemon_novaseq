### Author: Lucas C. Wheeler ### 

# Import dependencies
from Bio import SeqIO
import glob
import numpy as np
import pandas as pd
import argparse


def filter_fastas(fpattern, rm_ext, save_ext, per_n, num_samples, id_split_string):

    """This function takes fasta files containing nucleotide sequences as input and 
    filters them based on the percentage of missing data. Sequences that don't pass
    the cutoff test (set as per_n) are removed from the fasta files. Data on the total
    number of removed sequences and the identities of removed sequences from each file
    are output as CSV files. Lists of files that pass and fail based on a cutoff (num_samples)
    of individual sequences that pass are also output."""
    
    # Define dictionaries to collect sequence stats
    rmlist = {} 
    ncount = {}
    ldict = {}
    pcount = {}
    
    per_n = float(per_n)
    num_samples = float(num_samples) 

    # Iterate over all files with the defined filename pattern
    for file in glob.glob(fpattern):

        # Open the fasta file to read with BioPython
        with open(file, "r") as fin:
        
            name = file.replace(rm_ext, "")    
            
            # Counter to keep track of the number of seqs per file that don't pass
            ncount[file] = 0
            ldict[file] = 0
            pcount[name+save_ext] = 0
            
            # Open a new file to store the seqs that pass the missing data test
            with open(name+save_ext, "w") as faout:
            
                rmlist[file] = {}

                # Iterate over all the individual sequences in the file
                for seq_record in SeqIO.parse(file, "fasta"):
                
                    sequence = str(seq_record.seq).upper()
                    sequence = sequence.replace("W", "N").replace("R", "N").replace("Y", "N").replace("S", "N").replace("K", "N").replace("M", "N")
                    ID = str(seq_record.id)
                    #print(ID.split(id_split_string))
                    ID = ID.split(id_split_string)[1]
                    
                    # Calculate the percentage of the seq that is missing data
                    pc = float(sequence.count("N")) / float(len(sequence)) * 100
                    
                    # Test whether the % missing data passes the cutoff test
                    if pc <= per_n:
                        faout.write(">"+ID+"\n")
                        faout.write(sequence+"\n") 
                        ldict[file] += float(len(sequence))    
                        pcount[name+save_ext] += 1                    
                    else:
                        rmlist[file][ID] = pc
                        ncount[file] += 1
                    
                    
    # Store the stats on removed seqs for the dataset                
    rmlist_df = pd.DataFrame.from_dict(rmlist)
    rmlist_df.to_csv("removed-seqs-stats.csv", sep=",")

    ncount_df = pd.DataFrame.from_dict(ncount, orient="index")
    ncount_df.columns = ["num_above_N_cutoff"]
    ncount_df.to_csv("missing-data-counts.csv", sep=",")

    l_df = pd.DataFrame.from_dict(ldict, orient="index")
    l_df.columns = ["total_bp"]
    l_df.to_csv("total-sequence-length-data.csv", sep=",")


    # Save a list of the files where more than the allowed % cutoff of seqs have passed
    
    with open("passed-files.txt", "w") as pout:
        for key in pcount.keys():
            if float(pcount[key]) >= num_samples:
                pout.write(key + "\n")
            
        
# Add argparse main function to take command line arguments (input and output file names)
def main():

    """
        Run the filter_fastas function on input files using argparse arguments.
    """    
    
    
    # Set up for using argparse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-per_n", help="Percentage missing data cutoff (below this seqs are kept)")
    parser.add_argument("-fpattern", help="File patttern to generate list of fasta input files")
    parser.add_argument("-rm_ext", help="File extension to remove from name")
    parser.add_argument("-save_ext", help="File extension to use for saving cleaned files")
    parser.add_argument("-num_samples", help="Number of samples/seqs in each rearranged consensus file or number of samples cutoff for failing a file")    
    parser.add_argument("-id_split_string", help="Str in seq ID on which to split out species label")    
    
    args = parser.parse_args()


    # Run the function to filter the fasta files
    filter_fastas(args.fpattern, args.rm_ext, args.save_ext, args.per_n, args.num_samples, args.id_split_string)

        
if __name__ == "__main__":
    main()

    
