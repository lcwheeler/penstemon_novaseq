### Author: Lucas C. Wheeler ### 

# Import dependencies
from Bio import SeqIO
import glob


# Specificy the pattern to recognize input fasta files
fpattern = "CDS_full*.fa"

# Specify the reference fasta file to pull the seq IDs
reffile = "CDS_full_cardBL_1_S443.fa"

scaf_dict = {}

with open(reffile, "r") as ref:

    # Iterate over all the individual sequences in the file
    for seq_record in SeqIO.parse(reffile, "fasta"):
        scaf_dict[str(seq_record.id)] = str(seq_record.id) + ".rearranged.fa"

# Keep a list of recorded IDs to prevent saving duplicates
recorded = []
        
# Iterate over all files with the defined filename pattern
for file in glob.glob(fpattern):

    sample_name = file.replace(".fa", "")

    with open(file, "r") as con_fa:

        # Iterate over all the individual sequences in the file
        for seq_record in SeqIO.parse(con_fa, "fasta"):
            
            name_tracker = str(seq_record.id) + "_" + sample_name

            if name_tracker not in recorded:
                recorded.append(name_tracker)

                fout = open(scaf_dict[str(seq_record.id)], "a")
                fout.write(">" + str(seq_record.id) + "_" + sample_name + "\n") 
                fout.write(str(seq_record.seq) + "\n")
                fout.close()

            else:
                pass
