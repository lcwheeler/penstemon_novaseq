import pandas as pd
import argparse


# Define a function that replaces tip labels
def replace_tip_labels(keyfile, treefile):
     
    """
    Replace the tip labels of a Newick tree with new labels using a dictionary of
    keys stored in a csv file.
    """
    
    tk = pd.read_csv(keyfile, sep=",", index_col=0, header=None).to_dict()[1]

    with open(treefile, "r") as fin:
        ts = str(fin.readline())
        for key in tk.keys():
            if str(key) in ts:
                #print(str(key))
                #print(str(tk[key]))
                ts = ts.replace(str(key), str(tk[key]))
            else:
                pass
        
    return(ts)
    
    
# Add argparse main function to take command line arguments (input and output file names)
def main():

    """
        Run the replace_tip_labels function on input files using argparse arguments.
    """    
    
    # Set up for using argparse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-k", "--keyfile", help="name of the input CSV key file")
    parser.add_argument("-t", "--treefile", help="name of the input Newick tree file")
    parser.add_argument("-o", "--outfile", help="name of the output relabeled tree file")

    args = parser.parse_args()


    # Run the function to replace names
    ts = replace_tip_labels(args.keyfile, args.treefile)
    ts = ts.replace("'", "")

    # Save the tree with relabeled tips as a new file
    with open(args.outfile, "w") as fout:
        fout.write(ts)
        
if __name__ == "__main__":
    main()

