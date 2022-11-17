### Author: Lucas C. Wheeler ### 

# Import dependencies
import argparse
import re

def convert_tree(treefile, pattern=r'\[.*?\]', support_name="pp1=", outfile="reformatted-tree.tre"):

    """This function converts a newick tree with bracketed lists of support values output from Astral 
    into a file with a single selected support value that can be opened in figtree. 
    """

    
    #ts = None
    #with open(treefile, "r") as fin:
        #for line in fin:
            #ts = line
            
    fin = open(treefile, "r")
    ts = str(fin.readline())
    ts = ts.replace("'", "")
    
    # Generate a list of all the hits to the support value list pattern
    mlist = re.findall(pattern, ts)
    print(mlist)
            
    # Replace each list of supports with the single pp1 value
    for m in mlist:
        rs = m.split(support_name)[1].split(";")[0]
        rs = str(float(rs))
        ts = ts.replace(m, rs)
        
    ts = ts.replace("'", "")
    
    
    # Save the reformatted tree file 
    
    with open(outfile, "w") as fout:
        fout.write(ts)
        
        
# Add argparse main function to take command line arguments (input and output file names)
def main():

    """
        Run the filter_fastas function on input files using argparse arguments.
    """    
    
    
    # Set up for using argparse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-treefile", help="Newick tree file from Astral")
    parser.add_argument("-pattern", help="Regex pattern to replace in tree string")
    parser.add_argument("-support_name", help="The selected support value to replace supports list")
    parser.add_argument("-outfile", help="Name of output reformatted tree file")

    args = parser.parse_args()


    # Run the function to convert the tree format
    convert_tree(args.treefile, args.pattern, args.support_name, args.outfile)

        
if __name__ == "__main__":
    main()

