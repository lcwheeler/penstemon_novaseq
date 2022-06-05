stat_tree = 'no_namemap.cf.stat_tree' #stat_tree output from iqtree --gcf --cf-verbose
treepaths = 'numbered_treepaths.txt' #key with numbered treepaths
nodenumber = 29 #node id as labeled by gcf
disc_tree = 1 #this is either 1 or 2: (gD1 or gD2)
outfile_name = 'discordant_trees.txt'


with open(stat_tree, 'r') as infile1:
    with open(treepaths, 'r') as infile2:
        incols = infile1.read().splitlines()
        intrees = infile2.read().splitlines()
        with open(outfile_name, 'w') as output:
            for i in range(0,len(incols)):
                if incols[i].split()[0] == str(nodenumber) and incols[i].split()[disc_tree+2] == str(1):
                    treenumber = int(incols[i].split()[1])
                    output.write(intrees[treenumber-1] + '\n')

