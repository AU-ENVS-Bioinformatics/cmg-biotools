import sys, os.path
#sys.path.append("c:\\usr\\biotools\\system\\PYTHON")
# Authors: Karin Lagesen
# For license see /usr/biotools/CMG-biotools.license
inputfile = open(sys.argv[1], "r")
inlines = inputfile.readlines()
inputfile.close()

for line in inlines:
    source = line.strip("\n")
    name = os.path.basename(line)  
    name = name.strip("\n")
    name = name.replace(".proteins.fna", "")
    name = name.replace(".proteins.fsa", "")
    name = name.replace(".orf.fna", "")
    name = name.replace(".orf.fsa", "")
    name = name.replace(".fna.BSU03.orfs.fsa","")
    name = name.replace("_supercontigs.fasta.BSU03.orfs.fsa", "")
    name = name.replace(".fasta.BSU03.orfs.fsa", "")
    name = name.replace(".refseq.proteins.fsa", "")
    name = name.replace("../data/", "")

    print "<entry>\n\t<length>unused</length>\n\t<source>" + source + "</source>\n\t<title>" + name + "</title>\n\t<subtitle></subtitle>\n</entry>"

print "</sources>\n</compare>"
