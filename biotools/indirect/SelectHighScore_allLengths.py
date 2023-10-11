import sys, re
#sys.path.append("c:\\usr\\biotools\\system\\PYTHON")
from Bio import SeqIO

# Authors: Karin Lagesen
# For license see /usr/biotools/CMG-biotools.license

def getFirstOK(infilename):
    printrec = "NoneFound"
    handle = open(infilename, "rU")
    maxlength = 0
    for seq_record in SeqIO.parse(handle, "fasta"):
        if printrec == "NoneFound":
            if len(seq_record.seq) > maxlength:
                maxlength = len(seq_record.seq)
                printrec = seq_record
    handle.close()
    return(printrec,maxlength)

def changeName(seq_record, infilename, uniqid):
    orgnameA = re.sub(".rrna", "", infilename)
    orgname = re.sub("^[A-Z]([a-z]+)_", orgnameA[0]+ "_", orgnameA)
    seq_record.id = uniqid + "_" + orgname
    return seq_record

if __name__ == "__main__":
    infilename = sys.argv[1]
    uniqid = sys.argv[2]
    chosenseq, maxlength = getFirstOK(infilename)
    if chosenseq != "NoneFound":
        newnameseq = changeName(chosenseq, infilename, uniqid)
        print newnameseq.format("fasta")
    else:
        print >> sys.stderr, infilename + " has a max length of " + str(maxlength) + ", and is thus not included"
