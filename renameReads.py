import sys
import sys
import re
import os
from Bio import SeqIO
from bisect import bisect
import operator


class SAMEntry:
    def __init__(self, props):
        self.qname = props[0]
        self.flag = props[1]
        self.flagBitRev = str(
            "{000000000000000:b}".format(int(self.flag)))[::-1]
        self.flagBitRev = self.flagBitRev + '0'*(15 - len(self.flagBitRev))
        self.firstInPair = 0
        if (self.flagBitRev[6] == '1'):
            self.firstInPair = 1
        self.rcmp = 0
        if (self.flagBitRev[4] == '1'):
            self.rcmp = 1
        self.rname = props[2]
        self.pos = int(props[3]) - 1  # change to 0-based
        self.mapq = props[4]
        self.rnext = props[6]
        self.pnext = props[7]
        self.tlen = props[8]
        self.seq = props[9]
        self.qual = props[10]
        self.optional = {}
        for i in range(11, len(props)):
            field = props[i].split(':')
            if field[1] == 'i':
                field[2] = int(field[2])
            elif field[1] == 'f':
                field[2] = float(field[2])
            self.optional[field[0]] = [field[1], field[2]]
        self.len = len(self.seq)

    # how to interpret flag in sam file
    # 00000000000
    # However SAM uses this single number as a series of boolean (true false) flags where each position in the
    # array of bits represents a different sequence attribute.
    # Bit 0 = This query (read) was part of a pair during sequencing
    # Bit 1 = Both this query (read) and the mate of the pair are mapped
    # Bit 2 = The query sequence is unmapped
    # Bit 3 = The mate is unmapped
    # Bit 4 = Strand of query (0=forward 1=reverse)
    # Bit 5 = Strand of mate (0=forward 1=reverse)
    # Bit 6 = The query is the first sequence of a pair
    # Bit 7 = The query is the second sequence of a pair
    # Bit 8 = The alignment reported is not the best alignment of this sequence
    # Bit 9 = The sequence or alignment does not pass quality controls
    # Bit 10 = The sequence is considered to be a PCR or optical duplicate


def sam_parse(ifname):
    infile = open(ifname)
    sam = []
    i = 0
    for line in infile:
        if line.startswith('@'):  # header
            continue
        else:  # entry
            props = line.strip().split('\t')
            if props[5] != '*':
                sam += [SAMEntry(props)]
                print(sam[i].flagBitRev)
        i = i + 1
    infile.close()
    return sam


def RC(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return "".join([seq_dict[base] for base in reversed(seq)])


def findCorrectReads(Reads, sam_file, outCorrectedReads):
    from Bio import SeqIO
    infile = open(sam_file, "r")
    oufFastq = open(outCorrectedReads, 'w')
    first = 0
    with open(Reads, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if (first == 1):
                first = 0
            else:
                first = 1

            samLine = next(infile)
            while(samLine.startswith('@')):  # header
                samLine = next(infile)
            props = samLine.strip().split('\t')
            samObj = SAMEntry(props)
            while ((samObj.firstInPair == 0 and first == 1)):
                samLine = next(infile)
                props = samLine.strip().split('\t')
                samObj = SAMEntry(props)
            while ((samObj.firstInPair == 1 and first == 0)):
                samLine = next(infile)
                props = samLine.strip().split('\t')
                samObj = SAMEntry(props)

            readID = record.id.split("/")[0]
            if (first == 1):
                readID = readID + ".1"
            else:
                readID = readID + ".2"
            seq = str(record.seq)
            if (samObj.rcmp):
                readID = str(readID) + "\tRC"
                seq = RC(seq)
            oufFastq.write("@"+readID+"\n")
            oufFastq.write(seq+"\n")
            oufFastq.write("+\n")
            oufFastq.write(str(record.format("fastq")).split("\n")[3]+"\n")
    oufFastq.close()
    infile.close()


def main(argv=None):
    if argv == None:
        argv = sys.argv
    if len(argv) < 3:
        print('Usage: renameReads.py reads.fastq samFile.sam  output.fastq')
        exit()
    reads_input = argv[1]
    sam_file = argv[2]
    outCorrectedReads = argv[3]
    findCorrectReads(reads_input, sam_file, outCorrectedReads)


main()
