# Transform files into formats that other programs need

import string
import re

def format_1(ifile,ofile):

    f = open(ifile,'r')
    of = open(ofile,'w')
    for line in f:
        if '>' in line:
            fields = line.split()
            info = string.join(fields[:2],' ')
            of.write(info+'\n')
            of.write(fields[-1]+'\n')
        else:
            of.write(line)
    f.close()
    of.close()

def BreakFusion_clean(ifile,ofile):
    f = open(ifile,'r')
    of = open(ofile,'w')
    for line in f:
        columns = line.split('\t')
        if re.match('Gene',columns[14]):
            of.write(line)
        else:
            newline = '\t'.join(columns[:13])+'\t'+' '+'\t'+'\t'.join(columns[13:])
            of.write(newline)


def ChimeraScan_clean(ifile,ofile):
    f = open(ifile,'r')
    of = open(ofile,'w')
    line_1 = f.readline()
    if line_1[0] == '#':
        pass
    else:
        of.write(line_1)
    for line in f:
        of.write(line)
    


def overlapRatio(start,end,txStart,txEnd):
    if end<txStart or start>txEnd:
        return 0.0
    elif start<txStart:
        return float(end-txStart+1)/float(end-start+1)
    elif end>txEnd:
        return float(txEnd-start+1)/float(end-start+1)
    else:
        return 1.0

def mapGeneOnChrom(chrom,start,end,allGenes):
    score_high = 0
    mapGene = 'null'
    try:
        allGenesOnChrom = allGenes[chrom]
    except:
        # Chromosome not listed
        return ''
    for gene in allGenesOnChrom:
        txStart = allGenesOnChrom[gene]['start']
        txEnd = allGenesOnChrom[gene]['end']
        score = overlapRatio(start,end,txStart,txEnd)
        if score > 0.5:
            return gene
        elif score > score_high:
            score_high = score
            mapGene = gene
        else:
            pass
    return mapGene



def GeneAnnot(ifile,ofile,RefSeqFile):
    allGenes = {}
    f = open(RefSeqFile,'r')
    f.readline()
    chromList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
                 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
                 'chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    for c in chromList:
        allGenes[c] = {}

    for line in f:
        chrom = line.split()[2]
        name2 = line.split()[12]
        start = int(line.split()[4])
        end = int(line.split()[5])
        for c in chromList:
            if chrom == c:
                allGenesOnChrom = allGenes[chrom]
                if name2 not in allGenesOnChrom.keys():
                    allGenesOnChrom[name2] = {'start':start,'end':end}

    f.close()

    # Scan integrated gene fusion detection result for genes and replace with
    # official gene names
    f = open(ifile,'r')
    of = open(ofile,'w')
    of.write(f.readline())

    for line in f:
        cols = line.split('\t')
        chrom = cols[3]
        # deal with irregular ifile
        if chrom[0] == '#':
            continue
        start = int(cols[5])
        end = int(cols[6])
        mapGene = mapGeneOnChrom(chrom,start,end,allGenes)
        if mapGene:
            cols[4] = mapGene
        chrom = cols[9]
        start = int(cols[11])
        end = int(cols[12])
        mapGene = mapGeneOnChrom(chrom,start,end,allGenes)
        if mapGene:
            cols[10] = mapGene
        newline = '\t'.join(cols)
        of.write(newline)

    t = time()
    f.close()
    of.close()
