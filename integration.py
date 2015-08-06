import string

def inputFile(filename,tool):
    # define header columns
    sampleID = 0
    UNC_ID = 1
    chrom1 = 2
    gene1 = 3
    gene1_start = 4
    gene1_end = 5
    gene1_ori = 6
    chrom2 = 7
    gene2 = 8
    gene2_start = 9
    gene2_end = 10
    gene2_ori = 11
    Type = 12
    N_col = 13
    
    f = open(filename,'r')
    output = open('Rearrangement','a')
    
    if tool=='ChimeraScan':
        for line in f:
            line_out = ['N\A']*N_col
            columns = line.split('\t')
            line_out[sampleID] = columns[0]
            line_out[UNC_ID] = columns[1]
            line_out[chrom1] = columns[2]
            line_out[gene1] = columns[14]
            line_out[gene1_start] = columns[3]
            line_out[gene1_end] = columns[4]
            line_out[gene1_ori] = columns[10]
            line_out[chrom2] = columns[5]
            line_out[gene2] = columns[15]
            line_out[gene2_start] = columns[6]
            line_out[gene2_end] = columns[7]
            line_out[gene2_ori] = columns[11]
            line_out[Type] = columns[16]
            output.write(string.join(line_out,'\t')+'\n')

    elif tool=='BEDA':
        for line in f:
            line_out = ['N\A']*N_col
            columns = line.split('\t')
            line_out[UNC_ID] = columns[0]
            line_out[gene1] = columns[1]
            line_out[gene2] = columns[2]
            line_out[Type] = columns[7]
            output.write(string.join(line_out,'\t')+'\n')

    elif tool=='BreakFusion':
        for line in f:
            line_out = ['N\A']*N_col
            columns = line.split('\t')
            line_out[sampleID] = columns[0]
            line_out[chrom1] = columns[2]
            line_out[gene1_start] = columns[3]
            line_out[chrom2] = columns[5]
            line_out[gene2_start] = columns[6]
            line_out[Type] = columns[8]
            output.write(string.join(line_out,'\t')+'\n')

    else:
        print 'Gene fusion tool not supported!'
    
    f.close()
    output.close()
    
