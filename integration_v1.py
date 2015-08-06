import string

def inputFile(filename,tool):
    
    f = open(filename,'r')
    output = open('Rearrangement','a')
    header = f.readline()
    names = header.split('\t')
    if names[0]=='Sample ID from TCGA':
        col_5pChr = 1

    line = f.readline()
    if tool=='ChimeraScan':
        while line:
            words = line.split('\t')
            five_Chr = words[0].split(':')[1]
            five_start = words[1]
            five_end = words[2]
            five_gene = words[12]
            five_ori = words[8]
            three_Chr = words[3]
            three_start = words[4]
            three_end = words[5]
            three_gene = words[13]
            three_ori = words[9]
            fusion_type = words[14]
            spanning_reads = words[21]
            score = words[7]
            line_out = string.join([five_Chr,five_start,five_end,five_gene,
                                    five_ori,three_Chr,three_start,three_end,
                                    three_gene,three_ori,fusion_type],'\t')
            output.write(line_out+'\n')
            line = f.readline()
    elif tool=='BEDA':
        while line:
            words = line.split('\t')
            five_Chr = 'N/A'
            five_start = 'N/A'
            five_end = 'N/A'
            five_gene = words[1]
            five_ori = 'N/A'
            three_Chr = 'N/A'
            three_start = 'N/A'
            three_end = 'N/A'
            three_gene = words[2]
            three_ori = 'N/A'
            fusion_type = words[7]
            spanning_reads = words[11]
            
            line_out = string.join([five_Chr,five_start,five_end,five_gene,
                                    five_ori,three_Chr,three_start,three_end,
                                    three_gene,three_ori,fusion_type],'\t')
            output.write(line_out+'\n')
            line = f.readline()
    elif tool=='BreakFusion':
        while line:
            words = line.split('\t')
            five_Chr = words[0].split(':')[1]
            five_start = words[1]
            five_end = 'N/A'
            five_gene = 'N/A'
            five_ori = words[2]
            three_Chr = words[3]
            three_start = words[4]
            three_end = 'N/A'
            three_gene = 'N/A'
            three_ori = words[5]
            fusion_type = words[6]
            spanning_reads = 'N/A'
            score = words[8]
            line_out = string.join([five_Chr,five_start,five_end,five_gene,
                                    five_ori,three_Chr,three_start,three_end,
                                    three_gene,three_ori,fusion_type],'\t')
            output.write(line_out+'\n')
            line = f.readline()
    else:
        print 'The specified gene fusion tool is not supported!'
    
    f.close()
    output.close()


