from time import time
import pickle
import sys
sys.path.append('/Users/HSH/Documents')
import GeneFusion

t0 = time()

ifile = '/Users/HSH/Desktop/tmp_fusions_ov_chimerascan_integrate.txt'
ofile = '/Users/HSH/Desktop/tmp_fusions_ov_chimerascan_annotate.txt'

f = open(ifile,'r')
of = open(ofile,'w')
header = f.readline()
colName = header.split('\t')
colName_new = colName[:9] + \
              ['Protein 1 Uniprot Entry','Protein 1 Subcellular Location','Protein 1 Function'] + \
              colName[9:15] + \
              ['Protein 2 Uniprot Entry','Protein 2 Subcellular Location','Protein 2 Function'] + \
              colName[15:]
header_new = '\t'.join(colName_new)
of.write(header_new)

N = 0

# import saved Uniprot information database
pkl_file = open('/Users/HSH/Documents/protInfoDict.pkl','rb')
protInfoDict = pickle.load(pkl_file)

for line in f:
    col = line.split('\t')
    tool = col[0]
    if tool == 'ChimeraScan':
        gene1 = col[4].split(',')[0]
        gene2 = col[10].split(',')[0]
    else:
        gene1 = col[4]
        gene2 = col[10]
    # If the gene is searched before, no need to search again
    # saves run time on large input file
    b1 = (gene1 in protInfoDict)
    b2 = (gene2 in protInfoDict)
    if b1:
        protInfo1 = protInfoDict[gene1]
    else:
        protInfo1 = GeneFusion.UniprotInfo(gene1)
    if b2:
        protInfo2 = protInfoDict[gene2]
    else:
        protInfo2 = GeneFusion.UniprotInfo(gene2)
    
    # Some protein information does not exist
    p_keys = ['pEntry','pLoc','pFunction']
    p_value1 = []
    p_value2 = []
    for p_key in p_keys:
        try:
            p_value1.append(protInfo1[p_key])
        except KeyError:
            p_value1.append('')
        try:
            p_value2.append(protInfo2[p_key])
        except KeyError:
            p_value2.append('')
    if not b1:
        protInfoDict[gene1] = {'pEntry':p_value1[0],'pLoc':p_value1[1],'pFunction':p_value1[2]}
    if not b2:
        protInfoDict[gene2] = {'pEntry':p_value2[0],'pLoc':p_value2[1],'pFunction':p_value2[2]}
    
    col_new = col[:9] + p_value1 + col[9:15] + p_value2 + col[15:]
    line_new = '\t'.join(col_new)
    of.write(line_new)

    N+=1
    if N > 10:
        break
    
f.close()
of.close()

t = time()
