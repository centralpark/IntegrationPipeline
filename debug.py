import re
from Bio import Entrez
from StringIO import StringIO
gene = 'CSDE1'


Entrez.email = 'seas363@gmail.com'
query = '('+gene+'[Gene Name]) AND human[Organism]'
handle = Entrez.esearch(db = 'gene',term = query)
record = Entrez.read(handle)
handle.close()
ID = record['IdList'][0]
handle = Entrez.efetch(db = 'gene', id = ID, rettype = 'gene_table', retmode = 'text')
record = handle.read()
handle.close()

ind = record.find('annotated AA length:')
aa_len = 0
index = 0
while ind > 0:
    matchObj = re.search('annotated AA length: ([0-9]*)',record[ind:ind+30])
    temp = matchObj.group(1)
    if temp > aa_len:
        aa_len = temp
        index = ind
    ind = record.find('annotated AA length:',ind+30)

f = StringIO(record)
f.seek(index)
for i in range(5):
    f.readline()

line = f.readline()
while line and line!='\n':
    print line
    [genom_inter_exon,genom_inter_code,gene_inter_exon,gene_inter_code,exon_len,code_len,intron_len] = line.split('\t')
    print gene_inter_code
    line = f.readline()
    
f.close()
