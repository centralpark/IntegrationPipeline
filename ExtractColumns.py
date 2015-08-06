import string

output_col = open('/Users/HSH/Desktop/Output_ColNames.txt','w')

f = open('/Users/HSH/Desktop/output.txt','r')

line = f.readline()

while line:
    words = line.split(':')
    word = words[0]
    output_col.write(word.rstrip()+'\n')
    line = f.readline()

output_col.close()
f.close()
