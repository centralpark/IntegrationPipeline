import urllib
import string
import re
from bs4 import BeautifulSoup

def FindKinaseDomain(Entry):
    
    url = string.join(['http://www.uniprot.org/uniprot/',Entry,'.txt'],'')

    f = urllib.urlopen(url)

    load_file = open(Entry,'w')
    load_file.write(f.read())
    load_file.close()
    f.close()

    info = open(Entry,'r')
    text = info.read()
    loc = text.find('FT   DOMAIN')
    if loc>=0:    
        info.seek(loc)

        FT_end = True
        found = False
        while FT_end:
            line = info.readline()
            params = line.split()
            Description = string.join(params[4:],' ')
            if not (params[0]=='FT' and params[1]=='DOMAIN'):
                break
            if 'kinase' in Description:
                print 'kinase domain found!'
                print 'Description:',Description
                print 'Position(s):',params[2],'-',params[3]
                found = True

        if not found:
            print 'kinase domain not found!'

    else:
        print 'kinase domain not found!'


def UniprotInfo(GeneName):
    
    urls = 'http://www.uniprot.org/uniprot/?query='+GeneName+\
      '+AND+organism:Human+reviewed:yes&sort=score'
    query = urllib.urlopen(urls)
    info = {}
    
    try:
        results_page = BeautifulSoup(query.read())
        results = results_page.find('table',attrs={'id':'results'})
        result = results.tr.next_sibling
        entry = result.find('a',href=re.compile('/uniprot/')).string.encode('ascii')

        query.close()
        
        info['pEntry'] = entry

        url = 'http://www.uniprot.org/uniprot/'+entry
        html = urllib.urlopen(url)
        page = BeautifulSoup(html.read())
        html.close()

        comments = page.find('div', attrs={'id':'content-comments'})
        row_locs = comments.find_all('a',href=re.compile('/locations'))
        temp = []
        for loc in row_locs:
                temp.append(loc.string.encode('ascii'))
        info['pLoc'] = '.'.join(temp)
        ontologies = page.find('div',attrs={'id':'content-terms'})
        for row in ontologies.table.tbody.children:
                if not row.find('td'):
                        continue
                if re.search('Molecular function',row.td.string.encode('ascii','ignore')):
                        temp = []
                        for func in row.find_all('a',href=True):
                                temp.append(func.string.encode('ascii'))
                        info['pFunction'] = '.'.join(temp)
                        break

        
        
        return info
    except:
        info = {}
        return info

def integrateFile(filename,tool,ofile):
    # define header columns
    Tool = 0
    sampleID = 1
    UNC_ID = 2
    chrom1 = 3
    gene1 = 4
    gene1_start = 5
    gene1_end = 6
    exon1 = 7
    gene1_ori = 8
    chrom2 = 9
    gene2 = 10
    gene2_start = 11
    gene2_end = 12
    exon2 = 13
    gene2_ori = 14
    Type = 15
    total_fragments = 16
    spanning_fragments = 17
    chimera_cluster_id = 18
    score = 19
    transcript_id_5 = 20
    transcript_id_3 = 21

    N_col = 22


    f = open(filename,'r')
    output = open(ofile,'a')

    if tool=='ChimeraScan':
        for line in f:
            line_out = [' ']*N_col
            columns = line.split('\t')
            line_out[Tool] = 'ChimeraScan'
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
            if columns[16] =='Read_Through':
                line_out[Type] = 'Read-through'
            line_out[total_fragments] = columns[18]
            line_out[spanning_fragments] = columns[19]
            line_out[chimera_cluster_id] = columns[8]
            line_out[score] = columns[9]
            line_out[transcript_id_5] = columns[12]
            line_out[transcript_id_3] = columns[13]
            output.write('\t'.join(line_out)+'\n')

    elif tool=='BEDA':
        for line in f:
            line_out = [' ']*N_col
            columns = line.split('\t')
            line_out[Tool] = 'BEDA'
            line_out[UNC_ID] = columns[0]
            line_out[gene1] = columns[1]
            line_out[gene2] = columns[2]
            line_out[Type] = columns[7]
            if line_out[Type] == 'intra':
                line_out[Type] = 'Intrachromosomal'
            elif line_out[Type] =='inter':
                line_out[Type] = 'Interchromosomal'
            else:
                pass
            output.write('\t'.join(line_out)+'\n')

    elif tool=='BreakFusion':
        for line in f:
            line_out = [' ']*N_col
            columns = line.split('\t')
            line_out[Tool] = 'BreakFusion'
            line_out[sampleID] = columns[0]
            line_out[chrom1] = columns[2]
            line_out[gene1_start] = columns[3]
            line_out[chrom2] = columns[5]
            line_out[gene2_start] = columns[6]

            matchObj = re.match('Gene:(.*?)\|(.*?),.*?\|.*?:(.*?)-.*?\|.*?:(.*?),(.*)',columns[14])
            line_out[gene1] = matchObj.group(1)
            line_out[gene2] = matchObj.group(2)
            if matchObj.group(3) == 'NA':
                line_out[exon1] = 'NA'
                line_out[gene1_start] = 'NA'
                line_out[gene1_end] = 'NA'
            else:
                temp = re.split(':|/',matchObj.group(3))
                line_out[exon1] = temp[0]
                line_out[gene1_start] = temp[1]
                line_out[gene1_end] = temp[2]
            if matchObj.group(4) == 'NA':
                line_out[exon2] = 'NA'
                line_out[gene2_start] = 'NA'
                line_out[gene2_end] = 'NA'
            else:
                temp = re.split(':|/',matchObj.group(4))
                line_out[exon2] = temp[0]
                line_out[gene2_start] = temp[1]
                line_out[gene2_end] = temp[2]
            temp = re.match('Fusion,(.*)',matchObj.group(5))
            if temp:
                line_out[Type] = temp.group(1)
                if line_out[Type] == 'Readthrough':
                    line_out[Type] = 'Read-through'
            output.write('\t'.join(line_out)+'\n')
    else:
        print 'Gene fusion tool not supported!'

    f.close()
    output.close()



def integrate():
    ofile = '/Users/HSH/Desktop/Rearrangement'
    f = open(ofile,'w')
    f.write('Gene Fusion Tool\t'
                'Sample ID\t'
                'UNC_ID\t'
                'Chromosome 1\t'
                'Gene 1\t'
                'Gene 1 Start\t'
                'Gene 1 End\t'
                'Exon 1\t'
                'Gene 1 Orientation\t'
                'Chromosome 2\t'
                'Gene 2\t'
                'Gene 2 Start\t'
                'Gene 2 End\t'
                'Exon 2\t'
                'Gene 2 Orientation\t'
                'Fusion Type\t'
                'Total Fragments\t'
                'Spanning Fragments\t'
                'Chimera Cluster ID\t'
                'Score\t'
                '5\' transcript ids\t'
                '3\' transcript ids\n')

    GeneFusion.integrateFile('/Users/HSH/Desktop/Data/all_chimerascan.txt',\
                             'ChimeraScan',ofile)
    GeneFusion.integrateFile('/Users/HSH/Desktop/Data/beda_fusions.txt',\
                             'BEDA',ofile)
    GeneFusion.integrateFile('/Users/HSH/Desktop/Data/all_breakfuse.tsv',\
                             'BreakFusion',ofile)
