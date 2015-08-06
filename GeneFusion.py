import urllib
import string
import re
import pickle
from os import listdir
from os.path import join
from bs4 import BeautifulSoup
import fileFormat
from Bio import Entrez


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
    
    info = {}
    # Dealing with null gene name annotation
    if GeneName.lower()=='null':
        return info
    
    urls = 'http://www.uniprot.org/uniprot/?query='+GeneName+\
      '+AND+organism:Human+reviewed:yes&sort=score'
    query = urllib.urlopen(urls)

    try:
        results_page = BeautifulSoup(query.read())
        results = results_page.find('table',attrs={'id':'results'})
        if not results:
            # The protein may not be reviewed
            query.close()
            urls = 'http://www.uniprot.org/uniprot/?query='+GeneName+\
                   '+AND+organism:Human&sort=score'
            query = urllib.urlopen(urls)
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
        url = 'http://www.uniprot.org/uniprot/'+entry
        html = urllib.urlopen(url)
        page = BeautifulSoup(html.read())
        html.close()

        try:
            comments = page.find('div', attrs={'id':'content-comments'})
            row_locs = comments.find_all('a',href=re.compile('/locations'))
            temp = []
            for loc in row_locs:
                    temp.append(loc.string.encode('ascii'))
            info['pLoc'] = '.'.join(temp)
        except:
            info['pLoc'] = ''

        try:
            info['pFunction'] = ''
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
        except:
            pass
        
        return info
    except:
        # No Uniprot entry found, probably non-coding gene
        info['pEntry'] = 'NA'
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
            line_out[sampleID] = ''
            line_out[UNC_ID] = ''
            line_out[chrom1] = columns[0]
            line_out[gene1] = columns[12]
            line_out[gene1_start] = columns[1]
            line_out[gene1_end] = columns[2]
            line_out[gene1_ori] = columns[8]
            line_out[chrom2] = columns[3]
            line_out[gene2] = columns[13]
            line_out[gene2_start] = columns[4]
            line_out[gene2_end] = columns[5]
            line_out[gene2_ori] = columns[9]
            line_out[Type] = columns[14]
            if columns[14] =='Read_Through':
                line_out[Type] = 'Read-through'
            line_out[total_fragments] = columns[16]
            line_out[spanning_fragments] = columns[17]
            line_out[chimera_cluster_id] = columns[6]
            line_out[score] = columns[7]
            line_out[transcript_id_5] = columns[10]
            line_out[transcript_id_3] = columns[11]
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



def integrate(directory,ofile):
    # integrate all gene fusion results file in the directory
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
    f.close()

    for ifile in listdir(directory):
        if re.search('chimera',ifile,re.I):
            integrateFile(join(directory,ifile),'ChimeraScan',ofile)
        elif re.search('beda',ifile,re.I):
            integrateFile(join(directory,ifile),'BEDA',ofile)
        elif re.search('break',ifile,re.I):
            filename = ifile+'_CLEAN.txt'
            fileFormat.BreakFusion_clean(join(directory,ifile),join(directory,filename))
            integrateFile(join(directory,filename),'BreakFusion',ofile)
        else:
            pass
    


def EncodeAnnot(gene):
    try:
        Entrez.email = 'seas363@gmail.com'
        query = '('+gene+'[Gene Name]) AND human[Organism]'
        handle = Entrez.esearch(db = 'gene',term = query)
        record = Entrez.read(handle)
        handle.close()
        ID = record['IdList'][0]
        handle = Entrez.efetch(db = 'gene', id = ID, retmode = 'xml')
        record = Entrez.read(handle)
        handle.close()
        encode = record[0]['Entrezgene_type'].attributes['value'].encode('ascii')
        return encode
    except:
        return ''


    
def AddUniprotAnnotation(ifile,ofile,encodeInfoDict={}):
    f = open(ifile,'r')
    of = open(ofile,'w')
    header = f.readline()
    colName = header.split('\t')
    colName_new = colName[:9] + \
                  ['Encode 1','Protein 1 Uniprot Entry','Protein 1 Subcellular Location','Protein 1 Function'] + \
                  colName[9:15] + \
                  ['Encode 2','Protein 2 Uniprot Entry','Protein 2 Subcellular Location','Protein 2 Function'] + \
                  colName[15:]
    header_new = '\t'.join(colName_new)
    of.write(header_new)

    # Pre-fill for non-gene
    encodeInfoDict['null'] = ['',{}]

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
        b1 = (gene1 in encodeInfoDict)
        b2 = (gene2 in encodeInfoDict)
        if b1:
            encode1 = encodeInfoDict[gene1][0]
            protInfo1 = encodeInfoDict[gene1][1]
        else:
            encode1 = EncodeAnnot(gene1)
            protInfo1 = {}
            if encode1 == 'protein-coding':
                protInfo1 = UniprotInfo(gene1)
            encodeInfoDict[gene1] = [encode1,protInfo1]
        if b2:
            encode2 = encodeInfoDict[gene2][0]
            protInfo2 = encodeInfoDict[gene2][1]
        else:
            encode2 = EncodeAnnot(gene2)
            protInfo2 = {}
            if encode2 == 'protein-coding':
                protInfo2 = UniprotInfo(gene2)
            encodeInfoDict[gene2] = [encode2,protInfo2]
            
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
        
        col_new = col[:9] + [encode1] + p_value1 + col[9:15] + [encode2] + p_value2 + col[15:]
        line_new = '\t'.join(col_new)
        of.write(line_new)


    f.close()
    of.close()
    return encodeInfoDict
