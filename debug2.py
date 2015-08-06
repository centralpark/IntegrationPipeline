from time import time
from Bio import Entrez


Entrez.email = 'seas363@gmail.com'
handle = Entrez.esearch(db = 'pubmed', term = 'biopython')
record = Entrez.read(handle)
record['IdList']