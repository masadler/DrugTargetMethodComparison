import urllib
from urllib.request import urlopen

import io 

PATH = ""

filename = PATH + "genome_annotation.tsv"
genetype = "protein_coding"
genetype = genetype.replace(" ","")
        
cmd = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"/><Filter name = "biotype" value = "'+genetype+'"/><Attribute name = "ensembl_gene_id" /><Attribute name = "chromosome_name" /><Attribute name = "transcript_start" /><Attribute name = "transcript_end" /><Attribute name = "strand" /><Attribute name = "external_gene_name" /><Attribute name = "band"/></Dataset></Query>'
url = 'http://grch37.ensembl.org/biomart/martservice?query='

with urlopen(url+urllib.parse.quote(cmd)) as f:
          
    charset = f.info().get_content_charset()
    
    html = f.read().decode(charset)

    with open(filename, 'w') as f:
        f.write(html)