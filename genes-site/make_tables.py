#!/usr/bin/env python3

# phenotypes.json - [{phecode:'008', category:'infectious diseases', best_pval:1e-5, best_assoc:'KEGG_FOO', num_sig_assocs:47}, ...]
# genes.json - [{name:'PCSK9', best_pval:1e-5, best_assoc:'008', num_sig_assocs:24}, ...]

# MAYBE: use sqlite instead of json so that the browser can page?

import sqlite3, json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

conn = sqlite3.connect('assoc_new.db')
#conn = sqlite3.connect('assoc.db')
conn.row_factory = sqlite3.Row

def get_samplesize(num_case, num_control):
    samplesize = str(num_case)
    if( num_control != None ):
        samplesize += ':' + str(num_control)
    return samplesize
def get_samplesizes(num_cases, num_controls):
    samplesizes = []
    for n0,n1 in zip(num_cases,num_controls):
        sample_item = get_samplesize(n0, n1)
        samplesize.append(sample_item)
    return samplesize
    
pheno_by_id = {}
for row in conn.execute('SELECT * FROM pheno'):
    num_cases=row['num_cases']
    num_controls=row['num_controls']
    sample_size = get_samplesize(num_cases, num_controls)
    pheno_by_id[row['id']] = dict(
        phecode=row['phecode'], phenostring=row['phenostring'].strip(), category=row['category'], num_cases=row['num_cases'], num_controls=row['num_controls'],
        #Change by SLEE
        #phecode_exclude_range=row['phecode_exclude_range'], phecode_exclude_description=row['phecode_exclude_description'], 
        #Add sample size, for quantitative traits. For binary traits, samplesize=case:control
        sample_size=sample_size,
        sex=row['sex'], best_pval=999, best_assoc=None, num_sig_assocs=0)

gene_by_id = {}
for row in conn.execute('SELECT * FROM gene'):
    gene_by_id[row['id']] = dict(
        name=row['name'], chrom=row['chrom'],
        best_pval=999, best_assoc=None, num_sig_assocs=0)

for i, row in enumerate(conn.execute('SELECT * FROM assoc')):
    if i % 1_000_000 == 0: print(i)
    pval = row['pval']
    pheno = pheno_by_id[row['pheno_id']]
    gene = gene_by_id[row['gene_id']]
    if pval < pheno['best_pval']: pheno['best_pval'], pheno['best_assoc'] = pval, gene['name']
    if pval < gene['best_pval']: gene['best_pval'], gene['best_assoc'] = pval, pheno['phecode']
    if pval <= 2.5e-6: #SLEE, change to 2.5e-6
        pheno['num_sig_assocs'] += 1
        gene['num_sig_assocs'] += 1

if not os.path.exists('static'): os.mkdir('static')
with open('static/phenotypes.json', 'w') as f:
    json.dump(sorted((p for p in pheno_by_id.values() if p['best_pval'] <= 1), key=lambda x:x['best_pval']), f, separators=(',', ':'))
with open('static/genes.json', 'w') as f:
    json.dump(sorted((p for p in gene_by_id.values() if p['best_pval'] <= 2.5e-6), key=lambda x:x['best_pval']), f, separators=(',', ':')) #SLEE, change to 2.5e-6
