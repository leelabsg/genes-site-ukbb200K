#!/usr/bin/env python3

# phenotypes.json - [{phecode:'008', category:'infectious diseases', best_pval:1e-5, best_assoc:'KEGG_FOO', num_sig_assocs:47}, ...]
# genes.json - [{name:'PCSK9', best_pval:1e-5, best_assoc:'008', num_sig_assocs:24}, ...]

# MAYBE: use sqlite instead of json so that the browser can page?

import sqlite3, json, os
conn = sqlite3.connect('assoc.db')
conn.row_factory = sqlite3.Row

pheno_by_id = {}
for row in conn.execute('SELECT * FROM pheno'):
    pheno_by_id[row['id']] = dict(
        phecode=row['phecode'], phenostring=row['phenostring'], category=row['category'],
        best_pval=999, best_assoc=None, num_sig_assocs=0)

gene_by_id = {}
for row in conn.execute('SELECT * FROM gene'):
    gene_by_id[row['id']] = dict(
        name=row['name'],
        best_pval=999, best_assoc=None, num_sig_assocs=0)

for i, row in enumerate(conn.execute('SELECT * FROM assoc')):
    if i % 1_000_000 == 0: print(i)
    pval = row['pval']
    pheno = pheno_by_id[row['pheno_id']]
    gene = gene_by_id[row['gene_id']]
    if pval < pheno['best_pval']: pheno['best_pval'], pheno['best_assoc'] = pval, gene['name']
    if pval < gene['best_pval']: gene['best_pval'], gene['best_assoc'] = pval, pheno['phecode']
    if pval < 1e-4:
        pheno['num_sig_assocs'] += 1
        gene['num_sig_assocs'] += 1

if not os.path.exists('static'): os.mkdir('static')
with open('static/phenotypes.json', 'w') as f:
    json.dump(sorted((p for p in pheno_by_id.values() if p['best_pval'] <= 1), key=lambda x:x['best_pval']), f, separators=(',', ':'))
with open('static/genes.json', 'w') as f:
    json.dump(sorted((p for p in gene_by_id.values() if p['best_pval'] <= 1e-4), key=lambda x:x['best_pval']), f, separators=(',', ':'))
