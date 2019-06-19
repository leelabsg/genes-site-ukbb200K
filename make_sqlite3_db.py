#!/usr/bin/env python3

import os, re, gzip, sqlite3, csv, math

def round_sig(x, digits):
    if x == 0:
        return 0
    elif abs(x) == math.inf or math.isnan(x):
        raise ValueError("Cannot round infinity or NaN")
    else:
        log = math.log10(abs(x))
        digits_above_zero = int(math.floor(log))
        return round(x, digits - 1 - digits_above_zero)
assert round_sig(0.00123, 2) == 0.0012
assert round_sig(1.59e-10, 2) == 1.6e-10

# make list of phenotypes
filenames = os.listdir('input_data/gene')
for fname in filenames: assert re.match(r'result_gene_([0-9]{3,4}(?:\.[0-9]{1,2})?).txt.gz$', fname), fname
phecodes = set(re.match(r'result_gene_([0-9]{3,4}(?:\.[0-9]{1,2})?).txt.gz$', fname).group(1) for fname in filenames)
print(len(phecodes), 'phecodes')
# also check variant files
filenames_for_variants = os.listdir('input_data/variant')
for fname in filenames_for_variants: assert re.match(r'result_singlevariant_([0-9]{3,4}(?:\.[0-9]{1,2})?).txt.gz$', fname), fname
assert phecodes == set(re.match(r'result_singlevariant_([0-9]{3,4}(?:\.[0-9]{1,2})?).txt.gz$', fname).group(1) for fname in filenames_for_variants)
phenos = {}
for row in csv.DictReader(open('input_data/phenotype-info.csv')):
    if row['phecode'] not in phecodes: continue
    phenos[row['phecode']] = {
        'phenostring': row['description'],
        'category': row['group'],
        # 'num_cases': row['n_cases'], # wrong.
        # 'num_controls': row['n_controls'], # wrong.
        'phecode_exclude_range': row['phecode_exclude_range'].strip(),
        'phecode_exclude_description': row['phecode_exclude_description'].strip(),
        'sex': row['sex'],
    }
for phecode in phecodes: assert phecode in phenos, phecode
assert 2 <= len(set(p['category'] for p in phenos.values())) < 100
for phecode, pheno in phenos.items():
    with gzip.open('input_data/gene/result_gene_{}.txt.gz'.format(phecode), 'rt') as f:
        for row in csv.DictReader(f, delimiter=' '):
            pheno['num_cases'] = int(row['Case'])
            pheno['num_controls'] = int(row['Control'])
            break

# make list of genes
# file are each missing a few genes.  we get the full 18336 after reading ~5 files.
genes = {}
for phecode in phecodes:
    with gzip.open('input_data/gene/result_gene_{}.txt.gz'.format(phecode), 'rt') as f:
        f.readline() # eat the header
        for line in f:
            genename, start, _ = line.split(' ',2)
            chrom = start.split(':',1)[0]
            assert chrom.isdigit() or chrom
            if genename not in genes: genes[genename] = {'chrom':chrom}
            else: assert genes[genename]['chrom'] == chrom
print(len(genes), 'genes')

# make mapping of (phecode, genename) -> pval
# store in sqlite3 as table (pheno, genename, pval)
phecode_ids = {phecode: id_ for id_, phecode in enumerate(sorted(phecodes))}
gene_ids = {genename: id_ for id_, genename in enumerate(sorted(genes))}
def pheno_row_generator():
    for phecode,id_ in phecode_ids.items():
        p = phenos[phecode]
        yield (id_, phecode, p['phenostring'], p['category'], p['num_cases'], p['num_controls'],
               p['phecode_exclude_range'], p['phecode_exclude_description'], p['sex'])

def gene_row_generator():
    for genename, id_ in gene_ids.items():
        g = genes[genename]
        yield (id_, genename, g['chrom'])

def assoc_row_generator(): # doesn't output primary key, let's sqlite3 auto-increment instead
    for i,(phecode,phecode_id) in enumerate(phecode_ids.items()):
        pheno = phenos[phecode]
        filename = 'result_gene_{}.txt.gz'.format(phecode)
        print(i, filename)
        with gzip.open('input_data/gene/' + filename, 'rt') as f:
            assert f.readline().strip().split() == 'GeneName Start_Pos End_Pos NumofRareVariants MAC_Case MAC_Control Case Control p.value'.split()
            for line in f:
                try:
                    genename, start, end, num_rare_str, mac_case_str, mac_control_str, num_cases_str, num_controls_str, pval_str = line.strip().split(' ')
                    assert genename in genes, line
                    chrom = start.split(':',1)[0]
                    startpos = int(start.split(':',2)[1])
                    assert end.split(':',1)[0] == chrom
                    endpos = int(end.split(':',2)[1]) # sometimes startpos > endpos
                    num_rare = int(num_rare_str); assert 0 < int(num_rare) < 1e6, line
                    mac_case = round(float(mac_case_str), 2); assert 0 <= mac_case < 1e6, line
                    mac_control = round(float(mac_control_str), 2); assert 0 <= mac_control < 1e6, line
                    assert 0 < mac_case + mac_control, line
                    num_cases = int(num_cases_str); assert num_cases == pheno['num_cases'], line
                    num_controls = int(num_controls_str); assert num_controls == pheno['num_controls'], line
                    pval = 1 if pval_str=='NA' else round_sig(float(pval_str), 2); assert 0 < pval <= 1, line
                    yield (phecode_id,gene_ids[genename], pval,num_rare, startpos,endpos, mac_case,mac_control)
                except Exception: raise Exception(line, filename, pheno)

db_fname = 'assoc.db'
db_tmp_fname = db_fname + '.tmp.db'
if os.path.exists(db_tmp_fname): os.unlink(db_tmp_fname)
conn = sqlite3.connect(db_tmp_fname)
with conn: # this commits insertions
    conn.execute('create table pheno (id INTEGER PRIMARY KEY, phecode VARCHAR, phenostring VARCHAR, category VARCHAR, num_cases INTEGER, num_controls INTEGER, phecode_exclude_range VARCHAR, phecode_exclude_description VARCHAR, sex VARCHAR)')
    conn.execute('create table gene (id INTEGER PRIMARY KEY, name VARCHAR, chrom VARCHAR)')
    conn.execute('create table assoc (id INTEGER PRIMARY KEY, pheno_id INTEGER, gene_id INTEGER, pval REAL, startpos INTEGER, endpos INT, num_rare INTEGER, mac_case REAL, mac_control REAL, FOREIGN KEY(pheno_id) REFERENCES pheno(id), FOREIGN KEY(gene_id) REFERENCES gene(id))')

    conn.executemany('INSERT INTO pheno VALUES (?,?,?,?,?,?,?,?,?)', pheno_row_generator())
    conn.executemany('INSERT INTO gene VALUES (?,?,?)', gene_row_generator())
    conn.executemany('INSERT INTO assoc (pheno_id,gene_id, pval,num_rare, startpos,endpos, mac_case,mac_control) VALUES (?,?,?,?,?,?,?,?)', assoc_row_generator())

    conn.execute('CREATE INDEX idx_assoc_pheno_id ON assoc (pheno_id)')
    conn.execute('CREATE INDEX idx_assoc_gene_id ON assoc (gene_id)')

if os.path.exists(db_fname): os.unlink(db_fname)
os.rename(db_tmp_fname, db_fname)
