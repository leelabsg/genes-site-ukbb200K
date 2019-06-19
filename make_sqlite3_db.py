#!/usr/bin/env python3

import os, re, gzip, sqlite3, csv, math, random

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
for row in csv.DictReader(open('/Users/peter/PROJECTS/gauss-site/pathways-diptavo/input_data/phenotype-colors.csv')):
    def _int(s): num = int(s.replace(',','')); assert num >= 0, s; return num
    if row['PheCode'] not in phecodes: continue
    phenos[row['PheCode']] = {
        'phenostring': row['Phenotype.Description'],
        'category': row['Phenotype.Category'],
    }
for phecode in phecodes:
    if phecode not in phenos:
        phenos[phecode] = {
            'phenostring': 'no-name-{}'.format(random.random()),
            'category': 'missing-category',
        }
for phecode in phecodes: assert phecode in phenos, phecode
assert 2 <= len(set(p['category'] for p in phenos.values())) < 100


# make list of genes
# file are each missing a few genes.  we get the full 18336 after reading ~5 files.
genenames = set()
for fname in filenames[:50]:
    with gzip.open('input_data/gene/'+fname, 'rt') as f:
        f.readline() # eat the header
        for i,line in enumerate(f):
            genename = line[:line.index(' ')]
            genenames.add(genename)
    #print(len(genenames), i, fname)
print(len(genenames), 'genes')

# make mapping of (phecode, genename) -> pval
# store in sqlite3 as table (pheno, genename, pval)
phecode_ids = {phecode: id_ for id_, phecode in enumerate(sorted(phecodes))}
gene_ids = {genename: id_ for id_, genename in enumerate(sorted(genenames))}
def pheno_row_generator():
    for phecode,id_ in phecode_ids.items():
        p = phenos[phecode]
        yield (id_, phecode, p['phenostring'], p['category'])

def gene_row_generator():
    for genename, id_ in gene_ids.items():
        yield (id_, genename)

def assoc_row_generator(): # doesn't output primary key, let's sqlite3 auto-increment instead
    for i,(phecode,phecode_id) in enumerate(phecode_ids.items()):
        filename = 'result_gene_{}.txt.gz'.format(phecode)
        print(i, filename)
        with gzip.open('input_data/gene/' + filename, 'rt') as f:
            assert f.readline().strip().split() == 'GeneName Start_Pos End_Pos NumofRareVariants MAC_Case MAC_Control Case Control p.value'.split()
            for line in f:
                try:
                    genename, start, end, num_rare_str, mac_case_str, mac_control_str, num_cases_str, num_controls_str, pval_str = line.strip().split(' ')
                    assert genename in genenames, line
                    num_rare = int(num_rare_str); assert 0 < int(num_rare) < 1e6, line
                    mac_case = round(float(mac_case_str), 2); assert 0 <= mac_case < 1e6, line
                    mac_control = round(float(mac_control_str), 2); assert 0 <= mac_control < 1e6, line
                    assert 0 < mac_case + mac_control, line
                    num_cases = int(num_cases_str); assert 0 < num_cases < 1e6, line
                    num_controls = int(num_controls_str); assert 0 < num_controls < 1e6, line
                    if pval_str == 'NA': print('    - warning: pval=NA on line [[{}]]'.format(line)); continue
                    pval = round_sig(float(pval_str), 2); assert 0 < pval <= 1, line
                    yield (phecode_id,gene_ids[genename], pval,num_rare, start,end, mac_case,mac_control, num_cases,num_controls)
                except Exception: raise Exception(line, filename)

db_fname = 'assoc.db'
#if os.path.exists(db_fname): raise Exception(db_fname + ' already exists, please delete')
if os.path.exists(db_fname): os.unlink(db_fname)
conn = sqlite3.connect(db_fname)
with conn: # this commits insertions
    conn.execute('create table pheno (id INTEGER PRIMARY KEY, phecode VARCHAR, phenostring VARCHAR, category VARCHAR)')
    conn.execute('create table gene (id INTEGER PRIMARY KEY, name VARCHAR)')
    conn.execute('create table assoc (id INTEGER PRIMARY KEY, pheno_id INTEGER, gene_id INTEGER, pval REAL, start VARCHAR, end VARCHAR, num_rare INTEGER, mac_case REAL, mac_control REAL, num_cases INTEGER, num_controls INTEGER, FOREIGN KEY(pheno_id) REFERENCES pheno(id), FOREIGN KEY(gene_id) REFERENCES gene(id))')

    conn.executemany('INSERT INTO pheno VALUES (?,?,?,?)', pheno_row_generator())
    conn.executemany('INSERT INTO gene VALUES (?,?)', gene_row_generator())
    conn.executemany('INSERT INTO assoc (pheno_id,gene_id, pval,num_rare, start,end, mac_case,mac_control, num_cases,num_controls) VALUES (?,?,?,?,?,?,?,?,?,?)', assoc_row_generator())

    conn.execute('CREATE INDEX idx_assoc_pheno_id ON assoc (pheno_id)')
    conn.execute('CREATE INDEX idx_assoc_gene_id ON assoc (gene_id)')
