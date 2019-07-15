#!/usr/bin/env python3
# TODO: use BLOB for variant_df(df)
'''
this script creates `variant.db` which contains all the variants for each (phecode, gene) pair.
the variant data is encoded in dataframe (dictionary of lists) format.
variants are sorted and positions are replaced by pos_delta to save space.
'''

'''
18600 genes
variant/ = 6.7GB in 793 phenos -> 8.5MB each -> result_singlevariant_969.txt.gz is typical
result_singlevariant_969.txt.gz is 8,5M gzipped (6B/line, 600B/gene), 59MB raw (42B/line, 3KB/gene)
result_singlevariant_969.txt.gz has 1.4M lines -> 75 variants/gene
file size comparison:
  969.txt.gz=8.5MB
  969.db=13MB
  969_nochrom.db=13MB
  969_nochrom_posdelta.db=12MB
  969_nochrom_posdelta_zstddict_level4.db=11MB in 14s
  969_nochrom_posdelta_zstddict_level10.db=9.7MB in 18s
  969_nochrom_posdelta_zstdicit_level19.db=8.7MB in 29s
'''

import sqlite3, gzip, csv, os, itertools, json
import zstandard
from boltons.iterutils import pairwise_iter
os.chdir(os.path.dirname(os.path.abspath(__file__)))

db_filepath = 'variant.db'

phecodes = sorted(row[0] for row in sqlite3.connect('assoc.db').execute('SELECT phecode FROM pheno'))
print('found', len(phecodes), 'phecodes')

def get_genes_variantdata():
    for i,phecode in enumerate(phecodes):
        print(' - reading pheno#{}: {}'.format(i, phecode))
        with gzip.open('../input_data/variant/result_singlevariant_{}.txt.gz'.format(phecode), 'rt') as f:
            rows = csv.DictReader(f, delimiter=' ')
            for genename,rowgroup in itertools.groupby(rows, key=lambda r:r['GeneName']):
                rows = sorted(rowgroup, key=lambda r:int(r['SNP'].split(':',2)[1])) # sort by pos
                df = {key:[] for key in 'pos base maf mac_case mac_control pval'.split()}
                chrom = rows[0]['SNP'].split(':')[0]
                for row in rows:
                    chrom_, pos_str, base = row['SNP'].split(':', 2)
                    assert chrom_ == chrom
                    df['pos'].append(int(pos_str))
                    df['base'].append(base)
                    df['maf'].append(float(row['MAF']))
                    df['mac_case'].append(int(row['MAC_Case']))
                    df['mac_control'].append(int(row['MAC_Control']))
                    df['pval'].append(float(row['p.value']))
                # To get better compression, replace df['pos'] with df['pos_delta'], which contains offsets from the previous value.
                # The first position keeps its actual value (ie, its offset from zero).
                df['pos_delta'] = [pos - previous_pos for previous_pos, pos in pairwise_iter([0]+df.pop('pos'))]
                variant_data = json.dumps(df, separators=(',',':')).encode('utf8')
                yield (phecode, genename, chrom, variant_data)

samples = [variant_data for phecode,genename,chrom,variant_data in itertools.islice(get_genes_variantdata(), 0, 18700*4)] # ~4 phenos
# similar to `zstd --train variant-zstd-training/* -o zstd-variant-dictionary`
print('collected samples for zstd dict training')
zstd_dict = zstandard.train_dictionary(131072, samples) # docs use 131072 and `zstd --train` produced 112KB
#with open(zstd_dict_filepath, 'wb') as f: f.write(zstd_dict.as_bytes()) # usable by `zstd -D dict < text > compressed`
print('trained zstd_dict on samples')

zstd_compressor = zstandard.ZstdCompressor(level=8, dict_data=zstd_dict) # level=8 is a decent middleground
def variant_df_gen():
    for phecode, genename, chrom, variant_data in get_genes_variantdata():
        compressed_data = zstd_compressor.compress(variant_data)
        yield (phecode, genename, chrom, compressed_data)
db_tmp_filepath = db_filepath + '.tmp.db'
if os.path.exists(db_tmp_filepath): os.unlink(db_tmp_filepath)
conn = sqlite3.connect(db_tmp_filepath)
with conn:
    conn.execute('CREATE TABLE compression_dict (id INT PRIMARY KEY, data BLOB)')
    conn.execute('INSERT INTO compression_dict (data) VALUES (?)', (zstd_dict.as_bytes(),))
    conn.execute('CREATE TABLE variant_df (id INT PRIMARY KEY, phecode TEXT, genename TEXT, chrom TEXT, df BLOB)')
    conn.executemany('INSERT INTO variant_df (phecode, genename, chrom, df) VALUES (?,?,?,?)', variant_df_gen())
    conn.execute('CREATE INDEX idx_variantdf_phecode_genename ON variant_df (phecode,genename)')
print('finished ', db_tmp_filepath)
if os.path.exists(db_filepath): os.unlink(db_filepath)
os.rename(db_tmp_filepath, db_filepath)
print('moved to ', db_filepath)
