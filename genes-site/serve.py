#!/usr/bin/env python3

from io import BytesIO
from gzip import GzipFile
import sqlite3, re, itertools, json, math
import zstandard
from flask import g, Flask, jsonify, abort, render_template, request, url_for, redirect, send_from_directory
app = Flask(__name__)

app.config['LZJS_VERSION'] = '0.9.1'
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 60*5

# change by SLEE
DATABASE = 'assoc_new.db'
#
def get_db():
    db = getattr(g, '_database', None)
    if db is None: db = g._database = sqlite3.connect(DATABASE)
    return db
@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None: db.close()
def get_df(query, args=()):
    '''get a dataframe (eg, `{phecode:['008','008.5'],pval:[0.3,0.01]}`) from the database'''
    cur = get_db().execute(query, args)
    colnames = [x[0] for x in cur.description]
    rows = cur.fetchall()
    cur.close()
    return {colname: [row[i] for row in rows] for i, colname in enumerate(colnames)}
def get_samplesize(num_case, num_control):
    samplesize = str(num_case)
    if( num_control != None ):
        samplesize += ':' + str(num_control)
    return samplesize
    
def get_samplesizes(num_cases, num_controls):
    samplesizes = []
    for n0,n1 in zip(num_cases,num_controls):
        sample_item = get_samplesize(n0, n1)
        samplesizes.append(sample_item)
    return samplesizes
        
@app.route('/')
def index_page():
    return render_template('index.html')
@app.route('/about')
def about_page():
    return render_template('about.html')
@app.route('/phenotypes')
def phenotypes_page():
    return render_template('phenotypes.html')
@app.route('/genes')
def genes_page():
    return render_template('genes.html')

@app.route('/gene/<genename>')
def gene_page(genename):
    matches = list(get_db().execute('SELECT chrom FROM gene WHERE name = ?', (genename,)))
    if not matches: return abort(404)
    chrom = matches[0][0]
    return render_template('gene.html', genename=genename, chrom=chrom)

@app.route('/pheno/<phecode>')
def pheno_page(phecode):
    matches = list(get_db().execute('SELECT phenostring,category,num_cases,num_controls FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    phenostring, category, num_cases, num_controls = matches[0]
    sample_size = get_samplesize(num_cases,num_controls)
    return render_template('pheno.html', phecode=phecode, phenostring=phenostring, category=category, sample_size=sample_size)

@app.route('/assoc/<genename>/<phecode>')
def assoc_page(genename, phecode):
    matches = list(get_db().execute('SELECT id,phenostring,category,num_cases,num_controls FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    pheno_id, phenostring, category, num_cases, num_controls = matches[0]
    sample_size = get_samplesize(num_cases,num_controls)
    
    matches = list(get_db().execute('SELECT id,chrom FROM gene WHERE name = ?', (genename,)))
    if not matches: return abort(404)
    gene_id,chrom = matches[0]

    matches = list(get_db().execute('SELECT pval,num_rare, startpos,endpos, assoc.id FROM assoc '
                                    'LEFT JOIN gene ON assoc.gene_id=gene.id '
                                    'WHERE pheno_id=? AND gene_id=?', (pheno_id, gene_id)))
    if not matches: return abort(404)
    m = matches[0]
    return render_template('assoc.html',
                           phecode=phecode, phenostring=phenostring, category=category, num_cases=num_cases, num_controls=num_controls,
                           sample_size=sample_size,
                           genename=genename,
                           pval=m[0],num_rare=m[1], chrom=chrom,startpos=m[2],endpos=m[3], associd =m[4])

@app.route('/api/assocgroup/<associd>')
def assocgroup_api(associd):
    df = get_df('SELECT description, pval, mac, rarevariants, ultra_rarevariants, pval_collapsed_ultrarare  FROM assoc_group WHERE assoc_id=?', (associd,))
    return jsonify(dict(assocgroup=df))

#added by WZ to have a seperate page for association rsuults by annotations
@app.route('/assocgroup/<associd>')
def assocgroup_page(associd):
    matches = list(get_db().execute('SELECT description, pval, mac, rarevariants, ultra_rarevariants, pval_collapsed_ultrarare  FROM assoc_group WHERE assoc_id=?', (associd,)))
    if not matches: return abort(404)
    return render_template('assocgroup.html', associd = associd)

@app.route('/download/pheno/<phecode>')
def download_pheno(phecode):
    matches = list(get_db().execute('SELECT id FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    return send_from_directory('../input_data/gene', 'result_gene_{}.txt.gz'.format(phecode), as_attachment=True)

@app.route('/api/gene/<genename>')
def gene_api(genename):
    matches = list(get_db().execute('SELECT id FROM gene WHERE name = ?', (genename,)))
    if not matches: return abort(404)
    gene_id = matches[0][0]
    df = get_df('SELECT phecode,phenostring,category,pval,startpos,endpos,num_cases,num_controls FROM assoc '
                'LEFT JOIN pheno ON assoc.pheno_id=pheno.id '
                'WHERE gene_id=? '
                'ORDER BY phecode', (gene_id,))
    sample_size = get_samplesizes(df['num_cases'], df['num_controls'])
    df['sample_size'] = sample_size
    return jsonify(dict(assocs=df))

@app.route('/api/pheno/<phecode>')
def pheno_api(phecode):
    matches = list(get_db().execute('SELECT id FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    pheno_id = matches[0][0]
    num_genes = list(get_db().execute('SELECT COUNT(*) FROM assoc WHERE pheno_id=? ', (pheno_id,)))[0][0]
    # change by SLEE, remove num_rare, mac_case, mac_control..
    df = get_df('SELECT name,pval,startpos,endpos,chrom FROM assoc '
                'LEFT JOIN gene ON assoc.gene_id=gene.id '
                'WHERE pheno_id=? '
                'ORDER BY pval', (pheno_id,))
                
    ##SLEE remove character after _ (in refseq)
    chrom_new =[]
    for chrom_item in df['chrom']:
        chrom_new.append(chrom_item.split("_",1)[0])
    df['chrom'] = chrom_new
    ####
    #print(df['pval'][999:1002])
    unbinned, manhattan_bins = make_manhattan_bins(df)
    return jsonify(dict(assocs=unbinned, manhattan_bins=manhattan_bins, num_genes=num_genes))
def make_manhattan_bins(df):
    #SLEE, remove num_rare, mac_case, mac_control 
    #assert set(df.keys()) == set('name pval startpos endpos num_rare mac_case mac_control chrom'.split())
    assert set(df.keys()) == set('name pval startpos endpos chrom'.split())
    unbinned = {key:values[:1000] for key,values in df.items()}
    # objs_to_bin = [{key:df[key][i] for key in df.keys()} for i in range(1000, len(df['pval']))]
    # objs_to_bin.sort(key=lambda x:x['
    pos_bin_size = int(4e6)
    nlpval_bin_size = 0.05
    nlpvals_by_chrompos = {} # {('1',0): [0, 0.05, 0.1, 0.45], ('1',3e6): [2.4]}
    for i in range(1000, len(df['pval'])):
        chrom = df['chrom'][i]
        binned_pos = int(df['startpos'][i] // pos_bin_size * pos_bin_size + pos_bin_size/2)
        nlpval = 999 if df['pval'][i]==0 else -math.log10(df['pval'][i])
        binned_nlpval = nlpval // nlpval_bin_size * nlpval_bin_size + nlpval_bin_size/2
        binned_nlpval = round(binned_nlpval, 3) # trim `0.35000000000000003` to `0.35` for convenience and network request size
        nlpvals_by_chrompos.setdefault((chrom, binned_pos), set()).add(binned_nlpval)
    bins = []
    for (chrom,pos), nlpvals in nlpvals_by_chrompos.items():
        nlpvals = sorted(nlpvals)
        extents = [[nlpvals[0], nlpvals[0]]]
        for nlpval in nlpvals:
            if nlpval < extents[-1][1] + nlpval_bin_size * 1.1: extents[-1][1] = nlpval
            else: extents.append([nlpval, nlpval])
        result_nlpvals, result_nlpval_extents = [], []
        for (start, end) in extents:
            if start == end: result_nlpvals.append(start)
            else: result_nlpval_extents.append([start,end])
        bins.append(dict(chrom=chrom, pos=pos, nlpvals=result_nlpvals, nlpval_extents=result_nlpval_extents))
    bins.sort(key=lambda x:(int(x['chrom']) if x['chrom'].isdigit() else 99+hash(x['chrom']), x['pos']))
    return (unbinned, bins)


@app.route('/api/variants/<genename>/<phecode>')
def variants_api(genename, phecode):
    matches = list(get_db().execute('SELECT id FROM gene WHERE name = ?', (genename,)))
    if not matches: return abort(404)
    matches = list(get_db().execute('SELECT id FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    variant_fetcher = getattr(g, '_variant_fetcher', None)
    if variant_fetcher is None: variant_fetcher = g._variant_fetcher = VariantFetcher()
    x = variant_fetcher.get(phecode, genename)
    return jsonify(x)
class VariantFetcher:
    def __init__(self):
        self.db = sqlite3.connect('variant.db')
        self.db.row_factory = sqlite3.Row
        zstd_dict = list(self.db.execute('SELECT data FROM compression_dict'))[0]['data']
        self.zstd_decompressor = zstandard.ZstdDecompressor(dict_data=zstandard.ZstdCompressionDict(zstd_dict))
    def get(self, phecode, genename):
        matches = list(self.db.execute('SELECT phecode,genename,chrom,df FROM variant_df WHERE phecode=? AND genename=?', (phecode,genename)))
        if len(matches) != 1: raise Exception('VariantFetcher got {} matches: {}'.format(len(matches), repr(matches)))
        m = matches[0]
        decompressed = self.zstd_decompressor.decompress(m['df'])
        #df = json.loads(decompressed)
        df = json.loads(decompressed.decode('ascii'))
        return dict(phecode=m['phecode'], genename=m['genename'], chrom=m['chrom'], df=df)


class Autocompleter:
    def __init__(self):
        phenos_df = get_df('SELECT phecode,phenostring FROM pheno')
        self.phenos = [{key: phenos_df[key][i] for key in phenos_df} for i in range(len(next(iter(phenos_df.values()))))]
        #self.phenos.sort(key=lambda p:float(p['phecode']))
        for p in self.phenos: p['phenostring--processed'] = self.process_string(p['phenostring'])
        genenames = sorted(get_df('SELECT name FROM gene')['name'])
        self.genes = [{'genename': name} for name in genenames]
        for p in self.genes: p['genename--processed'] = self.process_string(p['genename'])
    non_word_regex = re.compile(r"(?:_|[^\w\.])") # Most of the time we want to include periods in words but not underscores
    def process_string(self, string):
        # Cleaning inspired by <https://github.com/seatgeek/fuzzywuzzy/blob/6353e2/fuzzywuzzy/utils.py#L69>
        return ' ' + self.non_word_regex.sub(' ', string).lower().strip()
    def get_completions(self, query):
        processed_query = self.process_string(query) # replace junk with spaces and use lower-case
        for f in [self.get_completions_on_phecode, self.get_completions_on_phenostring, self.get_completions_on_genename]:
            results = list(itertools.islice(f(processed_query), 0, 10))
            if results: return results
        return []
    def get_best_completion(self, query):
        completions = self.get_completions(query)
        if not completions: return None
        return completions[0] # TODO: make this better, perhaps by checking for exact matches.
    def get_completions_on_phecode(self, processed_query):
        processed_query = processed_query.strip()
        if not re.match(r'^[0-9]+(?:\.[0-9]*)?$', processed_query): return
        for p in self.phenos:
            if p['phecode'].startswith(processed_query):
                yield {
                    'value': p['phecode'],
                    'display': '{} ({})'.format(p['phecode'], p['phenostring']),
                    'url': url_for('pheno_page', phecode=p['phecode'])
                }
    def get_completions_on_phenostring(self, processed_query):
        if len(processed_query) == 0: return
        for p in self.phenos:
            if processed_query in p['phenostring--processed']:
                yield {
                    'value': p['phenostring'],
                    'display': '{} ({})'.format(p['phenostring'], p['phecode']),
                    'url': url_for('pheno_page', phecode=p['phecode'])
                }
    def get_completions_on_genename(self, processed_query):
        if len(processed_query) == 0: return
        for p in self.genes:
            if processed_query in p['genename--processed']:
                yield {
                    'value': p['genename'],
                    'display': p['genename'],
                    'url': url_for('gene_page', genename=p['genename']),
                }

def get_autocompleter():
    a = getattr(g, '_autocompleter', None)
    if a is None: a = g._autocompleter = Autocompleter()
    return a
@app.route('/api/autocomplete')
def autocomplete_api():
    '''generate suggestions for the searchbox'''
    query = request.args.get('query', '')
    suggestions = get_autocompleter().get_completions(query)
    if suggestions:
        return jsonify(sorted(suggestions, key=lambda sugg: sugg['display']))
    return jsonify([])
@app.route('/go')
def go():
    '''attempt to send the user to a relevant page after they hit enter on the searchbox'''
    query = request.args.get('query', None)
    if query is None: return redirect(url_for('index_page'))
    best_suggestion = get_autocompleter().get_best_completion(query)
    if not best_suggestion: return redirect(url_for('index_page'))
    return redirect(best_suggestion['url'])


class Compress(object):
    '''A copy of Flask-Compress to avoid its packaging issues'''
    def __init__(self, app):
        self.app = app
        app.config.setdefault('COMPRESS_MIMETYPES', ['text/html', 'text/css', 'text/xml', 'application/json', 'application/javascript'])
        app.config.setdefault('COMPRESS_LEVEL', 6)
        app.config.setdefault('COMPRESS_MIN_SIZE', 500)
        app.after_request(self.after_request)
    def after_request(self, response):
        accept_encoding = request.headers.get('Accept-Encoding', '')
        if (response.mimetype not in self.app.config['COMPRESS_MIMETYPES'] or
            'gzip' not in accept_encoding.lower() or
            not 200 <= response.status_code < 300 or
            (response.content_length is not None and response.content_length < self.app.config['COMPRESS_MIN_SIZE']) or
            'Content-Encoding' in response.headers):
            return response
        response.direct_passthrough = False
        gzip_content = self.compress(self.app.config['COMPRESS_LEVEL'], response)
        response.set_data(gzip_content)
        response.headers['Content-Encoding'] = 'gzip'
        response.headers['Content-Length'] = response.content_length
        vary = response.headers.get('Vary')
        if not vary or 'accept-encoding' not in vary.lower():
            response.headers['Vary'] = (vary+', ' if vary else '') + 'Accept-Encoding'
        return response
    def compress(self, level, response):
        gzip_buffer = BytesIO()
        with GzipFile(mode='wb', compresslevel=level, fileobj=gzip_buffer) as gzip_file:
            gzip_file.write(response.get_data())
        return gzip_buffer.getvalue()
Compress(app)


if __name__ == '__main__':
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 1
    app.run(port=5000, debug=True)
