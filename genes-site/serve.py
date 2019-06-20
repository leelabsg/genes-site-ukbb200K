#!/usr/bin/env python3

import sqlite3, re, itertools, json
import zstandard
from flask import g, Flask, jsonify, abort, render_template, request, url_for, redirect
app = Flask(__name__)

app.config['LZJS_VERSION'] = '0.9.0'
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 60*5

DATABASE = 'assoc.db'
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
    return render_template('pheno.html', phecode=phecode, phenostring=phenostring, category=category, num_cases=num_cases, num_controls=num_controls)

@app.route('/assoc/<genename>/<phecode>')
def assoc_page(genename, phecode):
    matches = list(get_db().execute('SELECT id,phenostring,category,num_cases,num_controls FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    pheno_id, phenostring, category, num_cases, num_controls = matches[0]

    matches = list(get_db().execute('SELECT id,chrom FROM gene WHERE name = ?', (genename,)))
    if not matches: return abort(404)
    gene_id,chrom = matches[0]

    matches = list(get_db().execute('SELECT pval,num_rare, startpos,endpos, mac_case,mac_control FROM assoc '
                                    'LEFT JOIN gene ON assoc.gene_id=gene.id '
                                    'WHERE pheno_id=? AND gene_id=?', (pheno_id, gene_id)))
    if not matches: return abort(404)
    m = matches[0]
    return render_template('assoc.html',
                           phecode=phecode, phenostring=phenostring, category=category, num_cases=num_cases, num_controls=num_controls,
                           genename=genename,
                           pval=m[0],num_rare=m[1], chrom=chrom,startpos=m[2],endpos=m[3], mac_case=m[4],mac_control=m[5])


@app.route('/api/gene/<genename>')
def gene_api(genename):
    matches = list(get_db().execute('SELECT id FROM gene WHERE name = ?', (genename,)))
    if not matches: return abort(404)
    gene_id = matches[0][0]
    df = get_df('SELECT phecode,phenostring,category,pval,startpos,endpos,num_rare,mac_case,mac_control,num_cases,num_controls FROM assoc '
                'LEFT JOIN pheno ON assoc.pheno_id=pheno.id '
                'WHERE gene_id=? '
                'ORDER BY phecode', (gene_id,))
    return jsonify(dict(assocs=df))

@app.route('/api/pheno/<phecode>')
def pheno_api(phecode):
    matches = list(get_db().execute('SELECT id FROM pheno WHERE phecode=?', (phecode,)))
    if not matches: return abort(404)
    pheno_id = matches[0][0]
    num_genes = list(get_db().execute('SELECT COUNT(*) FROM assoc WHERE pheno_id=? ', (pheno_id,)))[0][0]
    df = get_df('SELECT name,pval,startpos,endpos,num_rare,mac_case,mac_control,chrom FROM assoc '
                'LEFT JOIN gene ON assoc.gene_id=gene.id '
                'WHERE pheno_id=? '
                'ORDER BY pval LIMIT 2000', (pheno_id,))
    return jsonify(dict(assocs=df, num_genes=num_genes))


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
        df = json.loads(decompressed)
        return dict(phecode=m['phecode'], genename=m['genename'], chrom=m['chrom'], df=df)


class Autocompleter:
    def __init__(self):
        phenos_df = get_df('SELECT phecode,phenostring FROM pheno')
        self.phenos = [{key: phenos_df[key][i] for key in phenos_df} for i in range(len(next(iter(phenos_df.values()))))]
        self.phenos.sort(key=lambda p:float(p['phecode']))
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
