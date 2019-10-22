### Usage

1. Clone this repository (with `git clone`)

2. Put data in `input_data/`. There are 3 type of files:
   - `input_data/gene/result_gene_<phecode>.txt.gz`
      - format: space-delimited, including a header row
         - columns: `GeneName`, `Start_Pos`, `End_Pos`, `NumofRareVariants`, `MAC_Case`, `MAC_Control`, `Case`, `Control`, `p.value`
         - eg, `A3GALT2 1:33306784:I:6 1:33312867:C:T 44 11 732.000087726994 114 11286 0.365819488062416`
   - `input_data/variant/result_singlevariant_<phecode>.txt.gz`
      - format: space-delimited, including a header row
         - columns: `GeneName`, `SNP`, `MAF`, `MAC_Case`, `MAC_Control`, `p.value`
         - eg, `A3GALT2 1:33306784:I:6 1.1e-04 0 5 8.3e-01`
   - `input_data/phenotype-info.csv` (already included in repo)

3. run `python3 -m pip install -r requirements.txt` (or, if running for only your own user, `python3 -m pip --user install -r requirements.txt`)

4. run `python3 make_sqlite3_db.py` to produce `assoc.db`

5. run `python3 make_tables.py` to produce `static/phenotypes.json` and `static/genes.json`

6. run `python3 make_variant_sqlite3_db.py` to produce `variant.db`

7. run the server using one of these:
   - `./serve.py` (insecure and slow, for development/debugging)
   - `gunicorn serve:app -k gevent -w4 --bind 0.0.0.0:8000` (fast, for production)
