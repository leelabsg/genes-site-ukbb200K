### Usage
Modified for UKBB200K SAIGE-GENE+ output. Database files are directly generated using data_process.R file. 

1. Run data_process.R to generate assoc_new.db file. For this, you need SAIGE-GENE+ output file, previous assoc.db (UKBB 50K analysis results) and refseq annotation file. 

2. Run `python3 make_tables.py` to produce `static/phenotypes.json` and `static/genes.json`

3. Run the server (see below)

### OLD Usage

1. Clone this repository (with `git clone`)

2. Put data in `input_data/`. There are 3 type of files:
   - `input_data/gene/result_gene_<phecode>.txt.gz`
      - format: tab-delimited, including a header row
         - columns: `GeneName`, `Start_Pos`, `End_Pos`, `NumofRareVariants`, `MAC_Case`, `MAC_Control`, `Case`, `Control`, `p.value`
         - eg, `A3GALT2 1:33306784:I:6 1:33312867:C:T 44 11 732.000087726994 114 11286 0.365819488062416`
   - `input_data/variant/result_singlevariant_<phecode>.txt.gz`
      - format: tab-delimited, including a header row
         - columns: `GeneName`, `SNP`, `MAF`, `MAC_Case`, `MAC_Control`, `p.value`
         - eg, `A3GALT2 1:33306784:I:6 1.1e-04 0 5 8.3e-01`
   - `input_data/phenotype-info.csv` (already included in repo)

3. Run `python3 -m pip install -r requirements.txt` (or, if running for only your own user, `python3 -m pip --user install -r requirements.txt`)

4. Run `python3 make_sqlite3_db.py` (inside the `genes-site` subdirectory) to produce `assoc.db`

5. Run `python3 make_tables.py` to produce `static/phenotypes.json` and `static/genes.json`

6. Run `python3 make_variant_sqlite3_db.py` to produce `variant.db`

7. Run the server using one of these:
   - `python3 serve.py` (insecure and slow, for development/debugging)
   - `gunicorn serve:app -k gevent -w4 --bind 0.0.0.0:5000` (fast, for production)

   Then open [http://localhost:5000](http://localhost:5000) in your web browser to access the site. If you are running the site on a remote computer, first start an ssh tunnel (like `ssh -L 5000:localhost:5000 <user>@<server>`) if on Mac/Linux or follow directions if on Windows.

8. To deploy your changes to the server:
   - If you changed the data in `input_data`, copy them to the server's directory `/home/pjvandehaar/genes-site/input_data`.
   - If you changed the files `assoc.db`, `static/phenotypes.json`, `static/genes.json`, or `variant.db`, copy them to `/home/pjvandehaar/genes-site/genes-site` on the server (and replace the files currently there), or just run the commands to re-generate them on the server.
   - Run `git push` to push your code to GitHub, and then on the server run `cd /home/pjvandehaar/genes-site && sudo git pull && sudo systemctl restart gunicorn-genes-site`.
