#!/bin/bash

# This script attempts to do all the work to host the site.  It expects to be on Ubuntu 18.04+ but likely also works on 16.04.

set -euo pipefail # exit if an error occurs rather than ignoring it
_readlinkf() { perl -MCwd -le 'print Cwd::abs_path shift' "$1"; } # cross-platform version of `readlink -f`
cd "$(dirname "$(_readlinkf "${BASH_SOURCE[0]}")")" # `cd` to the directory holding this script (which is the root of this git repo)

# Check that needed data is present.  If a missing file can be generated from other files, do that.
# if ! [ -e input_data/gene ]; then
    # echo "please populate input_data/gene/"
    # exit 1
# fi

# if ! [ -e genes-site/assoc.db ]; then
    # if [ -e input_data/variant ]; then
       # python3 genes-site/make_sqlite3_db.py
    # else
        # echo "either populate input_data/variant/ and run ./make_sqlite_db.py or copy assoc.db here"
        # exit 1
    # fi 
# fi

# if ! [ -e genes-site/static/phenotypes.json ] || ! [ -e genes-site/static/genes.json ]; then
    # python3 genes-site/make_tables.py
# fi

# if ! [ -e genes-site/variant.db ]; then
    # if [ -e input_data/variant ]; then
       # python3 genes-site/make_variant_sqlite3_db.py
    # else
        # echo "either populate input_data/variant/ and run ./make_variant_sqlite_db.py or copy variant.db here"
        # exit 1
    # fi
# fi

# Install dependencies
if ! [ -e venv ]; then
    sudo apt update && sudo apt install python3-pip python3-venv nginx
    python3 -m venv venv
    ./venv/bin/pip3 install -r requirements.txt
fi

# Make a Systemd Unit file that runs gunicorn to host the site (available only locally on this machine)
if ! [ -e /etc/systemd/system/gunicorn-genes-site-ukbb200K.service ]; then
    sudo tee /etc/systemd/system/gunicorn-genes-site-ukbb200K.service >/dev/null <<END
[Unit]
Description=Gunicorn instance to serve genes-site-ukbb200K
After=network.target
[Service]
User=nobody
Group=nogroup
WorkingDirectory=$PWD/genes-site/
ExecStart=$PWD/venv/bin/gunicorn -k gevent -w4 --bind localhost:8800 serve:app
[Install]
WantedBy=multi-user.target
END
    sudo systemctl daemon-reload
    sudo systemctl start gunicorn-genes-site-ukbb200K
    sudo systemctl enable gunicorn-genes-site-ukbb200K
fi

# Make nginx reverse-proxy the local-only gunicorn port to an externally-accessible subdomain
if ! [ -e /etc/nginx/sites-enabled/genes-site-ukbb200K ]; then
    sudo tee /etc/nginx/sites-available/genes-site-ukbb200K >/dev/null <<END
server {
    listen 80;
    server_name ukb-200kexome.leelabsg.org;
    location / {
        include proxy_params;
        proxy_pass http://localhost:8800;
    }
}
END
    sudo ln -s /etc/nginx/sites-available/genes-site-ukbb200K /etc/nginx/sites-enabled/
    sudo nginx -t # check that the file is good
    sudo systemctl restart nginx
fi

# Restart gunicorn to apply any changes
sudo systemctl restart gunicorn-genes-site-ukbb200K

# https cert
# sudo certbot --nginx -d ukb-200kexome.leelabsg.org

echo SUCCESS
