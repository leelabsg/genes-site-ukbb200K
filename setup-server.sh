#!/bin/bash
# This script attempts to do all the work to host the site.  It expects to be on Ubuntu 18.04+ but likely also works on 16.04.
set -euo pipefail # notify of errors rather than ignoring them
_readlinkf() { perl -MCwd -le 'print Cwd::abs_path shift' "$1"; }
cd "$(dirname "$(_readlinkf "${BASH_SOURCE[0]}")")"

if ! [ -e genes-site/assoc.db ]; then
    if [ -e input_data/gene ] && [ -e input_data/variant ]; then
       python3 genes-site/make_sqlite3_db.py
    else
        echo "either populate input_data and run ./make_sqlite_db.py or copy assoc.db here"
        exit 1
    fi
fi

if ! [ -e genes-site/static/phenotypes.json ] || ! [ -e genes-site/static/genes.json ]; then
    python3 genes-site/make_tables.py
fi

if ! [ -e genes-site/variant.db ]; then
    if [ -e input_data/variant ]; then
       python3 genes-site/make_variant_sqlite3_db.py
    else
        echo "either populate input_data and run ./make_variant_sqlite_db.py or copy variant.db here"
        exit 1
    fi
fi

if ! [ -e venv ]; then
    sudo apt update && sudo apt install python3-pip python3-venv nginx
    python3 -m venv venv
    ./venv/bin/pip3 install -r requirements.txt
fi

if ! [ -e /etc/systemd/system/gunicorn-genes-site.service ]; then
    sudo tee /etc/systemd/system/gunicorn-genes-site.service >/dev/null <<END
[Unit]
Description=Gunicorn instance to serve genes-site
After=network.target
[Service]
User=nobody
Group=nogroup
WorkingDirectory=$PWD/genes-site/
ExecStart=$PWD/venv/bin/gunicorn -k gevent -w4 --bind localhost:8897 serve:app
[Install]
WantedBy=multi-user.target
END
    sudo systemctl daemon-reload
    sudo systemctl start gunicorn-genes-site
    sudo systemctl enable gunicorn-genes-site
fi

if ! [ -e /etc/nginx/sites-enabled/genes-site ]; then
    sudo tee /etc/nginx/sites-available/genes-site >/dev/null <<END
server {
    listen 443;
    location / {
        include proxy_params;
        proxy_pass http://localhost:8897;
    }
}
END
    sudo ln -s /etc/nginx/sites-available/genes-site /etc/nginx/sites-enabled/
    sudo nginx -t # check that the file is good
    sudo systemctl restart nginx
fi

sudo systemctl restart gunicorn-genes-site

echo SUCCESS
