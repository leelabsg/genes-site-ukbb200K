#!/bin/bash
set -euo pipefail
_readlinkf() { perl -MCwd -le 'print Cwd::abs_path shift' "$1"; }
cd "$(dirname "$(_readlinkf "${BASH_SOURCE[0]}")")"

git pull --ff-only
echo

sudo systemctl restart genes-site-ukbb200K
echo

echo SUCCESSFUL UPDATE
