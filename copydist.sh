#!/bin/bash
# Make the dist directory
rst2html README README.html
git2cl > ChangeLog
./setup.py sdist
cp dist/vasputil-*.tar.gz ~/public_html/vasputil
cp README.html ~/public_html/vasputil
chmod a+r ~/public_html/vasputil/*
