#!/bin/bash

make dist
make doc
cp vasputil-*.tar.gz ~/public_html/vasputil
cp doc/README.html ~/public_html/vasputil
chmod a+r ~/public_html/vasputil/*
