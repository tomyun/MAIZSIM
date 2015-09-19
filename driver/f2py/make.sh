#!/bin/bash
find ../soil/ | grep -e ".*\.\(for\|FOR\)" | sed "/2DMAIZSIM.for/ d" | sed "/solupt.for/ d" | xargs f2py3.5 -c -m soil --f77flags="--std=legacy -fno-align-commons"
ln -sf f2py/soil.so ../
