#!/bin/bash
swig -python -py3 -c++ -o crop.cc crop.i
python3 setup.py build_ext --inplace
ln -sf swig/_crop.so ../
ln -sf swig/crop.py ../
