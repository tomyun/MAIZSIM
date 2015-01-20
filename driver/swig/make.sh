#!/bin/bash
swig -python -c++ -o crop.cc crop.i
python setup.py build_ext --inplace
ln -sf swig/_crop.so ../
ln -sf swig/crop.py ../
