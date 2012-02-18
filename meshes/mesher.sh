#!/bin/sh

# ./compile.sh
./biventricle
perl -pi -e 's/Dimension 3/Dimension\n3/g' biventricle.mesh
dolfin-convert biventricle.mesh biventricle.xml
python pvdconvert.py