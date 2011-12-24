#!/bin/sh

g++-mp-4.4 -frounding-math -I/Users/harish/Work/FEniCS/dev/include  -c leftventricle.cpp
g++-mp-4.4 -L/Users/harish/Work/FEniCS/dev/lib -lcgal -L/opt/local/lib -lmpfr -lgmp -lboost_thread-mt -frounding-math -o leftventricle leftventricle.o