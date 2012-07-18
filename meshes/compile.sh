#!/bin/sh

#g++-mp-4.6 -frounding-math -I/Users/harish/Work/FEniCS/dev/include  -c leftventricle.cpp
#g++-mp-4.6 -L/Users/harish/Work/FEniCS/dev/lib -lcgal -L/opt/local/lib -lmpfr -lgmp -lboost_thread-mt -frounding-math -o leftventricle leftventricle.o

g++-mp-4.6 -frounding-math -I/Users/harish/Work/FEniCS/dev/include  -c biventricle.cpp
g++-mp-4.6 -L/Users/harish/Work/FEniCS/dev/lib -lcgal -L/opt/local/lib -lmpfr -lgmp -lboost_thread-mt -frounding-math -o biventricle biventricle.o
