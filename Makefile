# Voro++ makefile
#
# Original Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011
#
# Modified by Paula C. Sanematsu (Syracuse University, Manning Lab)
# Email  : pcsanema@syr.edu
# Date   : May 2021

# Load the common configuration file
include voro++-0.4.6/config.mk

# List of executables
EXECUTABLES=random_points_vtk 

# Makefile rules
#all: $(EXECUTABLES)
all:
	$(MAKE) -C voro++-0.4.6/src
	$(MAKE) -C voro++-0.4.6/examples
	$(MAKE) random_points_vtk

random_points_vtk: random_points_vtk.cc
	$(CXX) $(CFLAGS) $(E_INC_VTK) $(E_LIB_VTK) -o random_points_vtk random_points_vtk.cc -lvoro++ -lvtkCommonCore-7.1

clean:
	$(MAKE) -C voro++-0.4.6/src clean
	$(MAKE) -C voro++-0.4.6/examples clean
	rm -f $(EXECUTABLES)

.PHONY: all clean
