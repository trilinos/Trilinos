# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_ARCH)

MY_COMPILER_FLAGS = -DHAVE_CONFIG_H $(CXXFLAGS) -c -g\
	            -I$(TRILINOS_LIB)/include/
		    
# Adjust to your installation of Zoltan!
ZOLTAN_LINKER_FLAGS = \
	-L/home/u/eboman/zoltan/Zoltan.bugfix/Obj_linux -lzoltan \
	-L/Net/local/proj/zoltan/arch/linux/lib/ParMETIS3 -lparmetis -lmetis

MY_LINKER_FLAGS = $(LDFLAGS) $(TEST_C_OBJ) \
	-L$(TRILINOS_LIB)/lib/ \
	-ltriutils -lepetra -lepetraext $(LIB_PATHS) $(ARCH_LIBS) \
	-llapack -lf77blas -latlas \
	$(ZOLTAN_LINKER_FLAGS)

all: ex1 

ex1: ex1.cpp
	$(CXX)     ex1.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex1.o   $(MY_LINKER_FLAGS)    -o ex1.exe
	/bin/rm -r ex1.o

ex2: ex2.cpp
	$(CXX)     ex2.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex2.o   $(MY_LINKER_FLAGS)    -o ex2.exe
	/bin/rm -r ex2.o

clean:
	/bin/rm -f *.o ex?.exe
