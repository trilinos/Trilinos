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

MY_LINKER_FLAGS = $(LDFLAGS) $(TEST_C_OBJ) \
	-L$(TRILINOS_LIB)/lib/ \
	-ltriutils -lepetra $(LIB_PATHS) $(ARCH_LIBS) -llapack -lblas

all: ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 ex10 \
     ex11 ex12 ex13 ex14 ex15 ex16 ex17 ex18 ex19 ex20\
     ex21 ex22 ex23 ex24 ex25
     
ex1: ex1.cpp
	$(CXX)     ex1.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex1.o   $(MY_LINKER_FLAGS)    -o ex1.exe
	/bin/rm -r ex1.o

ex2: ex2.cpp
	$(CXX)     ex2.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex2.o   $(MY_LINKER_FLAGS)    -o ex2.exe
	/bin/rm -r ex2.o

ex3: ex3.cpp
	$(CXX)     ex3.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex3.o   $(MY_LINKER_FLAGS)    -o ex3.exe
	/bin/rm -r ex3.o

ex4: ex4.cpp
	$(CXX)     ex4.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex4.o   $(MY_LINKER_FLAGS)    -o ex4.exe
	/bin/rm -r ex4.o

ex5: ex5.cpp
	$(CXX)     ex5.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex5.o   $(MY_LINKER_FLAGS)    -o ex5.exe
	/bin/rm -r ex5.o

ex6: ex6.cpp
	$(CXX)     ex6.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex6.o   $(MY_LINKER_FLAGS)    -o ex6.exe
	/bin/rm -r ex6.o

ex7: ex7.cpp
	$(CXX)     ex7.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex7.o   $(MY_LINKER_FLAGS)    -o ex7.exe
	/bin/rm -r ex7.o

ex8: ex8.cpp
	$(CXX)     ex8.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex8.o   $(MY_LINKER_FLAGS)    -o ex8.exe
	/bin/rm -r ex8.o

ex9: ex9.cpp
	$(CXX)     ex9.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex9.o   $(MY_LINKER_FLAGS)    -o ex9.exe
	/bin/rm -r ex9.o

ex10: ex10.cpp
	$(CXX)     ex10.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex10.o   $(MY_LINKER_FLAGS)    -o ex10.exe
	/bin/rm -r ex10.o

ex11: ex11.cpp
	$(CXX)     ex11.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex11.o   $(MY_LINKER_FLAGS)    -o ex11.exe
	/bin/rm -r ex11.o

ex12: ex12.cpp
	$(CXX)     ex12.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex12.o   $(MY_LINKER_FLAGS)    -o ex12.exe
	/bin/rm -r ex12.o

ex13: ex13.cpp
	$(CXX)     ex13.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex13.o   $(MY_LINKER_FLAGS)    -o ex13.exe
	/bin/rm -r ex13.o

ex14: ex14.cpp
	$(CXX)     ex14.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex14.o   $(MY_LINKER_FLAGS)    -o ex14.exe
	/bin/rm -r ex14.o

ex15: ex15.cpp
	$(CXX)     ex15.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex15.o   $(MY_LINKER_FLAGS)    -o ex15.exe
	/bin/rm -r ex15.o

ex16: ex16.cpp
	$(CXX)     ex16.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex16.o   $(MY_LINKER_FLAGS)    -o ex16.exe
	/bin/rm -r ex16.o

ex17: ex17.cpp
	$(CXX)     ex17.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex17.o   $(MY_LINKER_FLAGS)    -o ex17.exe
	/bin/rm -r ex17.o

ex18: ex18.cpp
	$(CXX)     ex18.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex18.o   $(MY_LINKER_FLAGS)    -o ex18.exe
	/bin/rm -r ex18.o

ex19: ex19.cpp
	$(CXX)     ex19.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex19.o   $(MY_LINKER_FLAGS)    -o ex19.exe
	/bin/rm -r ex19.o

ex20: ex20.cpp
	$(CXX)     ex20.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex20.o   $(MY_LINKER_FLAGS)    -o ex20.exe
	/bin/rm -r ex20.o

ex21: ex21.cpp
	$(CXX)     ex21.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex21.o   $(MY_LINKER_FLAGS)    -o ex21.exe
	/bin/rm -r ex21.o

ex22: ex22.cpp
	$(CXX)     ex22.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex22.o   $(MY_LINKER_FLAGS)    -o ex22.exe
	/bin/rm -r ex22.o

ex23: ex23.cpp
	$(CXX)     ex23.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex23.o   $(MY_LINKER_FLAGS)    -o ex23.exe
	/bin/rm -r ex23.o

ex24: ex24.cpp
	$(CXX)     ex24.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex24.o   $(MY_LINKER_FLAGS)    -o ex24.exe
	/bin/rm -r ex24.o

ex25: ex25.cpp
	$(CXX)     ex25.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex25.o   $(MY_LINKER_FLAGS)    -o ex25.exe
	/bin/rm -r ex25.o

clean:
	/bin/rm -f *.o ex?.exe ex??.exe
