# Please add "export TRILINOS_INSTALL_DIR=/home/login/path-to-your-trilinos-install-dir"
# somewhere in you local configuration or use the following line but in this case, don't commit your Makefile !
# TRILINOS_INSTALL_DIR = 

EXE_ROOT = 	Cthulhu_Example_Comm

CLIENT_EXTRA_INCLUDES = 

CLIENT_EXTRA_LIBS =

#
# Include the Trilinos export makefile from package=Tpetra
#
# Note that there seems to be a bug in the export makefile system.  Sometimes
# an extra "-l" appears in package foo's macro foo_TPL_LIBRARIES.
# Reported in Trilinos bug #4946.
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Tpetra

LIB_PATH = $(TRILINOS_INSTALL_DIR)/lib

INCLUDE_PATH = $(TRILINOS_INSTALL_DIR)/include -I./ -I../src/Def -I../src/Platform -I../src/Map -I../src/Comm

CXX = $(Tpetra_CXX_COMPILER)

CXXFLAGS = $(Tpetra_CXX_FLAGS)

LIBS = $(CLIENT_EXTRA_LIBS) $(SHARED_LIB_RPATH_COMMAND)  $(Tpetra_LIBRARIES)  $(Tpetra_TPL_LIBRARIES) $(Tpetra_EXTRA_LD_FLAGS) 

EXTRA_OBJS = ../src/Platform/Cthulhu_DefaultPlatform.o

#
# Rules for building executables and objects. 
#

all: ${EXE_ROOT}.exe
	@echo "done"

%.exe : %.o $(EXTRA_OBJS) $(CLIENT_EXTRA_LIBS)
	$(CXX) -o $@ $(LDFLAGS) $(CXXFLAGS) -g -O0 -Wall $< $(EXTRA_OBJS) -L$(LIB_PATH) $(LIBS)

%.o : %.cpp ${CLIENT_EXTRA_INCLUDES}
	$(CXX) -c -o $@ $(CXXFLAGS) -g -O0 -Wall -I$(INCLUDE_PATH) $(Tpetra_TPL_INCLUDES) ${CLIENT_EXTRA_INCLUDES} $<

matrix_driver.o: matrix_driver.cpp matrix_factory.hpp matrix_types.hpp
	$(CXX) -c -o $@ $(CXXFLAGS) -g -O0 -I$(INCLUDE_PATH) $(Tpetra_TPL_INCLUDES) $<

clean:
	rm -f ${EXE_ROOT}.exe ${EXE_ROOT}.o
