## Set the Trilinos install directory 
TRILINOS_INSTALL_DIR = /home/jhu/Trilinos/development-clean/sandbox

## Include any direct Trilinos library dependencies - in this case only ML
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.ML

COMPILE_FLAGS = $(ML_CXX_FLAGS) $(ML_INCLUDE_DIRS)
# TPL => "Third Party Library"
LINK_FLAGS = $(ML_LIBRARY_DIRS) $(ML_LIBRARIES) $(ML_TPL_LIBRARY_DIRS) $(ML_TPL_LIBRARIES) $(ML_EXTRA_LD_FLAGS)

## ## Build your application code ## 
ml_read_MatrixMarket.cpp.exe: ml_read_MatrixMarket.o
	$(ML_CXX_COMPILER) $(ML_CXX_FLAGS) -o ml_read_MatrixMarket.exe ml_read_MatrixMarket.o $(LINK_FLAGS) 

ml_read_MatrixMarket.o: ml_read_MatrixMarket.cpp TumiStuff.h
	$(ML_CXX_COMPILER) $(COMPILE_FLAGS) -c ml_read_MatrixMarket.cpp

clean:
	\rm -f *.o *.exe *~
