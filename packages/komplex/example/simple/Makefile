TRILINOS_INSTALLATION=..
PUC=KOMPLEX
PLC=komplex
ROOT=cxx_main
 
include $(TRILINOS_INSTALLATION)/include/Makefile.export.$(PLC).macros
include $(TRILINOS_INSTALLATION)/include/Makefile.export.$(PLC)
include $(TRILINOS_INSTALLATION)/include/Makefile.export.epetraext.macros
include $(TRILINOS_INSTALLATION)/include/Makefile.export.epetraext
include $(TRILINOS_INSTALLATION)/include/Makefile.export.galeri.macros
include $(TRILINOS_INSTALLATION)/include/Makefile.export.galeri
 
INCDIR=$($(PUC)_INCLUDES) $(GALERI_INCLUDES)
LIBDIR= $(GALERI_LIBS) $($(PUC)_LIBS) 
 
CXXFLAGS=$(INCDIR) $($(PLC)_CXXFLAGS) $($(PLC)_DEFS) $($(PUC)_CXXFLAGS) $($(PUC)_DEFS)
CXX=$($(PLC)_CXX) $($(PUC)_CXX)
 
TARGET=$(ROOT).exe
SRC=$(ROOT).cpp
 
OBJ=$(SRC:.cpp=.o)
 
$(TARGET): $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ) $(LIBDIR)
 
clean:
	rm -rf $(OBJ) $(TARGET) *~
 
again: clean $(TARGET)
