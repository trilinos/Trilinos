# Root PyTrilinos directory
ROOT       := ../

# Source Files
# SRC        := Callback.cxx            \
#               NumPyArray.cxx          \
#               NumPyWrapper.cxx        \
# 	      PyInterface.cxx

# The interface files
#INTERFACES := RawEpetra.i EpetraExt.i NOX_Epetra.i
INTERFACES := RawEpetra.i EpetraExt.i

# Sub-directories
SUBDIRS    := NOX

# Common Makefile macros
include $(ROOT)common.mk

# Default targets
TARGET     := $(WRAPPERS) NOX_Wrappers

# Phony targets
.PHONY: all clean clobber

# Secondary targets
.SECONDARY: $(DEPEND)

# First target: default make
all: $(TARGET)

# Target dependency on common library
$(TARGET): $(SRC)

# NOX Wrapper target
NOX_Wrappers:
	cd NOX && $(MAKE)

# Remove object and wrapper files
clean:
	rm -f *.o *_wrap.cxx

# Remove dependency files, proxies, static library and proxies
clobber: clean
	rm -f $(DEPDIR)*.d $(PROXIES)
