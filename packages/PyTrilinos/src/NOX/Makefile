# Root PyTrilinos directory
ROOT       := ../../

# The interface files
INTERFACES := Abstract.i Parameter.i Solver.i StatusTest.i

# Common Makefile macros
include $(ROOT)common.mk

# Phony targets
.PHONY: all clean clobber

# Secondary targets
.SECONDARY: $(WRAPPERS) $(DEPEND)

# First target: default make
all: $(WRAPPERS)

# Install target
install: $(DYNAMIC_LIBS) $(PYC_FILES)
	install -d $(INSTALL_DIR)
	install $^ $(INSTALL_DIR)

# Remove object and wrapper files
clean:
	rm -f *.o *_wrap.cxx

# Remove dependency files, proxies, static and dynamic libraries,
# byte-compiled python files and documentation files
clobber: clean
	rm -f $(DEPDIR)*.d $(PROXIES) *.pyc
