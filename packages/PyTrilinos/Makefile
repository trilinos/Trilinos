# Root PyTrilinos directory
ROOT       := ./

# Tarball name
TARBALL    := PyTrilinos.tgz

# Source directory
SRC_DIR    := src/

# Sub-directories
SUBDIRS    := $(SRC_DIR) $(SRC_DIR)NOX

# Common Makefile macros
include $(ROOT)common.mk

# Phony targets
.PHONY: all build clean clobber sweep all_clean all_clobber all_sweep

# First target: default make
all: build

#build: $(WRAPPERS) Common/libcommon.a NOX_Module
build: wrappers
	./setup.py build

wrappers:
	cd $(SRC_DIR) && make

install: build
	./setup.py install

# Remove object and wrapper files
clean:
	./setup.py clean

# Remove dependency files, proxies, static and dynamic libraries,
# byte-compiled python files and documentation files
clobber: clean
	rm -rf build

#
# Additional cleanup targets
all_clean: clean
	@for d in $(SUBDIRS); \
	do                    \
            echo $$d:;        \
	    cd $$d;           \
	    $(MAKE) clean;    \
	    cd ..;            \
	done

all_clobber: clobber
	@for d in $(SUBDIRS); \
	do                    \
            echo $$d:;        \
	    cd $$d;           \
	    $(MAKE) clobber;  \
	    cd ..;            \
	done

all_sweep: sweep
	@for d in $(SUBDIRS); \
	do                    \
            echo $$d:;        \
	    cd $$d;           \
	    $(MAKE) sweep;    \
	    cd ..;            \
	done

# Tarball target
tarball: all_sweep
	cd ..; \
	tar cvzf $(TARBALL) $(BASE_PWD)
	@echo "Tarball in ../$(TARBALL)"

# Target for listing all the tests
list_tests:
	grep "test:" test/Makefile

test: all
	cd test; $(MAKE) $@

test%: all
	cd test; $(MAKE) $@
