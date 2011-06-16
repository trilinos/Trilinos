# Define the dependency files
DEPDIR   = .deps
DEPFILES = $(shell ls $(DEPDIR)/*.Pi 2> /dev/null; ls $(DEPDIR)/dummy)

# Define the list of packages that are supported
PACKAGES  = Teuchos Epetra TriUtils EpetraExt Isorropia AztecOO \
	    Galeri Amesos IFPACK Komplex Anasazi ML NOX

# Define the list of swig output files to be generated, the generated
# dependency files, and the generated XML files
SWIGFILES   = $(patsubst %, %_dox.i, $(PACKAGES))
GENDEPFILES = $(patsubst %, $(DEPDIR)/%_dox.Pi, $(PACKAGES))
GENXMLFILES = $(patsubst %, %/*.x*, $(PACKAGES))

# Default target
.PHONY: all
all: $(SWIGFILES)

# Include the dependency files
include $(DEPFILES)
.PHONY: depend
depend: $(GENDEPFILES)

# Debugging help
.PHONY: check
check:
	@echo DEPFILES    = $(DEPFILES)
	@echo PACKAGES    = $(PACKAGES)
	@echo SWIGFILES   = $(SWIGFILES)
	@echo GENDEPFILES = $(GENDEPFILES)
#	@echo GENXMLFILES = $(GENXMLFILES)

# Generation rule for swig output files
%_dox.i: Doxyfile_%
	doxygen $^ 2> /dev/null; \
	./doxy2swig.py $*/index.xml $@;

# Generate the dependency files
.SECONDARY: $(GENDEPFILES)
$(DEPDIR)/%_dox.Pi:
	./makedeps.py $*_dox.i Doxyfile_$* > $@

# Clean the directory
.PHONY: clean
clean:
	$(RM) -r $(PACKAGES)
	$(RM) $(DEPDIR)/*.Pi

# Clobber the directory, including generated *_dox.i files
.PHONY: clobber
clobber: clean
	$(RM) *_dox.i
