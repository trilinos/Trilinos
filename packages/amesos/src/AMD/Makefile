#------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#------------------------------------------------------------------------------

# Compile everything, including the C-callable routine and the mexFunctions.
# Do not compile the FORTRAN versions.
all:
	( cd Source ; make )
	( cd Demo   ; make )
	( cd MATLAB ; make )
	# - cat Doc/License

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd Source ; make )
	( cd Demo   ; make )
	# - cat Doc/License

# compile the FORTRAN libraries and demo programs (not compiled by "make all")
fortran:
	( cd Source ; make fortran )
	( cd Demo   ; make fortran )
	# - cat Doc/License

# compile a FORTRAN demo program that calls the C version of AMD
# (not compiled by "make all")
cross:
	( cd Demo   ; make cross )
	# - cat Doc/License

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Source ; make clean )
	( cd Demo   ; make clean )
	( cd MATLAB ; make clean )
	( cd Doc    ; make clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Source ; make purge )
	( cd Demo   ; make purge )
	( cd MATLAB ; make purge )
	( cd Doc    ; make purge )

# create PDF documents for the original distribution
doc:
	( cd Doc    ; make )

# get ready for distribution
dist: purge
	( cd Demo   ; make dist )
	( cd Doc    ; make )

ccode: lib
