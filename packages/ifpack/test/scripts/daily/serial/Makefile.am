#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for an mpi build
if HAVE_MPI
else

noinst_PROGRAMS = ifpackTestAll

ifpackTestAll_SOURCES = $(srcdir)/ifpackTestAll

ifpackTestAll$(EXEEXT): $(srcdir)/ifpackTestAll
	cp $(srcdir)/ifpackTestAll ifpackTestAll$(EXEEXT) 

endif
