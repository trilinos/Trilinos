#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for an mpi build
if HAVE_MPI
else

noinst_PROGRAMS = autotoolsTestAll

autotoolsTestAll_SOURCES = $(srcdir)/autotoolsTestAll

autotoolsTestAll$(EXEEXT): $(srcdir)/autotoolsTestAll
	cp $(srcdir)/autotoolsTestAll autotoolsTestAll$(EXEEXT) 

endif
