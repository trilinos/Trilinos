#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

noinst_PROGRAMS = autotoolsTestAllMpi

autotoolsTestAllMpi_SOURCES = $(srcdir)/autotoolsTestAllMpi

autotoolsTestAllMpi$(EXEEXT): $(srcdir)/autotoolsTestAllMpi
	cp $(srcdir)/autotoolsTestAllMpi autotoolsTestAllMpi$(EXEEXT) 

endif

