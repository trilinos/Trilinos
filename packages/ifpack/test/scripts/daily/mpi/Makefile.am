#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

noinst_PROGRAMS = ifpackTestAllMpi

ifpackTestAllMpi_SOURCES = $(srcdir)/ifpackTestAllMpi

ifpackTestAllMpi$(EXEEXT): $(srcdir)/ifpackTestAllMpi
	cp $(srcdir)/ifpackTestAllMpi ifpackTestAllMpi$(EXEEXT) 
else

EXTRA_DIST = ifpackTestAllMpi

endif

