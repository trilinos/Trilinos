#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

noinst_PROGRAMS = epetraTestAllMpi

epetraTestAllMpi_SOURCES = $(srcdir)/epetraTestAllMpi

epetraTestAllMpi$(EXEEXT): $(srcdir)/epetraTestAllMpi
	cp $(srcdir)/epetraTestAllMpi epetraTestAllMpi$(EXEEXT) 

endif

