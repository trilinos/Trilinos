#@HEADER
#@HEADER
# Scripts in this directory are automatically detected by the test harness
# if set up to build properly below.  Scripts are executed by the test harness.

#np# This Makefile.am basically copies the "autotoolsTestAllMpi" script from 
#np# the source tree to the build tree and gives it a .exe extension.  To set
#np# up multiple scripts, list all scripts in the "PROGRAMS" primary and create
#np# a (script_name)_SOURCES primary and make rule for each script.  Don't 
#np# forget to get rid of the "autotoolsTestAllMpi" example.  This 
#np# directory is for mpi scripts only.  Keep all rules and primaries within 
#np# the "if HAVE_MPI" statement.

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

noinst_PROGRAMS = autotoolsTestAllMpi

autotoolsTestAllMpi_SOURCES = $(srcdir)/autotoolsTestAllMpi

autotoolsTestAllMpi$(EXEEXT): $(srcdir)/autotoolsTestAllMpi
	cp $(srcdir)/autotoolsTestAllMpi autotoolsTestAllMpi$(EXEEXT) 

endif

