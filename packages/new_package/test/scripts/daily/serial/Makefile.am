#@HEADER
#@HEADER
# Scripts in this directory are automatically detected by the test harness
# if set up to build properly below.  Scripts are executed by the test harness.

#np# This Makefile.am basically copies the "autotoolsTestAll" script from the 
#np# source tree to the build tree and gives it a .exe extension.  To set up
#np# multiple scripts, list all scripts in the "PROGRAMS" primary and create
#np# a (script_name)_SOURCES primary and make rule for each script.  Don't 
#np# forget to get rid of the "autotoolsTestAll" example.  This  directory 
#np# is for serial scripts only.  Keep all rules and primaries within 
#np# the "else" portion of the "if HAVE_MPI" statement.

EXEEXT = .exe
# Nothing should be done for an mpi build
if HAVE_MPI
else

noinst_PROGRAMS = autotoolsTestAll


autotoolsTestAll_SOURCES = $(srcdir)/autotoolsTestAll

autotoolsTestAll$(EXEEXT): $(srcdir)/autotoolsTestAll
	cp $(srcdir)/autotoolsTestAll autotoolsTestAll$(EXEEXT) 

endif
