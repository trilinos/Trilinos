#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

noinst_PROGRAMS = TestSuperLUdist

TestSuperLUdist_SOURCES = $(srcdir)/TestSuperLUdist

TestSuperLUdist$(EXEEXT): $(srcdir)/TestSuperLUdist
	cp $(srcdir)/TestSuperLUdist TestSuperLUdist$(EXEEXT) 

endif

