#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

if HAVE_AMESOS_SLUD
SLUDIST_TEST = TestSuperLUdist 
else
SLUDIST_TEST = 
endif

if HAVE_AMESOS_DSCPACK
DSCPACK_TEST = TestDscpack
else
DSCPACK_TEST = 
endif


noinst_PROGRAMS = $(SLUDIST_TEST) $(DSCPACK_TEST)

if HAVE_AMESOS_SLUD
TestSuperLUdist_SOURCES = $(srcdir)/TestSuperLUdist

TestSuperLUdist$(EXEEXT): $(srcdir)/TestSuperLUdist
	cp $(srcdir)/TestSuperLUdist TestSuperLUdist$(EXEEXT) 
endif

if HAVE_AMESOS_DSCPACK
TestDscpack_SOURCES = $(srcdir)/TestDscpack

TestDscpack$(EXEEXT): $(srcdir)/TestDscpack
	cp $(srcdir)/TestDscpack TestDscpack$(EXEEXT) 
endif

endif

