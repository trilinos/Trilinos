#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

if HAVE_AMESOS_SLUD
SLUDIST_TEST = TestOldSuperludist 
else
SLUDIST_TEST = 
endif

if HAVE_AMESOS_DSCPACK
DSCPACK_TEST = TestDscpack
else
DSCPACK_TEST = 
endif

if HAVE_AMESOS_UMFPACK
UMFPACK_TEST = TestUmfpack
else
UMFPACK_TEST = 
endif

if HAVE_AMESOS_KLU
KLU_TEST = TestKlu
else
KLU_TEST = 
endif

if HAVE_AMESOS_SUPERLUDIST
SUPERLUDIST_TEST = TestSuperlud
else
SUPERLUDIST_TEST = 
endif

if HAVE_AMESOS_MUMPS
MUMPS_TEST = TestMumps
else
MUMPS_TEST = 
endif


noinst_PROGRAMS = $(SLUDIST_TEST) $(DSCPACK_TEST) $(UMFPACK_TEST) \
	$(MUMPS_TEST) $(KLU_TEST) $(SUPERLUDIST_TEST) 

if HAVE_AMESOS_SLUD
TestOldSuperludist_SOURCES = $(srcdir)/TestOldSuperludist

TestOldSuperludist$(EXEEXT): $(srcdir)/TestOldSuperludist
	cp $(srcdir)/TestOldSuperludist TestOldSuperludist$(EXEEXT) 
endif

if HAVE_AMESOS_DSCPACK
TestDscpack_SOURCES = $(srcdir)/TestDscpack

TestDscpack$(EXEEXT): $(srcdir)/TestDscpack
	cp $(srcdir)/TestDscpack TestDscpack$(EXEEXT) 
endif

if HAVE_AMESOS_UMFPACK
TestUmfpack_SOURCES = $(srcdir)/TestUmfpack

TestUmfpack$(EXEEXT): $(srcdir)/TestUmfpack
	cp $(srcdir)/TestUmfpack TestUmfpack$(EXEEXT) 
endif

if HAVE_AMESOS_KLU
TestKlu_SOURCES = $(srcdir)/TestKlu

TestKlu$(EXEEXT): $(srcdir)/TestKlu
	cp $(srcdir)/TestKlu TestKlu$(EXEEXT) 
endif

if HAVE_AMESOS_SUPERLUDIST
TestSuperlud_SOURCES = $(srcdir)/TestSuperlud

TestSuperlud$(EXEEXT): $(srcdir)/TestSuperlud
	cp $(srcdir)/TestSuperlud TestSuperlud$(EXEEXT) 
endif

if HAVE_AMESOS_MUMPS
TestMumps_SOURCES = $(srcdir)/TestMumps

TestMumps$(EXEEXT): $(srcdir)/TestMumps
	cp $(srcdir)/TestMumps TestMumps$(EXEEXT) 
endif

endif

