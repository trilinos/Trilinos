#@HEADER
#@HEADER

EXEEXT = .exe

if HAVE_AMESOS_SLUS
SLUS_TEST = TestAmesosSerial
else
SLUS_TEST = 
endif

if HAVE_AMESOS_UMFPACK
UMFPACK_TEST = TestUmfpackSerial
else
UMFPACK_TEST = 
endif

noinst_PROGRAMS = $(SLUS_TEST) $(UMFPACK_TEST)

if HAVE_AMESOS_SLUS
TestAmesosSerial_SOURCES = $(srcdir)/TestAmesosSerial

TestAmesosSerial$(EXEEXT): $(srcdir)/TestAmesosSerial
	cp $(srcdir)/TestAmesosSerial TestAmesosSerial$(EXEEXT) 
endif

if HAVE_AMESOS_UMFPACK
TestUmfpackSerial_SOURCES = $(srcdir)/TestUmfpackSerial

TestUmfpackSerial$(EXEEXT): $(srcdir)/TestUmfpackSerial
	cp $(srcdir)/TestUmfpackSerial TestUmfpackSerial$(EXEEXT) 
endif


