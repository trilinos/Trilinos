#@HEADER
#@HEADER

EXEEXT = .exe

if HAVE_AMESOS_SLUS
SLUS_TEST = TestAmesosSerial
else
SLUS_TEST = 
endif

noinst_PROGRAMS = $(SLUS_TEST)

if HAVE_AMESOS_SLUS
TestAmesosSerial_SOURCES = $(srcdir)/TestAmesosSerial

TestAmesosSerial$(EXEEXT): $(srcdir)/TestAmesosSerial
	cp $(srcdir)/TestAmesosSerial TestAmesosSerial$(EXEEXT) 
endif


