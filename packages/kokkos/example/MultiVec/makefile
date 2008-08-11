FC=g77
FFLAGS = -O5 -funroll-all-loops -malign-double
BLASF=dgemmf.o
TIMER=second_linux.o

TARGET1=tgemmV1
TARGET2=tgemmV2

LIBS = -lblas

all: $(TARGET1) $(TARGET2)

$(TARGET1): tgemmV1.o $(TIMER) $(BLASF)
	$(FC) $(FFLAGS) -o $(TARGET1) tgemmV1.o $(TIMER) $(BLASF) $(LIBS)

$(TARGET2): tgemmV2.o $(TIMER) $(BLASF)
	$(FC) $(FFLAGS) -o $(TARGET2) tgemmV2.o $(TIMER) $(BLASF) $(LIBS)
clean:
	rm -rf *.o
