
FC=g77
FFLAGS = -O5 -funroll-all-loops -malign-double

#CFLAGS = -DWALL

LD = g77 -O5 -funroll-all-loops -malign-double

LIBS = 

OBJS= tsploops.o csrcsc.o  scscrd.o scscmm1.o scscmm2a.o scsrmm1.o scsrmm2a.o \
      snrm2.o srelerr.o swaxpy.o timer.o

TARGET=tsploops

$(TARGET): $(OBJS)
	$(LD) -o $(TARGET)  $(LIBS) $(OBJS)

clean:
	rm -rf $(TARGET) $(OBJS)
