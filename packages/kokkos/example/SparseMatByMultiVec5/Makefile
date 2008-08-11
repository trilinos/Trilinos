
FC=g77
#FFLAGS = -O5 -funroll-all-loops -malign-double
FFLAGS = -O5 -funroll-all-loops

#CFLAGS = -DWALL

LD = g77 -O5 -funroll-all-loops

LIBS = 

OBJS= tsploops.o csrcsc.o  scscrd.o scscmm1a.o scscmm2a.o scsrmm1a5.o scsrmm2a5.o \
      snrm2.o srelerr.o swaxpy.o timer.o

TARGET=tsploops

$(TARGET): $(OBJS)
	$(LD) -o $(TARGET)  $(LIBS) $(OBJS)

clean:
	rm -rf $(TARGET) $(OBJS)
