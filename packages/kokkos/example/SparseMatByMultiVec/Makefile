FFLAGS = -O3
LD = f77 -O3
#FFLAGS = -g -O0
#LD = g77 -g -O0
#FC=g77
LIBS = 

OBJS= tsploops.o csrcsc.o  scscrd.o scscmm1.o scscmm2a.o scsrmm1.o scsrmm2a.o \
      snrm2.o srelerr.o swaxpy.o second_sparc.o

TARGET=tsploops

$(TARGET): $(OBJS)
	$(LD) -o $(TARGET)  $(LIBS) $(OBJS)

clean:
	rm -rf $(TARGET) $(OBJS)
