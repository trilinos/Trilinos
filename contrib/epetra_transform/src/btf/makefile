INCLUDE	=
OBJS =	concmp.o genbtf.o corsum.o finsum.o fnrsum.o izero.o \
	maxmatch.o rectblk.o mmc13e.o prtivs.o strchk.o main.o mattrans.o
LDFLAGS	= 
LIBS = 

btf: $(OBJS) $(INCLUDE)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

clean:
	/bin/rm -f core *.o *.out
