# Used for switching on memory checking
#C_pre          = purify -best-effort -cache-dir=/tmp

#Compiler
C_comp		= cc -DHTEST -DZOLTAN_NO_MPI

# Assume that environment OSTYPE is set to 'solaris' or 'linux'!
 OSTYPE        = solaris
# OSTYPE        = linux
C_flag_solaris  = -fast
#C_flag_linux    = -g -Wall

# Use to switch off time measurement
DTIME		= -DWITHTIME

VPATH           = ../Utilities/Memory
I		= -I../Utilities/Memory
C		= $(C_pre) $(C_comp) $(C_flag_$(OSTYPE)) $(DTIME) $(I)
L		= -L. -lhg -lm

src.c		= hg_test.c hg_readfile.c hg_util.c \
		  hg_part.c hg_match.c hg_pack.c hg_group.c hg_scale.c \
		  hg_coarse.c hg_global.c hg_local.c hg_sort.c \
		  zoltan_hg_readfile.c mem.c

all:            	./hg_test
./hg_test:		./libhg.a hg_test.o
			$(C) hg_test.o $(L) -o $@
./libhg.a:		$(src.c:.c=.o)
			ar rv $@ $?
.c.o:			$<
	  		$(C) -c $<

clean:	
	  		rm -f ./hg_test ./libhg.a $(src.c:.c=.o)
