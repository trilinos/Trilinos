CC = mpicc
CFLAGS = -Wall -Wstrict-prototypes -O3 -m64
# CFLAGS = -Wall -Wstrict-prototypes -g -m64 -DZOLTANV31
# CFLAGS = -Wall -Wstrict-prototypes -g -m64 -D_DEBUG

BUILD = linux

ZOLTAN31 = ../../../Zoltan_v3.1_Release
ZOLTAN32 = ../../../Zoltan_v3.2_buggy

INC31 = -I$(ZOLTAN31)/src/include
LIBP31 = $(ZOLTAN31)/src/Obj_$(BUILD)
LIBS31 = -L$(LIBP31) -lzoltan -lm
INC32 = -I$(ZOLTAN32)/src/include
LIBP32 = $(ZOLTAN32)/src/Obj_$(BUILD)
LIBS32 = -L$(LIBP32) -lzoltan -lm

all: zcol zcoldriver zcoldriver-v3.1

zcol : zcol.c $(LIBP31)/libzoltan.a makefile
	$(CC) $(CFLAGS) -DZOLTANV31 -o $@ zcol.c $(INC31) $(LIBS31)


zcoldriver-v3.1 : zcoldriver.c $(LIBP31)/libzoltan.a makefile
	$(CC) $(CFLAGS) -DZOLTANV31 -o $@ zcoldriver.c $(INC31) $(LIBS31)

zcoldriver : zcoldriver.c $(LIBP32)/libzoltan.a makefile
	$(CC) $(CFLAGS) -o $@ zcoldriver.c $(INC32) $(LIBS32)

clean:
	rm -f *.o  zcoldriver zcoldriver-v3.1 zcol



