# $Id$
#

EXAMPLE_ZOLTAN_LIB_OBJECTS = exmesh.o exmisc.o

all:	libexzoltan.a

libexzoltan.a: $(EXAMPLE_ZOLTAN_LIB_OBJECTS)
	$(AR) libexzoltan.a $(EXAMPLE_ZOLTAN_LIB_OBJECTS)

%.o : %.c exzoltan.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

clean:
	rm -f *.o libexzoltan.a 


