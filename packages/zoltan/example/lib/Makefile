# $Id$
#

EXAMPLE_ZOLTAN_LIB_OBJECTS = exmesh.o exmisc.o

all:	libexzoltan.a

libexzoltan.a: $(EXAMPLE_ZOLTAN_LIB_OBJECTS)
	ar r libexzoltan.a $(EXAMPLE_ZOLTAN_LIB_OBJECTS)

%.o : %.c exzoltan.h
	$(CC) $(INCLUDE) -c $<

clean:
	rm -f *.o libexzoltan.a 


