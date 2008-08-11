# $Id$
#

EXAMPLE_ZOLTAN_LIB_OBJECTS = exmesh.o exmisc.o exphg.o

all:	libexzoltan.a

libexzoltan.a: $(EXAMPLE_ZOLTAN_LIB_OBJECTS)
	@echo "Building $@"
	$(AR) $@ $(EXAMPLE_ZOLTAN_LIB_OBJECTS)

%.o : %.c exzoltan.h
	@echo "Building $@"
	@$(CC) $(CFLAGS) $(INCLUDE) -c $<

clean:
	rm -f *.o libexzoltan.a 


