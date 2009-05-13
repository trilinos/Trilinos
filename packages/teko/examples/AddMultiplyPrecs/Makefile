include ../../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE) -I/home/eccyr/Utils/EcUtils++/include
LDDIR = -L$(TRILINOS_LIB) -L$(PB_LIBDIR) -L/home/eccyr/Utils/EcUtils++/lib
LDLIB = -lpb -lecutils++ $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = 

all: addprecs multprecs

multprecs: Driver.cpp $(OBJS) $(PB_LIBDIR)/libpb.a
	$(CC) $(INCLUDE) $(CPPFLAGS) $(OBJS) $(LDDIR) Driver.cpp -o multprecs $(LDLIB) 

addprecs: Driver.cpp $(OBJS) $(PB_LIBDIR)/libpb.a 
	$(CC) -DADD_PREC $(INCLUDE) $(CPPFLAGS) $(OBJS) $(LDDIR) Driver.cpp -o addprecs $(LDLIB) 

clean:
	-rm multprecs addprecs
