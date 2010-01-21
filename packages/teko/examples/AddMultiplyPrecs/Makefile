include ../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(Teko_INCLUDE) -I/home/eccyr/Utils/EcUtils++/include
LDDIR = -L$(TRILINOS_LIB) -L$(Teko_LIBDIR) -L/home/eccyr/Utils/EcUtils++/lib
LDLIB = -lteko -lecutils++ $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = 

all: addprecs multprecs

multprecs: Driver.cpp $(OBJS) $(Teko_LIBDIR)/libteko.a
	$(CC) $(INCLUDE) $(CPPFLAGS) $(OBJS) $(LDDIR) Driver.cpp -o multprecs $(LDLIB) 

addprecs: Driver.cpp $(OBJS) $(Teko_LIBDIR)/libteko.a 
	$(CC) -DADD_PREC $(INCLUDE) $(CPPFLAGS) $(OBJS) $(LDDIR) Driver.cpp -o addprecs $(LDLIB) 

clean:
	-rm multprecs addprecs
