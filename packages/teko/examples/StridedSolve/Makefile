include ~/PB/make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE) -I/home/eccyr/Utils/EcUtils++/include
LDDIR = -L$(TRILINOS_LIB) -L$(PB_LIBDIR) -L/home/eccyr/Utils/EcUtils++/lib
LDLIB = -lpb -lecutils++ $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = 

all: strided_solve

strided_solve: strided_solve.cpp $(OBJS) $(PB_LIBDIR)/libpb.a ~/PB/make.include
	$(CC) $(INCLUDE) $(CPPFLAGS) $(OBJS) $(LDDIR) strided_solve.cpp -o strided_solve $(LDLIB) 

clean:
	-rm strided_solve
