include ../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(Teko_INCLUDE) -I/home/eccyr/Utils/EcUtils++/include
LDDIR = -L$(TRILINOS_LIB) -L$(Teko_LIBDIR) -L/home/eccyr/Utils/EcUtils++/lib
LDLIB = -lteko -lecutils++ $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = 

all: strided_solve

strided_solve: strided_solve.cpp $(OBJS) $(Teko_LIBDIR)/libteko.a
	$(CC) $(INCLUDE) $(CPPFLAGS) $(OBJS) $(LDDIR) strided_solve.cpp -o strided_solve $(LDLIB) 

clean:
	-rm strided_solve
