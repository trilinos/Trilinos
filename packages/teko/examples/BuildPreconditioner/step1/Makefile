include ../../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE)
LDDIR = -L$(TRILINOS_LIB) -L$(PB_LIBDIR)
LDLIB = -lpb $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = ExamplePreconditionerFactory.o

all: $(OBJS) example-test

example-test: example-test.cpp
	$(CC) $(INCLUDE) $(CPPFLAGS) example-test.cpp $(LDDIR) -o $@ $(LDLIB)

$(OBJS): %.o : %.cpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@

clean:
	-rm $(OBJS)
