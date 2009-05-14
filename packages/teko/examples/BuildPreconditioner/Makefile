include ../../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(PB_INCLUDE) -I/home/eccyr/Utils/EcUtils++/include
LDDIR = -L$(TRILINOS_LIB) -L$(PB_LIBDIR) -L/home/eccyr/Utils/EcUtils++/lib
LDLIB = -lpb $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = 

all: $(OBJS) example-driver

example-driver: example-driver.cpp
	$(CC) $(INCLUDE) $(CPPFLAGS) example-driver.cpp $(LDDIR) -o $@ $(LDLIB)

$(OBJS): %.o : %.cpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@

clean:
	-rm $(OBJS)
