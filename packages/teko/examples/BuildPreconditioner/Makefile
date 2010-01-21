include ../../make.include

INCLUDE = -I$(TRILINOS_INCLUDE) -I$(Teko_INCLUDE)
LDDIR = -L$(TRILINOS_LIB) -L$(Teko_LIBDIR) 
LDLIB = -lteko $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(DEBUG_FLAGS)

OBJS = 

all: $(OBJS) example-driver

example-driver: example-driver.cpp $(Teko_LIBDIR)/libteko.a
	$(CC) $(INCLUDE) $(CPPFLAGS) example-driver.cpp $(LDDIR) -o $@ $(LDLIB)

$(OBJS): %.o : %.cpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $< -o $@

clean:
	-rm $(OBJS) example-driver
