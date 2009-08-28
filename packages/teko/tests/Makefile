.PHONY: clean subdirs all preproc

TRILINOS_DIR=/usr/local/trilinos-9/

include ../make.include

INCLUDE = -I$(PB_INCLUDE) -I$(PB_TESTS_INCLUDE) -I$(TRILINOS_INCLUDE) 
LDDIR = -L$(TRILINOS_LIB) -L$(PB_LIBDIR)
LDLIB = -lpbtests -lpb $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = Test_Utils.o
PREC = $(OBJS:%.o=%._tmp)

DRIVERS = testdriver

all: subdirs $(OBJS) libpbtests.a $(DRIVERS)

preproc: $(PREC)

subdirs: $(OBJS)
	make -C src/
	ar rcs libpbtests.a $(PB_DIR)/tests/objs/*.o
	cp libpbtests.a $(PB_DIR)/lib/libpbtests.a

exercise_Block2x2PrecFact: exercise_Block2x2PrecFact.cpp $(OBJS) $(PB_LIBDIR)/libpb.a
	$(CC) $(INCLUDE) $(CPPFLAGS) -L./ $(LDDIR) exercise_Block2x2PrecFact.cpp -o $@ $(LDLIB)

testdriver: testdriver.cpp $(OBJS) $(PB_LIBDIR)/libpb.a $(PB_LIBDIR)/libpbtests.a
	$(CC) $(INCLUDE) $(CPPFLAGS) -L./ $(LDDIR) testdriver.cpp -o $@ $(LDLIB)

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $<
	cp $@ $(PB_DIR)/tests/objs/$@

$(PREC): %._tmp : %.cpp %.hpp
	cpp $(INCLUDE) $<

libpbtests.a:
	ar rcs libpbtests.a $(PB_DIR)/tests/objs/*.o
	cp libpbtests.a $(PB_DIR)/lib/libpbtests.a

clean:
	-rm failure.log
	-make -C src/ clean
	-rm -f libpbtests.a $(DRIVERS) $(OBJS) $(PB_DIR)/lib/libpbtests.a $(OBJS:%=$(PB_DIR)/tests/objs/%)
	-rm -f objs/*.o
