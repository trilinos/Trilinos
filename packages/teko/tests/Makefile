.PHONY: clean subdirs all preproc

TRILINOS_DIR=/usr/local/trilinos-9/

include ../make.include

INCLUDE = -I$(Teko_INCLUDE) -I$(Teko_TESTS_INCLUDE) -I$(TRILINOS_INCLUDE) 
LDDIR = -L$(TRILINOS_LIB) -L$(Teko_LIBDIR) $(LIBDIR)
LDLIB = -ltekotests -lteko $(TRILINOS_LIBRARIES)
CPPFLAGS=-DHAVE_CONFIG_H $(THYRA_CXXFLAGS) $(DEBUG_FLAGS)

OBJS = Test_Utils.o
PREC = $(OBJS:%.o=%._tmp)

DRIVERS = testdriver

all: subdirs $(OBJS) libtekotests.a $(DRIVERS)

preproc: $(PREC)

subdirs: $(OBJS)
	mkdir -p objs
	make -C src/
	ar rcs libtekotests.a $(Teko_DIR)/tests/objs/*.o
	cp libtekotests.a $(Teko_DIR)/lib/libtekotests.a

exercise_Block2x2PrecFact: exercise_Block2x2PrecFact.cpp $(OBJS) $(Teko_LIBDIR)/libteko.a
	$(CC) $(INCLUDE) $(CPPFLAGS) -L./ $(LDDIR) exercise_Block2x2PrecFact.cpp -o $@ $(LDLIB)

testdriver: testdriver.cpp $(OBJS) $(Teko_LIBDIR)/libteko.a $(Teko_LIBDIR)/libtekotests.a
	$(CC) $(INCLUDE) $(CPPFLAGS) -L./ $(LDDIR) testdriver.cpp -o $@ $(LDLIB)

$(OBJS): %.o : %.cpp %.hpp
	$(CC) -c $(INCLUDE) $(CPPFLAGS) $<
	cp $@ $(Teko_DIR)/tests/objs/$@

$(PREC): %._tmp : %.cpp %.hpp
	cpp $(INCLUDE) $<

libtekotests.a:
	ar rcs libtekotests.a $(Teko_DIR)/tests/objs/*.o
	cp libtekotests.a $(Teko_DIR)/lib/libtekotests.a

clean:
	-rm failure.log
	-make -C src/ clean
	-rm -f libtekotests.a $(DRIVERS) $(OBJS) $(Teko_DIR)/lib/libtekotests.a $(OBJS:%=$(Teko_DIR)/tests/objs/%)
	-rm -f objs/*.o
