
include ../make.platform
include ../config/vars.$(FEI_ARCH)

default:all

all: test0.x test1.x test2.x

test0.x: test0.cpp
	@echo "config.makefile: CXX=$(CXX)"
	$(CXX) test0.cpp -o test0.x

test3.o: test3.cpp
	@echo "config.makefile: CXX=$(CXX)"
	@echo "config.makefile: CXXFLAGS=$(CXXFLAGS)"
	@echo "config.makefile: CPPFLAGS=$(CPPFLAGS)"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c test3.cpp -o test3.o

test1.o: test1.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c test1.cpp -o test1.o

test1.x: test1.o
	@echo "config.makefile: LDFLAGS=$(LDFLAGS)"
	@echo "config.makefile: LIBS=$(LIBS)"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test1.o $(LDFLAGS) -o test1.x $(LIBS)

test2.o: test2.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c test2.cpp -o test2.o

test2.x: test2.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test2.o $(LDFLAGS) -o test2.x $(LIBS)

test_exception.o: test_exception.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c test_exception.cpp -o test_exception.o

test_exception.x: test_exception.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test_exception.o $(LDFLAGS) -o test_exception.x $(LIBS)

