TEST_CC = cc_main.cc
TEST_CC_OBJ          = $(TEST_CC:.cc=.o)

default: $(TEST_CC_OBJ)
	$(CXX) $(LDFLAGS) $(TEST_CC_OBJ)

clean:
	rm -rf *~ *.o a.out
