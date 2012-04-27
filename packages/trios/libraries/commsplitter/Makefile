#CC = cc
#FC = ftn
CC = mpicc
FC = mpif90
AR = ar
RM = rm -f
CFLAGS = -O2
OBJS = init.o logging.o wrappers.o
LIB = libcommsplitter.a
TEST = ctest.exe ftest.exe

CFLAGS = -O0 -g

all: $(LIB) $(TEST)

$(LIB): $(OBJS)
	$(AR) rcs $(LIB) $(OBJS)

test: $(TEST)
	@echo "******************************************************************"
	@echo "Each app should report that it belongs to a communicator of size 2"
	@echo "******************************************************************"
	mpirun -n 2 ./ctest.exe : -n 2 ./ftest.exe

debug: all

$(TEST): $(OBJS) ftest.f90 ctest.c
	$(FC) -g ftest.f90 -o ftest.exe -L. -lcommsplitter -lz
	$(CC) -g ctest.c -o ctest.exe -L. -lcommsplitter -lz

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(LIB) $(TEST) *~
