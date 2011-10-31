CC = cc
FC = ftn
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

debug: all

$(TEST): $(OBJS)
	$(FC) -g ftest.f -o ftest.exe -L. -lcommsplitter -lz
	$(CC) -g ctest.c -o ctest.exe -L. -lcommsplitter -lz

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(LIB) $(TEST) *~
