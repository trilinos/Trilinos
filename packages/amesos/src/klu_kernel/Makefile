#-------------------------------------------------------------------------------
# Makefile for klu kernel
#-------------------------------------------------------------------------------

all: klu.mexglx klu_u.mexglx klu_l.mexglx klu_p.mexglx klu_ut.mexglx \
    klu_lt.mexglx

CC = mex -inline

INC = klu.h klu_kernel.h tictoc.h

klu_kernel.o: klu_kernel.c $(INC)
	mex -c klu_kernel.c

klu.o: klu.c $(INC)
	mex -c klu.c

tictoc.o: tictoc.c $(INC)
	mex -c tictoc.c

klu.mexglx: klumex.c klu.o tictoc.o klu_kernel.o $(INC)
	mex -output klu klu.o klu_kernel.o klumex.c tictoc.o

klu_u.mexglx: klu_umex.c klu.o tictoc.o klu_kernel.o $(INC)
	mex -output klu_u klu.o klu_kernel.o klu_umex.c tictoc.o

klu_l.mexglx: klu_lmex.c klu.o tictoc.o klu_kernel.o $(INC)
	mex -output klu_l klu.o klu_kernel.o klu_lmex.c tictoc.o

klu_ut.mexglx: klu_utmex.c klu.o tictoc.o klu_kernel.o $(INC)
	mex -output klu_ut klu.o klu_kernel.o klu_utmex.c tictoc.o

klu_lt.mexglx: klu_ltmex.c klu.o tictoc.o klu_kernel.o $(INC)
	mex -output klu_lt klu.o klu_kernel.o klu_ltmex.c tictoc.o

klu_p.mexglx: klu_pmex.c klu.o tictoc.o klu_kernel.o $(INC)
	mex -output klu_p klu.o klu_kernel.o klu_pmex.c tictoc.o

clean:
	'rm' -f *.o

purge: clean
	'rm' *.mex*

.c.o:
	$(CC) -c $*.c

