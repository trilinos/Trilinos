all: exo2adams coords2rbm rigid_spread capt2_read sym2full az2bin merge2bin

CFLAGS=-I../src -I../src/Include
exo2adams: exo2adams.o
	cc -o exo2adams exo2adams.c
ml_rbm.o: ../src/Utils/ml_rbm.c
	cc -c ../src/Utils/ml_rbm.c
coords2rbm:  ml_coords2rbm.o ml_rbm.o
	cc -o coords2rbm ml_coords2rbm.o ml_rbm.o
rigid_spread: rigid_spread.c
	cc -o rigid_spread rigid_spread.c
capt2_read: az_capt2read.c
	cc -Dadams -o capt2_read az_capt2read.c
sym2full: sym2full.c
	cc -o sym2full sym2full.c
az2bin: az2bin.c
	cc -Dadams -o az2bin az2bin.c
merge2bin: merge2bin.c
	cc -Dadams -Dbinary -o merge2bin merge2bin.c
clean:
	rm -f *.o az2bin exo2adams coords2rbm rigid_spread capt2_read sym2full \
	          merge2bin
