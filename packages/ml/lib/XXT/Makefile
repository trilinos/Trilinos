OBJ= bit_mask.o bss_malloc.o comm.o error.o gs.o ivec.o queue.o stack.o xxt.o
CFLAGS = -g -Dr8 -I../../src  $(MPI_INC_PATH)

libxxt.a: $(OBJ)
	ar ruv libxxt.a *.o
	
clean:
	rm -f $(OBJ) *.a
