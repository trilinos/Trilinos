OBJ= bit_mask.o bss_malloc.o comm.o error.o gs.o ivec.o queue.o stack.o xyt.o
CFLAGS = -fast -Knoieee -Dr8 -DMPISRC -DMLSRC -DTIMINGS -I../../Obj $(MPI_INC_PATH) 
#CFLAGS = -fast -Knoieee -Dr8 -DMPISRC -DMLSRC -DTIMINGS -DBLAS -I../../Obj $(MPI_INC_PATH) 
#CFLAGS = -fast -Knoieee -Dr8 -DMPISRC -DMLSRC -DML_TIMING -I../../Obj $(MPI_INC_PATH) 
#CFLAGS = -fast -Knoieee -DDEBUG -DVV -Dr8 -DMPISRC -DMLSRC -I../../Obj $(MPI_INC_PATH)
CC=pgcc -cougar

libxyt.a: $(OBJ)
	xar ruv libxyt.a *.o

clean:
	rm -f $(OBJ) *.a
