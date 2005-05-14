#-------------------------------------------------------------------------------
# Makefile for the CColamd directory
#-------------------------------------------------------------------------------

# TODO: don't use this,  use Pk/Makefile instead.
#
include ../Make/Make.include

#-------------------------------------------------------------------------------

OF = -output

LIB = ../Utils/utils.a $(LIBS)

ALLMEX = ccolamd.$(MEXT)

INCLUDE = ccolamd_internal.h ccolamd.h

OBJ = ccolamd_postorder.o ccolamd_post_tree.o ccolamd_apply_order.o ccolamd.o \
	ccolamd_fsize.o

$(OBJ): $(INCLUDE)

.c.o:
	$(CC) $(ORDER_CFLAGS) -c $*.c

matlab: $(ALLMEX) ccolamd.a

ccode: ccolamd.a

clean:
	'rm' -f *.o *.ln

purge: clean
	'rm' *.mex* ccolamd.a

ccolamd.a: $(OBJ)
	ar cr ccolamd.a $(OBJ)
	ranlib ccolamd.a
