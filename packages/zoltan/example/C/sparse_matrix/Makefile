##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
##############################################################################
# CVS File Information
#    $RCSfile$
#    $Author$
#    $Date$
#    Revision: 1.107 $
##############################################################################

##############################################################################
CC = mpicc 
CFLAGS = 

# Link with Zoltan library 
LIBS = -L../Obj_linux -lzoltan 

# Include ParMetis libraries if necessary 
LIBS += -L/Net/proj/zoltan/arch/linux/lib/ParMETIS3 -lparmetis -lmetis
                         
# Path to zoltan.h
INC_DIRS = ../include

# All
all: ex1 

# Make .o files
ex1.o: ex1.c matrix.h
	$(CC) -c $(CFLAGS) -I$(INC_DIRS)  ex1.c
queries1.o: queries1.c matrix.h
	$(CC) -c $(CFLAGS) -I$(INC_DIRS) queries1.c

# Main program 
ex1: ex1.o queries1.o 
	$(CC) ex1.o queries1.o -o $@ $(LIBS) -lm 
