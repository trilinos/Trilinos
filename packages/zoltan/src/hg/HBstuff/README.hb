
 
                       Harwell-Boeing File I/O in C
                                V. 1.0
 
           National Institute of Standards and Technology, MD.
                             K.A. Remington
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                NOTICE
 
 Permission to use, copy, modify, and distribute this software and
 its documentation for any purpose and without fee is hereby granted
 provided that the above copyright notice appear in all copies and
 that both the copyright notice and this permission notice appear in
 supporting documentation.
 
 Neither the Author nor the Institution (National Institute of Standards
 and Technology) make any representations about the suitability of this
 software for any purpose. This software is provided "as is" without
 expressed or implied warranty.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-----------
DESCRIPTION
-----------

This package provides several I/O functions for Harwell-Boeing file,
and is particularly useful for C programs, since the Harwell-Boeing
format can rely heavily on Fortran specific formatting conventions.
A makefile and several driver programs demonstrating usage are also
provided.

------------------
NOTE ON INDEX BASE
------------------

The column pointer and row index vectors are by default 1-based (Fortran
compatible array indexing).  To read/write 0-based vectors (for C array
indexing compatibility), simply change the macro definition in the
makefile. (use -D_SP_base=0)

------------------
DISTRIBUTION FILES
-----------------

            Main source code:    iohb.c, iohb.h
 Matrix Market functionality:    mmio.c, mmio.h
        Makefile for drivers:    makefile
              Sample drivers:    hb2mtxstrm.c
                                 hbmat2hb.c
                                 hbmat2mtx.c
                                 hbrhs2mtx.c
                                 sample.c
                 Sample data:    data/*
                     
--------------
SAMPLE DRIVERS
--------------

To build the driver programs, edit the makefile to reflect your
compiler and desired compile-time flags. Then type 'make'.

The drivers included are:
sample     - Demonstrates basic usage, with many explanatory comments.
             The functionality is uninteresting (it reads a Harwell-Boeing
             file and writes out another Harwell-Boeing file), but it
             demonstrates the I/O capabilities typically required.

hbmat2hb   - Another example of usage, demonstrating how to read the    
             input data as characters so that output will exactly reflect 
             the accuracy of the input.  

hbmat2mtx  - Another example of usage, with more useful functionality.
             The program reads a Harwell-Boeing file and writes a commented
             Matrix Market formatted file reflecting the original data.

hbrhs2mtx  - This program will extract any available auxillary vectors
             from the Harwell-Boeing files (assuming they are in "full"
             form), and creates a Matrix Market array formatted file for
             each vector.  (The data is read as character input, so that 
             accuracy in the generated output is identical to that of 
             the input file.)

hb2mtxstrm - This program is similar to hbmat2mtx, but doesn't store
             the values of the matrix... instead it streams through the
             entries, writing out a Matrix Market formatted file along
             the way.  This saves considerably on translation time
             for large matrices.  This program is only intended for 
             translating files of REAL and COMPLEX data, not PATTERN files.
             PATTERN files can be translated with hbmat2mtx.
             (The data is read as character input, so that accuracy in
             the generated output is identical to that of the input file.)

-----------------
SAMPLE DATA FILES
-----------------

A data directory containing various Harwell-Boeing formatted files
is also provided.  The files in this directory reflect some of the
"legal", but very difficult to parse, incarnations of the Harwell-Boeing
format.  They can be used as arguments to the provided driver programs.

