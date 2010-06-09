/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*****************************************************************************
*
* testcp - copy file test.exo created by testwt
*
* author - Sandia National Laboratories
*          Larry A. Schoof - Original
*
*          
* environment - UNIX
*
* entry conditions - 
*   input parameters:
*
* exit conditions - 
*
* revision history - 
*
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "netcdf.h"
#include "exodusII.h"


int main (int argc, char **argv)
{
   int exoid, exoid1, error, idum;
   int CPU_word_size,IO_word_size;

   float version;

   char *cdum = 0;

   ex_opts(EX_VERBOSE | EX_ABORT);

/* open EXODUS II files */

   CPU_word_size = 0;                   /* sizeof(float) */
   IO_word_size = 0;                    /* use size in file */

   exoid = ex_open ("test.exo",         /* filename path */
                     EX_READ,           /* access mode = READ */
                     &CPU_word_size,    /* CPU word size */
                     &IO_word_size,     /* IO word size */
                     &version);         /* ExodusII library version */

   printf ("\nafter ex_open\n");
   if (exoid < 0) exit(1);

   printf ("test.exo is an EXODUSII file; version %4.2f\n",
            version);
   printf ("         CPU word size %1d\n",CPU_word_size);
   printf ("         I/O word size %1d\n",IO_word_size);
   ex_inquire(exoid,EX_INQ_API_VERS, &idum, &version, cdum);
   printf ("EXODUSII API; version %4.2f\n", version);

   CPU_word_size = 8;                   /* this really shouldn't matter for
                                           the copy but tests the conversion
                                           routines */
   IO_word_size = 4;

   exoid1 = ex_create ("testcp.exo",    /* filename */
                        EX_CLOBBER,     /* OK to overwrite */
                        &CPU_word_size, /* CPU float word size in bytes */
                        &IO_word_size); /* I/O float word size in bytes */

   printf ("\nafter ex_create, exoid = %3d\n",exoid1);
   if (exoid1 < 0) exit(1);

   printf ("         CPU word size %1d\n",CPU_word_size);
   printf ("         I/O word size %1d\n",IO_word_size);

   /* ncopts = NC_VERBOSE; */

   error = ex_copy (exoid, exoid1);
   printf ("\nafter ex_copy, error = %3d\n", error);

   error = ex_close (exoid);
   printf ("\nafter ex_close, error = %3d\n", error);

   error = ex_close (exoid1);
   printf ("\nafter ex_close, error = %3d\n", error);
   return 0;
}
