/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
#ifndef PS_PARIO_CONST_H
#define PS_PARIO_CONST_H

#include "rf_io_const.h"

/* Global variables. */
extern double   PIO_Time_Array[];		/* Vector for timings */

/*
 * The following variables are used when a single processor is to write
 * info for a processor other than itself.
 */


/* Structure used to store the information necessary for parallel I/O. */

struct Parallel_IO
{
  int     Dsk_List_Cnt;

  int    *Dsk_List;
  int   **RDsk_List;

  int     Num_Dsk_Ctrlrs;            /* The number of disk controllers.     */
  int     Num_Dsks_PCtrlr;           /* The number of disks per controller. */
  int     PDsk_Add_Fact;             /* The offset from zero used by the    */
                                     /* the target machine.                 */

  int     Zeros;        /* 1 - if the target machine uses leading zeros when */
                        /*     designating the disk number (eg - the paragon */
                        /*     uses /pfs/io_01)                              */
                        /* 0 - if it does not (eg - the tflop uses           */
                        /*     /pfs/tmp_1)                                   */

  /* 1 - don't create a subdirectory for the spread files; write
   * them in same directory as the mesh file
   */
  int     NoSubdirectory;

  /* The root location of the parallel disks */
  char    Par_Dsk_Root[MAX_FNL+1];

  /* The subdirectory to write files to */
  char    Par_Dsk_SubDirec[MAX_FNL+1];

  /* The name of the scalar mesh file. */
  char    Scalar_Exo_File_Name[MAX_FNL+1];

  /* The base name of the parallel results file. */
  char    Par_Exo_Res_File_Name[MAX_FNL+1];

  /* The filename extension for the parallel files */
  char    Exo_Extension[MAX_FNL+1];
  
  char    Staged_Writes[5];

};

extern struct Parallel_IO PIO_Info;

extern char    Par_Nem_File_Name[];        /* The parallel nemesis file name */

int which_file(char *, char *, char **, char **);
int read_pexoII_info(const char *);
void gen_disk_map(struct Parallel_IO *pio_info, int proc_info[], int proc, int nproc);
void gen_par_filename(char *scalar_fname, char *par_fname, int proc_for, int nprocs);
void add_fname_ext(char *cOrigFile, const char *cExt);
#endif
