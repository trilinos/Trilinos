/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef PS_PARIO_CONST_H
#define PS_PARIO_CONST_H

#include "rf_io_const.h"
#include <string>

/* Global variables. */
extern double PIO_Time_Array[]; /* Vector for timings */

/*
 * The following variables are used when a single processor is to write
 * info for a processor other than itself.
 */

/* Structure used to store the information necessary for parallel I/O. */

struct Parallel_IO
{
  int Dsk_List_Cnt;

  int * Dsk_List;
  int **RDsk_List;

  int Num_Dsk_Ctrlrs;  /* The number of disk controllers.     */
  int Num_Dsks_PCtrlr; /* The number of disks per controller. */
  int PDsk_Add_Fact;   /* The offset from zero used by the    */
                       /* the target machine.                 */

  int Zeros; /* 1 - if the target machine uses leading zeros when */
             /*     designating the disk number (eg - the paragon */
             /*     uses /pfs/io_01)                              */
             /* 0 - if it does not (eg - the tflop uses           */
             /*     /pfs/tmp_1)                                   */

  /* 1 - don't create a subdirectory for the spread files; write
   * them in same directory as the mesh file
   */
  int NoSubdirectory;

  /* The root location of the parallel disks */
  std::string Par_Dsk_Root;

  /* The subdirectory to write files to */
  std::string Par_Dsk_SubDirec;

  /* The filename extension for the parallel files */
  std::string Exo_Extension;

  bool Staged_Writes;
};

extern struct Parallel_IO PIO_Info;

extern std::string Par_Nem_File_Name; /* The parallel nemesis file name */

void        gen_disk_map(struct Parallel_IO *pio_info, int proc_info[], int proc, int nproc);
std::string gen_par_filename(const std::string &scalar_fname, int proc_for, int nprocs);
#endif
