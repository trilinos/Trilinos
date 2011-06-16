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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pe_common.h"
#include "rf_io_const.h"
#include "rf_message.h"
#include "rf_allo.h"
#include "rf_util.h"

#include "rf_mp_const.h"
#include "el_geom_const.h"
#include "ps_pario_const.h"

/*********** R O U T I N E S   I N    T H I S  F I L E ***********************

  Name_of_Routine            type             Called by
  ---------------------------------------------------------------------------
  pdisk_stage_begin()        void        load_mesh        (pe_exoII_io.c)
  pdisk_stage_end()          void        load_mesh        (pe_exoII_io.c)
  gen_par_filename()         void        load_lb_info     (rf_load_lb_info.c)
                                         load_mesh        (pe_exoII_io.c)
          add_fname_ext()    void        gen_par_filename (pe_util.c)

******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void pdisk_stage_begin(struct Parallel_IO *pio_info, int my_proc_id,
                       int proc_info[], int proc_ids[], int pindx,
                       int *glob_pindx)
{
  int    i=0, disk_id, iproc, itype, st, first_pos=0;

  if(strcmp(pio_info->Staged_Writes, "yes") == 0) {

    /* Get the disk this processor will be writing to */
    disk_id = pio_info->RDsk_List[proc_ids[pindx]][0];

    /* Find my index in the global list */
    for(iproc=pindx; iproc < Proc_Info[0]; iproc++) {
      if(my_proc_id == pio_info->RDsk_List[iproc][1]) {
        *glob_pindx = iproc;
        break;
      }
    }

    /*
     * Find out if this processor has first access to the required
     * RAID by scanning the list.
     */
    for(iproc=0; iproc < Proc_Info[0]; iproc++) {
      if(pio_info->RDsk_List[iproc][0] == disk_id) {
        first_pos = iproc;
        break;
      }
    }

    if(first_pos == proc_ids[pindx]) {
      if(Debug_Flag >= 2) {
        printf("staging: Processor %d proceeding to write for %d...\n",
               Proc, proc_ids[pindx]);
      }
      return;
    }
    else {

      /*
       * Find the processor ID which is writing to the disk before
       * this processor.
       */
      for(iproc=(proc_ids[pindx]-1); iproc >= 0; iproc--) {
        if(pio_info->RDsk_List[iproc][0] == disk_id) {
          break;
        }
      }

      iproc = pio_info->RDsk_List[iproc][1];
      itype = EXOII_MESH_INT;

      if(Debug_Flag >= 4)
        printf("staging: Processor %d waiting for Processor %d...\n",
               my_proc_id, iproc);

      nread_big((char *) &i, sizeof(int), &iproc, &itype, &st);

      if(Debug_Flag >= 4)
        printf("staging: Processor %d proceeding...\n", my_proc_id);
    }
  }

  return;

} /********************* END pdisk_stage_begin() *****************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void pdisk_stage_end(struct Parallel_IO *pio_info, int my_proc_id,
                     int proc_info[], int proc_ids[], int pindx,
                     int glob_pindx)
{
  int    i=0, iproc, itype, st;

  if(strcmp(pio_info->Staged_Writes, "yes") == 0) {

    for(iproc=glob_pindx+1; iproc < Proc_Info[0]; iproc++) {
      if(pio_info->RDsk_List[iproc][0] ==
         pio_info->RDsk_List[glob_pindx][0])
        break;
    }

    if(iproc < Proc_Info[0]) {
      iproc = pio_info->RDsk_List[iproc][1];
      itype = EXOII_MESH_INT;

      if(Debug_Flag > 3)
        printf("staging: Processor %d telling Processor %d to "
               "proceed...\n", my_proc_id, iproc);

      nwrite_big((char *) &i, sizeof(int), iproc, itype, &st);
    }
  }

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void gen_disk_map(struct Parallel_IO *pio_info, int proc_info[],
                  int proc, int nproc)
/*
 * This function generates a map of which processor ID writes to which
 * RAID. Note that this is for each processor in the list, not necessarily
 * the same as the number of processors being used during the run. Each
 * processor has an identical list.
 */
{
  char   yo[]="gen_disk_map";

  int    i1, iproc, proc_id, ctrl_id;
/*------------------------ EXECUTION BEGINS ------------------------------*/

  /* Allocate memory for the list */
  pio_info->RDsk_List = (int **)array_alloc(__FILE__, __LINE__, 2,
                                            proc_info[0], 2, sizeof(int));
  if(!(pio_info->RDsk_List)) {
    fprintf(stderr, "%s: ERROR, insufficient memory\n", yo);
    exit(1);
  }

  /* Generate the list of disks to which data will be written */
  if(pio_info->Dsk_List_Cnt <= 0) {
    for(iproc=0; iproc < proc_info[0]; iproc++) {
      ctrl_id = (iproc % pio_info->Num_Dsk_Ctrlrs);
      pio_info->RDsk_List[iproc][0] = ctrl_id + pio_info->PDsk_Add_Fact;
    }
  }
  else {
    for(iproc=0; iproc < proc_info[0]; iproc++) {
      pio_info->RDsk_List[iproc][0] = pio_info->Dsk_List[iproc%
                                                   pio_info->Dsk_List_Cnt];
    }
  }

  /* Generate the list of processors on which info is stored */
  for(iproc=0; iproc < proc_info[0]; iproc++) {
    proc_id = iproc;
    while(proc_id >= nproc)
      proc_id -= nproc;

    pio_info->RDsk_List[iproc][1] = proc_id;
  }

  /* With a high enough Debug_Flag write the disk list to the screen */
  if(Debug_Flag >= 5) {
    print_sync_start(Proc, Num_Proc, FALSE);
    if(Proc == 0) {
      print_line("=", 79);
      printf("\t\tTABLE OF RAIDS TO USE\n\n");
      printf("\tProcessor\t|\tDisk list \"diskID(processor)\"\n");
      print_line("-", 79);
    }
    for(iproc=0; iproc < proc_info[2]; iproc++) {
      printf("\t  %d\t\t|\t", Proc_Ids[iproc]);
      for(i1=0; i1 < proc_info[0]; i1++) {
        printf(" %d(%d)", pio_info->RDsk_List[i1][0],
               pio_info->RDsk_List[i1][1]);
      }

      printf("\n");
    }

    print_sync_end(Proc, Num_Proc, FALSE);
  }

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void gen_par_filename(char *scalar_fname, char *par_fname,
                      int proc_for, int nprocs)
/*----------------------------------------------------------------------------
 *
 *      Author(s):     Gary Hennigan (1421)
 *----------------------------------------------------------------------------
 *      Function which generates the name of a parallel file for a
 *      particular processor. The function does this by appending
 *      "N.p" to the end of the input parameter "scalar_fname", where:
 *
 *              N - The number of processors utilized
 *              p - The processor ID.
 *
 *      In addition, the location of the parallel disk system is prepended
 *      to each file name.
 *---------------------------------------------------------------------------
 *      Example:
 *
 *        scalar_fname = "Parallel-exoII-"   (Input)
 *        par_fname    = "/raid/io_01/tmp/rf_crew/Parallel-exoII-8.0" (Output)
 *
 *      where, for this example:
 *
 *              N = 8 processors
 *              p = 0 particular processor ID
 *---------------------------------------------------------------------------
 *      Revision History:
 *
 *              05 November 1993:    Date of Creation
 *---------------------------------------------------------------------------
 */
{

  /*      Local variables      */

  int i1, iTemp1, ctrlID;
  int iMaxDigit=0, iMyDigit=0;
  char cTemp[MAX_FNL];

/************************* EXECUTION BEGINS *******************************/

  /*
   * Find out the number of digits needed to specify the processor ID.
   * This allows numbers like 01-99, i.e., prepending zeros to the
   * name to preserve proper alphabetic sorting of the files.
   */

  iTemp1 = nprocs;
  do
  {
    iTemp1 /= 10;
    iMaxDigit++;
  }
  while(iTemp1 >= 1);

  iTemp1 = proc_for;
  do
  {
    iTemp1 /= 10;
    iMyDigit++;
  }
  while(iTemp1 >= 1);

  /*
   * Append the number of processors in this run to the scalar file name
   * along with a '.' (period).
   */
  par_fname[0] = 0x00;
  strcpy(par_fname, scalar_fname);
  strcat(par_fname, ".");
  sprintf(cTemp, "%d", nprocs);
  strcat(par_fname, cTemp);
  strcat(par_fname, ".");

  /*
   * Append the proper number of zeros to the filename.
   */
  for(i1=0; i1 < iMaxDigit-iMyDigit; i1++)
    strcat(par_fname, "0");

  /*
   * Generate the name of the directory on which the parallel disk
   * array resides. This also directs which processor writes to what
   * disk.
   */
  sprintf(cTemp, "%d", proc_for);
  strcat(par_fname, cTemp);
  strcpy(cTemp, par_fname);


  /*
   * Finally, generate the complete file specification for the parallel
   * file used by this processor.
   */
  if(PIO_Info.NoSubdirectory == 1) {
    sprintf(par_fname, "%s%s%s",
	    PIO_Info.Par_Dsk_Root, PIO_Info.Par_Dsk_SubDirec, cTemp);
  } else {
    if(PIO_Info.Zeros) {
      ctrlID = PIO_Info.RDsk_List[proc_for][0];
      if(ctrlID <= 9) {
	sprintf(par_fname, "%s%d%d/%s%s", PIO_Info.Par_Dsk_Root,0,
		ctrlID, PIO_Info.Par_Dsk_SubDirec, cTemp);
      }
      else {
	sprintf(par_fname, "%s%d/%s%s", PIO_Info.Par_Dsk_Root,
		ctrlID, PIO_Info.Par_Dsk_SubDirec, cTemp);
      }
  }
    else {
      ctrlID = PIO_Info.RDsk_List[proc_for][0];
      sprintf(par_fname, "%s%d/%s%s", PIO_Info.Par_Dsk_Root, ctrlID,
	      PIO_Info.Par_Dsk_SubDirec, cTemp);
    }
  }
/* not supporting ncubed right now ------>
  else if(strcmp(PIO_Info.Targ_Machine,"ncube") == 0) {
    diskID = proc_for % PIO_Info.Num_Dsks_PCtrlr;
    ctrlID = (proc_for / PIO_Info.Num_Dsks_PCtrlr) % PIO_Info.Num_Dsk_Ctrlrs;
    sprintf(par_fname, "%s%d%d/%s%s", PIO_Info.Par_Dsk_Root,
            ctrlID, diskID+PIO_Info.PDsk_Add_Fact, PIO_Info.Par_Dsk_SubDirec,
            cTemp);

  }
<-------------------*/

  if(Debug_Flag >= 4)
    printf("[%d]: Parallel file name: %s\n", Proc, par_fname);

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void add_fname_ext(char *cOrigFile, char *cExt)
/*
 *     This function removes the existing extension on a file name,
 *     if it exists, and adds the extension input to the function as
 *     the variable cExt. The function overwrites the original file
 *     name.
 *
 *     Note that it is assumed enough memory is allocate for the original
 *     string to handle it's new extension.
 */
{

  char *cPtr;
  int i1, iExtLen;

  cPtr = strrchr(cOrigFile, '.');

  iExtLen = strlen(cExt);

  if(cPtr != NULL)
  {
    for(i1=0; i1 < iExtLen; i1++)
    {
      cPtr[i1] = cExt[i1];
    }
    cPtr[i1] = '\0';
  }
  else
  {
    cPtr = cOrigFile;
    cPtr += strlen(cOrigFile);
    for(i1=0; i1 < iExtLen; i1++)
    {
      cPtr[i1] = cExt[i1];
    }
    cPtr[i1] = '\0';
  }

  return;
}
