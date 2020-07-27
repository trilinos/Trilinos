/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "fmt/ostream.h"
#include "ps_pario_const.h" // for Parallel_IO, PIO_Info
#include "rf_allo.h"        // for array_alloc
#include "rf_io_const.h"    // for Debug_Flag
#include <cstdlib>          // for exit
#include <cstring>          // for strlen, etc
#include <sstream>
#include <string>

/*********** R O U T I N E S   I N    T H I S  F I L E ***********************

  Name_of_Routine            type             Called by
  ---------------------------------------------------------------------------
  pdisk_stage_begin()        void        load_mesh        (pe_exoII_io.c)
  pdisk_stage_end()          void        load_mesh        (pe_exoII_io.c)
  gen_par_filename()         void        load_lb_info     (rf_load_lb_info.c)
                                         load_mesh        (pe_exoII_io.c)

******************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void gen_disk_map(struct Parallel_IO *pio_info, int proc_info[], int /*proc*/, int nproc)
/*
 * This function generates a map of which processor ID writes to which
 * RAID. Note that this is for each processor in the list, not necessarily
 * the same as the number of processors being used during the run. Each
 * processor has an identical list.
 */
{
  int iproc;
  int proc_id;
  int ctrl_id;
  /*------------------------ EXECUTION BEGINS ------------------------------*/

  /* Allocate memory for the list */
  pio_info->RDsk_List =
      reinterpret_cast<int **>(array_alloc(__FILE__, __LINE__, 2, proc_info[0], 2, sizeof(int)));
  if ((pio_info->RDsk_List) == nullptr) {
    fmt::print(stderr, "{}: ERROR, insufficient memory\n", __func__);
    exit(1);
  }

  /* Generate the list of disks to which data will be written */
  if (pio_info->Dsk_List_Cnt <= 0) {
    for (iproc = 0; iproc < proc_info[0]; iproc++) {
      ctrl_id                       = (iproc % pio_info->Num_Dsk_Ctrlrs);
      pio_info->RDsk_List[iproc][0] = ctrl_id + pio_info->PDsk_Add_Fact;
    }
  }
  else {
    for (iproc = 0; iproc < proc_info[0]; iproc++) {
      pio_info->RDsk_List[iproc][0] = pio_info->Dsk_List[iproc % pio_info->Dsk_List_Cnt];
    }
  }

  /* Generate the list of processors on which info is stored */
  for (iproc = 0; iproc < proc_info[0]; iproc++) {
    proc_id = iproc;
    while (proc_id >= nproc) {
      proc_id -= nproc;
    }

    pio_info->RDsk_List[iproc][1] = proc_id;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
std::string gen_par_filename(const std::string &scalar_fname, int proc_for, int nprocs)
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

  int         i1;
  int         iTemp1;
  int         ctrlID;
  int         iMaxDigit = 0;
  int         iMyDigit  = 0;
  std::string par_filename;

  /************************* EXECUTION BEGINS *******************************/

  /*
   * Find out the number of digits needed to specify the processor ID.
   * This allows numbers like 01-99, i.e., prepending zeros to the
   * name to preserve proper alphabetic sorting of the files.
   */

  iTemp1 = nprocs;
  do {
    iTemp1 /= 10;
    iMaxDigit++;
  } while (iTemp1 >= 1);

  iTemp1 = proc_for;
  do {
    iTemp1 /= 10;
    iMyDigit++;
  } while (iTemp1 >= 1);

  /*
   * Append the number of processors in this run to the scalar file name
   * along with a '.' (period).
   */
  par_filename = scalar_fname + std::string(".") + std::to_string(nprocs) + std::string(".");

  /*
   * Append the proper number of zeros to the filename.
   */
  for (i1 = 0; i1 < iMaxDigit - iMyDigit; i1++) {
    par_filename += std::string("0");
  }

  /*
   * Generate the name of the directory on which the parallel disk
   * array resides. This also directs which processor writes to what
   * disk.
   */
  par_filename += std::to_string(proc_for);

  /*
   * Finally, generate the complete file specification for the parallel
   * file used by this processor.
   */
  if (PIO_Info.NoSubdirectory == 1) {
    par_filename = PIO_Info.Par_Dsk_Root + PIO_Info.Par_Dsk_SubDirec + par_filename;
  }
  else {
    if (PIO_Info.Zeros != 0) {
      ctrlID = PIO_Info.RDsk_List[proc_for][0];
      if (ctrlID <= 9) {
        par_filename = PIO_Info.Par_Dsk_Root + "0" + std::to_string(ctrlID) + "/" +
                       PIO_Info.Par_Dsk_SubDirec + par_filename;
      }
      else {
        par_filename = PIO_Info.Par_Dsk_Root + std::to_string(ctrlID) + "/" +
                       PIO_Info.Par_Dsk_SubDirec + par_filename;
      }
    }
    else {
      ctrlID       = PIO_Info.RDsk_List[proc_for][0];
      par_filename = PIO_Info.Par_Dsk_Root + std::to_string(ctrlID) + "/" +
                     PIO_Info.Par_Dsk_SubDirec + par_filename;
    }
  }
  if (Debug_Flag >= 4) {
    fmt::print("Parallel file name: {}\n", par_filename.c_str());
  }

  return par_filename;
}
