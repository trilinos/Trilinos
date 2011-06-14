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
/*--------------------------------------------------------------------------*/
/* Purpose: Parallel I/O utility for SALSA.                                 */
/*                                                                          */
/*    This program takes an Exodus II FEM database file and a parallel      */
/* load-balance file and creates a set of parallel Exodus II FEM files for  */
/* use by SALSA.                                                            */
/*--------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rf_salsa.h"
#include "rf_comm.h"
#include "rf_allo.h"

#include "el_geom_const.h"
#include "el_geom.h"

#include "rf_io_const.h"
#include "rf_io.h"

#include "pe_init_struct.h"
#include "ps_pario_const.h"
#include "ps_pario.h"

#include "netcdf.h"
#include "exodusII.h"

#include "rf_mp.h"

#include "rf_fem_const.h"
#include "rf_fem.h"

#if defined(USE_MPI)
#include <mpi.h>
#endif

extern void brdcst_command_info(void);
extern void add_to_log(const char *name);

int main(int argc, char *argv[])
{
/* Local declarations. */
  double       g_start_t, g_end_t, start_t, end_t;
  char         *salsa_cmd_file;
  int          i1, io_ws;
  static char yo[] = "nem_spread";
/***************************** BEGIN EXECUTION ******************************/

#if defined(USE_MPI)
  MPI_Init(&argc, &argv);
#endif

  g_start_t = second();

  /* Determine Processor Number and size of Parallel Machine */
  get_parallel_info(&Proc, &Num_Proc, &Dim);

  /* Scan the command line arguments for a version flag */
  for(i1=0; i1 < argc; i1++) {
    if(strcmp(argv[i1],"-V") == 0) {
      printf("%s version %s\n", UTIL_NAME, VER_STR);
      exit(0);
    }
  }


  /* Interpret the command line */
  switch(argc)
  {
  case 1:
    salsa_cmd_file = "nem_spread.inp";
    break;

  case 2:
    salsa_cmd_file = argv[1];
    break;

  case 4:
    salsa_cmd_file = argv[1];
    if(strstr(argv[2], "-p")) {
      i1 = sscanf(argv[3], "%d", &Proc_For);
      if(i1 != 1) {
        fprintf(stderr, "ERROR: Incorrect command line!\n");
        fprintf(stderr, " usage:\n");
        fprintf(stderr, "\tnem_spread [salsa command file] [<-p Proc> ");
        fprintf(stderr, "<-r raid #>]\n");
        exit(1);
      }
    }
    else {
      fprintf(stderr, "ERROR: Incorrect command line!\n");
      fprintf(stderr, " usage:\n");
      fprintf(stderr, "\tnem_spread [salsa command file] ");
      fprintf(stderr, "[<-p Proc> <-r raid #>]\n");
      exit(1);
    }
    break;

  default:
    fprintf(stderr, "%s MAIN: ERROR in command line,", yo);
    if(Proc == 0)
    {
      fprintf(stderr, " usage:\n");
      fprintf(stderr, "\tnem_spread [salsa command file] [<-p Proc> ");
      fprintf(stderr, "<-r raid #>]");
    }
    exit(1);
    break;
  }

  if(Proc == 0) {
    printf("%s version %s\n", UTIL_NAME, VER_STR);
    if (Num_Proc > 1) {
      printf("Total number of Processors = %d\n", Num_Proc);
    }
  }

  /* initialize some variables */
  ExoFile[0]                    = '\0';
  Exo_LB_File[0]                = '\0';
  Exo_Res_File[0]               = '\0';
  Debug_Flag                    = -1;

  Restart_Info.Flag             = -1;
  Restart_Info.Num_Times        = -1;
  Restart_Info.Block_Size       = -1;

  Num_Nod_Var                   = -1;
  Num_Elem_Var                  = -1;
  Num_Glob_Var                  = -1;
  Num_Nset_Var                  = -1;
  Num_Sset_Var                  = -1;

  PIO_Info.Dsk_List_Cnt         = -1;
  PIO_Info.Num_Dsk_Ctrlrs       = -1;
  PIO_Info.PDsk_Add_Fact        = -1;
  PIO_Info.Zeros                = -1;
  PIO_Info.NoSubdirectory       =  0;
  PIO_Info.Par_Dsk_Root[0]      = '\0';
  PIO_Info.Par_Dsk_SubDirec[0]  = '\0';
  PIO_Info.Staged_Writes[0]     = '\0';


  /*
   * Read in the ASCII input file from the front end.
   * NOTE: In this case we read only the information needed to create the
   *       parallel exodus files.
   */
  if(Proc == 0)
  {
    printf("Reading the command file, %s\n", salsa_cmd_file);
    if(read_pexoII_info(salsa_cmd_file) < 0)
    {
      fprintf(stderr,"%s ERROR: Could not read in the the I/O command file"
              " \"%s\"!\n", yo, salsa_cmd_file);
      exit(1);
    }

    if (!check_inp()) {
      fprintf(stderr, "%s ERROR: Error in user specified parameters.\n", yo);
      exit(1);
    }
  }

  /*
   * If the number of processors is greater than 1 then broadcast the
   * needed information read by processor 0 from the ASCII command
   * files to the other processors.
   */
  if(Num_Proc > 1)
    brdcst_command_info();
  else {
    strcpy(PIO_Info.Scalar_LB_File_Name, Exo_LB_File);
    strcpy(PIO_Info.Scalar_Exo_File_Name, ExoFile);
  }

  /*
   * Exchange Global Information garnered by Proc = 0
   * with all other processors.
   */
  if (Num_Proc > 1) exch_init_info ();

  /* If debug is on the turn on netCDF/Exodus information as well */
  if(Debug_Flag > 0)
     ex_opts(EX_VERBOSE | EX_DEBUG);

  /*
   * Read initial information from the mesh file.
   *  - This provides error checking against the load balance file
   *  - Broadcast the information
   */
  read_mesh_param(&io_ws);


  /*
   * Process the load balance information
   *  - Read info from the lb file with Proc 0
   *  - Distribute the information to all processors
   */
  start_t = second ();
  load_lb_info ();
  end_t   = second () - start_t;
  if (Proc == 0) printf ("\nLoad load balance information time: %f (sec.)\n\n", end_t);


  /*
   * Get any restart parameter information
   *  - Read the parameters from the input ExodusII file
   */
  if(Restart_Info.Flag > 0) {
    if (Proc == 0)
      printf ("Load exoII restart param info to each proc.\n\n");
    start_t = second ();
    read_restart_params(io_ws);
    end_t   = second () - start_t;
    if (Proc == 0)
      printf ("Load restart parameters time: %f (sec.)\n\n", end_t);
  }


  /*
   * Read the ExodusII mesh file and distribute information
   * contained in it to all processors
   *  - Each processor only gets the information that it needs
   *    to solve its piece of the mesh
   */
  if (Proc == 0)
    printf ("Load exoII mesh info to each proc.\n\n");
  start_t = second ();
  load_mesh (io_ws);
  end_t   = second () - start_t;

  if (Proc == 0)
    printf ("Load mesh time: %f (sec.)\n\n", end_t);

  /*
   * Get any restart variable data
   *  - Read the restart data from the input ExodusII file, distribute
   *    it to the processors, and write it to the parallel files
   */
  if(Restart_Info.Flag > 0) {
    if (Proc == 0)
      printf ("Load exoII restart data info to each proc.\n\n");
    start_t = second ();
    read_restart_data(io_ws);
    end_t   = second () - start_t;
    if (Proc == 0)
      printf ("Load restart data time: %f (sec.)\n\n", end_t);
  }

  if (Proc == 0)
    printf ("Write of parallel exodus complete\n");

  g_end_t = second() - g_start_t;

  g_start_t = gsum_double(g_end_t, Proc, Num_Proc);

  g_end_t = g_start_t / (float)Num_Proc;

  if(Proc == 0)
    printf("The average run time was: %.2fs\n", g_end_t);

  safe_free((void**)&(Proc_Ids));
  safe_free((void**)&(PIO_Info.RDsk_List));
  
#if defined(USE_MPI)
  MPI_Finalize();
#endif
  add_to_log(argv[0]);
  return 0;
}
