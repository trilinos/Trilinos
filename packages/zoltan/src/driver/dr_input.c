/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include <mpi.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"
#include "ch_init_dist_const.h"

/*--------------------------------------------------------------------------*/
/* Purpose: Determine file types for command files and read in the parallel */
/*          ExodusII command file.                                          */
/*          Taken from nemesis utilites nem_spread and nem_join.            */
/*--------------------------------------------------------------------------*/

#define TLIST_CNT 5

/*****************************************************************************/
/*****************************************************************************/
int read_cmd_file(char *filename, PROB_INFO_PTR prob,
                  PARIO_INFO_PTR pio_info)

/*
 *          This function reads the ASCII parallel-exodus command file.
 *
 *   Input
 *   -----
 *   filename - The name of the command file.
 *   pio_info - parallel I/O information.
 */
{
/* local declarations */
  FILE  *file_cmd;
  char  *yo = "read_cmd_file";
  char   cmesg[256]; /* for error messages */
  char   inp_line[MAX_INPUT_STR_LN + 1];
  char   inp_copy[MAX_INPUT_STR_LN + 1];
  char  *cptr, *cptr2;
  int    i, icnt;
  int    num_params_alloc = 0, param_index;

/***************************** BEGIN EXECUTION ******************************/

  /* Open the file */
  if((file_cmd=fopen(filename, "r")) == NULL)
    return 0;

  /* Begin parsing the input file */
  while(fgets(inp_line, MAX_INPUT_STR_LN, file_cmd)) {
    /* skip any line that is a comment */
    if((inp_line[0] != '#') && (inp_line[0] != '\n')) {

      strcpy(inp_copy, inp_line);
      clean_string(inp_line, " \t");
      cptr = strtok(inp_line, "\t=");

      /****** The file type ******/
      if (token_compare(cptr, "file type")) {
        if(pio_info->file_type < 0) {
          cptr = strtok(NULL, "\t=,");
          strip_string(cptr, " \t\n");
          if (cptr == NULL || strlen(cptr) == 0) {
            Gen_Error(0, "fatal: must specify file type");
            return 0;
          }

          string_to_lower(cptr, '\0');
          if (strstr(cptr, "nemesisi")) {
            pio_info->file_type = NEMESIS_FILE;
            pio_info->init_dist_type = INITIAL_FILE;
          }
          else if (strstr(cptr, "chaco")) {
            pio_info->file_type = CHACO_FILE;
            /* Set default values for options here and change them later */
            pio_info->init_dist_type = INITIAL_LINEAR;
            pio_info->init_dist_procs = -1;
            while ((cptr = strtok(NULL, ",")) != NULL){ 
              strip_string(cptr, " \t\n");
              string_to_lower(cptr, '=');
              if (strstr(cptr, "initial distribution")) {
                cptr2 = strchr(cptr, '=');
                if (cptr2 == NULL) {
                  Gen_Error(0, "fatal: initial distribution type is not "
                               "specified");
                  return 0;
                }
                cptr2++;
                string_to_lower(cptr2, '\n');
                if (strstr(cptr2, "linear") || strstr(cptr2, "block")) {
                  pio_info->init_dist_type = INITIAL_LINEAR;
                }
                else if (strstr(cptr2, "cyclic")) {
                  pio_info->init_dist_type = INITIAL_CYCLIC;
                }
                else if (strstr(cptr2, "file")) {
                  pio_info->init_dist_type = INITIAL_FILE;
                }
                else {
                  Gen_Error(0, "Invalid Chaco initial distribution type.");
                  return 0;
                }
              }
              else if (strstr(cptr, "initial procs")) {
                cptr2 = strchr(cptr, '=');
                if (cptr2 == NULL) {
                  Gen_Error(0, "fatal: initial no. of procs is not "
                               "specified");
                  return 0;
                }
                cptr2++;
                if(sscanf(cptr2, "%d", &(pio_info->init_dist_procs)) != 1) {
                  Gen_Error(0, "fatal: initial no. of procs must be an integer.");
                  return 0;
                }
              }
            }
          }
          else if (strstr(cptr, "hypergraph")) {
            pio_info->file_type = HYPERGRAPH_FILE;
            /* Set default values for options here and change them later */
            pio_info->init_dist_type = INITIAL_LINEAR;
            pio_info->init_dist_procs = -1;
          }
          else {
            sprintf(cmesg, "fatal(%s): unknown file type, %s", yo, cptr);
            Gen_Error(0, cmesg);
            return 0;
          }
        }
      }

      /****** The file name ******/
      else if (token_compare(cptr, "file name")) {
        if(strlen(pio_info->pexo_fname) == 0)
        {
          cptr = strtok(NULL, "\t=");
          strip_string(cptr, " \t\n");
          strcpy(pio_info->pexo_fname, cptr);
        }
      }

      /****** The Debug reporting level ******/
      else if (token_compare(cptr, "zdrive debug level")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Debug_Driver) != 1) {
          Gen_Error(0, "fatal: zdrive debug level must be an integer.");
          return 0;
        }
      }

      /****** List-based (MULTI) callback function testing ******/
      else if (token_compare(cptr, "test multi callbacks")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Test.Multi_Callbacks) != 1) {
          Gen_Error(0, "fatal: test multi callbacks must be an integer.");
          return 0;
        }
      }

      /****** Null import lists to Help_Migrate testing******/
      else if (token_compare(cptr, "test null import lists")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Test.Null_Import_Lists) != 1) {
          Gen_Error(0, "fatal: test null import lists must be an integer.");
          return 0;
        }
      }

      /****** Box- and Point-drop testing flag ******/
      else if (token_compare(cptr, "test drops")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Test.Drops) != 1) {
          Gen_Error(0, "fatal: test drops must be an integer.");
          return 0;
        }
      }

      /****** DDirectory testing flag ******/
      else if (token_compare(cptr, "test ddirectory")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Test.DDirectory) != 1) {
          Gen_Error(0, "fatal: test ddirectory must be an integer.");
          return 0;
        }
      }

      /****** Unusual Partition generation testing flag ******/
      else if (token_compare(cptr, "test local partitions")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Test.Local_Partitions) != 1) {
          Gen_Error(0, "fatal: test partitions must be an integer.");
          return 0;
        }
      }

      /****** Generate GNUplot output? ******/
      else if (token_compare(cptr, "gnuplot output")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Output.Gnuplot) != 1) {
          Gen_Error(0, "fatal: gnuplot output indicator must be an integer.");
          return 0;
        }
      }

      /****** Plot processor numbers or partitions? ******/
      else if (token_compare(cptr, "plot partitions")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Output.Plot_Partitions) != 1) {
          Gen_Error(0, "fatal: plot partitions indicator must be an integer.");
          return 0;
        }
      }

      /****** Generate Nemesis output? ******/
      else if (token_compare(cptr, "nemesis output")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Output.Nemesis) != 1) {
          Gen_Error(0, "fatal: nemesis output indicator must be an integer.");
          return 0;
        }
      }

      /****** Generate ASCII mesh file? ******/
      else if (token_compare(cptr, "print mesh info file")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Output.Mesh_Info_File) != 1) {
          Gen_Error(0, "fatal: Print Mesh Info File indicator must be an integer.");
          return 0;
        }
      }

      /****** The Number of iterations of the balancer to perform ******/
      else if (token_compare(cptr, "number of iterations")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Number_Iterations) != 1) {
          Gen_Error(0, "fatal: number of iterations must be an integer.");
          return 0;
        }
      }

      /****** zdrive action: Do load-balancing or ordering? *****/
      else if (token_compare(cptr, "zdrive action")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Driver_Action) != 1) {
          Gen_Error(0, "fatal: zdrive action must be an integer.");
          return 0;
        }
      }

      /****** The Chaco Input Assignment Inverse flag ******/
      else if (token_compare(cptr, "chaco input assignment inverse")) {
        cptr = strtok(NULL, "\t=");
        strip_string(cptr, " \t\n");
        if(sscanf(cptr, "%d", &Chaco_In_Assign_Inv) != 1) {
          Gen_Error(0, "fatal: Chaco Input Assignment Inverse must be an integer.");
          return 0;
        }
      }

      /****** Parallel Disk Information ******/
      else if (token_compare(cptr, "parallel disk info")) {

        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr," \t\n=");
        cptr = strtok(cptr, ",");
        strip_string(cptr, " \t\n");
        string_to_lower(cptr, '=');

        /* the first sub-option must be "number" */
        if (!strstr(cptr, "number")) {
          Gen_Error(0, "fatal: First sup-option for disk info must be "
                       "\"number\"");
          return 0;
        }
        else {
          cptr2 = strchr(cptr, '=');
          if (cptr2 == NULL) {
            Gen_Error(0, "fatal: integer value must be specified for"
                         " reserve space.");
            return 0;
          }
          cptr2++;
          icnt = sscanf(cptr2, "%d", &(pio_info->num_dsk_ctrlrs));
          if ((icnt <= 0) || (pio_info->num_dsk_ctrlrs < 0)) {
            Gen_Error(0, "fatal: Invalid value for # of raid controllers.");
            return 0;
          }
        }

        cptr = strtok(NULL, ",");

        /*
         * if number = 0, then the input file(s) is in the
         * root directory given by the parallel disk info, or it
         * is in the same directory as the executable if nothing
         * is given for the root infomation. So, no other options
         * can be given when number = 0
         */
        if (pio_info->num_dsk_ctrlrs == 0 && cptr != NULL) {
          Gen_Error(0, "fatal: Other options not allowed if number = 0.");
          return 0;
        }

        while (cptr != NULL) {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '=');
          if (strstr(cptr, "list")) {
            /*
             * So, "number" references the length of the list, and
             * I need to do some shuffling to make the new form
             * work with the old code.
             */
            pio_info->dsk_list_cnt = pio_info->num_dsk_ctrlrs;
            pio_info->num_dsk_ctrlrs = -1;

            /* "{" defines the beginning of the list. } */
            cptr = strchr(cptr, '{');
            if (cptr == NULL) {
              Gen_Error(0, "fatal: disk list must be specified.");
              return 0;
            }
            cptr++;

            /* allocate memory for to hold the values */
            pio_info->dsk_list = (int *) malloc(pio_info->dsk_list_cnt *
                                               sizeof(int));
            for (i = 0; i < (pio_info->dsk_list_cnt - 1); i++) {
              sscanf(cptr, "%d", &(pio_info->dsk_list[i]));
              cptr = strtok(NULL, ", \t;");
            }
            /* last one is a special case */
            sscanf(cptr, "%d}", &(pio_info->dsk_list[i]));
          }
          else if (strstr(cptr, "offset")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == NULL) {
              Gen_Error(0, "fatal: value must be specified with the "
                           "\"offset\" option.");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &(pio_info->pdsk_add_fact));
            if ((icnt <= 0) || (pio_info->pdsk_add_fact < 0)) {
              Gen_Error(0, "fatal: Invalid value for offset.");
              return 0;
            }
          }
          else if (strcmp(cptr, "zeros") == 0) {
            pio_info->zeros = 1;
          }

          cptr = strtok(NULL, ",");
        }
      } /* End "else if (token_compare(cptr, "parallel disk info"))" */

      /****** Parallel File Location ******/
      else if (token_compare(cptr, "parallel file location")) {
        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr," \t\n=");
        cptr = strtok(cptr, ",");

        while (cptr != NULL) {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '=');
          if (strstr(cptr, "root")) {
            cptr2 = strchr(cptr, '=');
            if(cptr2 == NULL)
            {
              Gen_Error(0, "fatal: must specify a path with \"root\"");
              return 0;
            }
            cptr2++;
            if(strlen(cptr2) == 0)
            {
              Gen_Error(0, "fatal: invalid path name specified with \"root\"");
              return 0;
            }
            strcpy(pio_info->pdsk_root, cptr2);
          }
          if (strstr(cptr, "subdir")) {
            cptr2 = strchr(cptr, '=');
            if(cptr2 == NULL)
            {
              Gen_Error(0, "fatal: must specify a path with \"subdir\"");
              return 0;
            }
            cptr2++;
            if(strlen(cptr2) == 0)
            {
              Gen_Error(0, "fatal: invalid path name specified with "
                           "\"subdir\"");
              return 0;
            }
            strcpy(pio_info->pdsk_subdir, cptr2);
            if (pio_info->pdsk_subdir[strlen(pio_info->pdsk_subdir)-1] != '/')
              strcat(pio_info->pdsk_subdir, "/");
          }

          cptr = strtok(NULL, ",");
        }
      }

      /****** Decomposition Info ******/
      else if (token_compare(cptr, "decomposition method")) {
        /* The method to use for decomposing the graph */
        if(strlen(prob->method) == 0)
        {
          cptr = strtok(NULL, "\t=");
          strip_string(cptr, " \t\n");
          strcpy(prob->method, cptr);
        }
      }

      /****** Zoltan Parameters ******/
      else if (token_compare(cptr, "zoltan parameters")) {
        /* parameters to be passed to Zoltan */
        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr," \t\n=");
        cptr = strtok(cptr, ",");

        /* parameters should be designated by "<param string>=<value string>" */
        while (cptr != NULL) {
          param_index = prob->num_params;
          prob->num_params++;

          /* ensure there is enough memory to store the new parameters */
          if (prob->num_params > num_params_alloc) {
            if (num_params_alloc == 0) {
              num_params_alloc += 5;
              prob->params = (Parameter_Pair *) malloc(
                                      num_params_alloc*sizeof(Parameter_Pair));
            }
            else {
              num_params_alloc += 5;
              prob->params = (Parameter_Pair *) realloc(prob->params,
                                      num_params_alloc*sizeof(Parameter_Pair));
            }
            if (prob->params == NULL) {
              Gen_Error(0, "fatal: not enough memory for Zoltan parameters");
              return 0;
            }
          }

          /* copy parameter name */
          i = strcspn(cptr, "=");
          strncpy(prob->params[param_index][0], cptr, i);
          prob->params[param_index][0][i]='\0';
          strip_string(prob->params[param_index][0], " \t\n");

          /* now get the value */
          cptr2 = strchr(cptr, '=');
          if(cptr2 == NULL)
          {
            sprintf(cmesg, "fatal(%s): must specify value for parameter %s", 
                            yo, prob->params[param_index][0]);
            Gen_Error(0, cmesg);
            return 0;
          }
          cptr2++;
          if (strlen(cptr2) == 0) {
            sprintf(cmesg, "fatal(%s): must specify value for parameter %s", 
                            yo, prob->params[param_index][0]);
            Gen_Error(0, cmesg);
            return 0;
          }
          strcpy(prob->params[param_index][1], cptr2);
          strip_string(prob->params[param_index][1], " \t\n");
          cptr = strtok(NULL, ",");
        }
      }
    } /* End "if(inp_line[0] != '#')" */
  } /* End "while(fgets(inp_line, MAX_INPUT_STR_LN, file_cmd))" */


  /* Close the command file */
  fclose(file_cmd);

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int check_inp(PROB_INFO_PTR prob, PARIO_INFO_PTR pio_info)
{
/***************************** BEGIN EXECUTION ******************************/

  /* check for the parallel Nemesis file for proc 0 */
  if (strlen(pio_info->pexo_fname) <= 0) {
    Gen_Error(0, "fatal: must specify file base name");
    return 0;
  }

  /* default file type is nemesis */
  if (pio_info->file_type < 0) pio_info->file_type = NEMESIS_FILE;

#ifndef ZOLTAN_NEMESIS
  /* 
   * if not compiling with the ZOLTAN_NEMESIS flag (i.e., not linking with 
   * Nemesis library), can't use NEMESIS_FILE file type.
   */

  if (pio_info->file_type == NEMESIS_FILE) {
    Gen_Error(0, "fatal: must compile for and link with Nemesis "
                 "libraries for Nemesis file types");
    return 0;
  }
#endif /* !ZOLTAN_NEMESIS */

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                 Check the parallel IO specifications                      */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /* check that there is a list of disks, or a number of raids */
  if ((pio_info->dsk_list_cnt <= 0) && (pio_info->num_dsk_ctrlrs < 0))
    pio_info->num_dsk_ctrlrs = 0; /* default to single directory */

  /* default is not to have preceeding 0's in the disk names */
  if (pio_info->zeros < 0) pio_info->zeros = 0;

  /* most systems that we deal with start their files systems with 1 not 0 */
  if (pio_info->pdsk_add_fact < 0) pio_info->pdsk_add_fact = 1;

  /*
   * if there are parallel disks, then the root and subdir locations must
   * be specified
   */
  if (pio_info->num_dsk_ctrlrs > 0 || pio_info->dsk_list_cnt > 0) {
    if (strlen(pio_info->pdsk_root) == 0) {
      Gen_Error(0, "fatal: must specify parallel disk root name");
      return 0;
    }
    if (strlen(pio_info->pdsk_subdir) == 0) {
      Gen_Error(0, "fatal: must specify parallel disk subdirectory");
      return 0;
    }
  }
  else
    if (strlen(pio_info->pdsk_root) == 0)
      strcpy(pio_info->pdsk_root, "."); /* default is execution directory */

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                 Check the Zoltan specifications                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*
   * Make sure a load-balancing method was provided.
   */
  if (strlen(prob->method) == 0) {
    Gen_Error(0, "fatal: load balance method must be specified");
    return 0;
  }

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void brdcst_cmd_info(int Proc, PROB_INFO_PTR prob, PARIO_INFO_PTR pio_info,
                     MESH_INFO_PTR mesh)
{
/* local declarations */
  int ctrl_id;
  int size;
  int int_params[12];  /* Make sure this array is large enough */
/***************************** BEGIN EXECUTION ******************************/
  
  int j = 0;
  int_params[j++] = Debug_Driver;
  int_params[j++] = Test.DDirectory;
  int_params[j++] = Test.Local_Partitions;
  int_params[j++] = Test.Multi_Callbacks;
  int_params[j++] = Test.Null_Import_Lists;
  int_params[j++] = Output.Gnuplot;
  int_params[j++] = Output.Nemesis;
  int_params[j++] = Output.Plot_Partitions;
  int_params[j++] = Output.Mesh_Info_File;
  int_params[j++] = Number_Iterations;
  int_params[j++] = Driver_Action;
  int_params[j++] = Test.Drops;

  MPI_Bcast(int_params, j, MPI_INT, 0, MPI_COMM_WORLD);

  j = 0;
  Debug_Driver           = int_params[j++];
  Test.DDirectory        = int_params[j++];
  Test.Local_Partitions        = int_params[j++];
  Test.Multi_Callbacks   = int_params[j++];
  Test.Null_Import_Lists = int_params[j++];
  Output.Gnuplot         = int_params[j++];
  Output.Nemesis         = int_params[j++];
  Output.Plot_Partitions        = int_params[j++];
  Output.Mesh_Info_File   = int_params[j++];
  Number_Iterations      = int_params[j++];
  Driver_Action          = int_params[j++];
  Test.Drops             = int_params[j++];

  MPI_Bcast(pio_info, sizeof(PARIO_INFO), MPI_BYTE, 0, MPI_COMM_WORLD);

  switch (pio_info->file_type) {
  case CHACO_FILE:
    mesh->data_type = GRAPH;
    break;
  case NEMESIS_FILE:
    mesh->data_type = MESH;
    break;
  case HYPERGRAPH_FILE:
    mesh->data_type = HYPERGRAPH;
    break;
  }

  if(pio_info->dsk_list_cnt > 0) {
    if(Proc != 0)
      pio_info->dsk_list = (int *) malloc(pio_info->dsk_list_cnt*sizeof(int));

    MPI_Bcast(pio_info->dsk_list, pio_info->dsk_list_cnt, MPI_INT,
              0, MPI_COMM_WORLD);
  }

  /* and broadcast the problem specifications */
  MPI_Bcast(prob, sizeof(PROB_INFO), MPI_BYTE, 0, MPI_COMM_WORLD);
  if (prob->num_params > 0) {
    size = prob->num_params * sizeof(Parameter_Pair);
    if (Proc != 0)
      prob->params = (Parameter_Pair *) malloc(size);
    MPI_Bcast(prob->params, size, MPI_CHAR, 0, MPI_COMM_WORLD);
  }

  /* now calculate where the file for this processor is */
  if(pio_info->dsk_list_cnt <= 0) {
    if (pio_info->num_dsk_ctrlrs > 0) {
      ctrl_id = (Proc % pio_info->num_dsk_ctrlrs);
      pio_info->rdisk = ctrl_id + pio_info->pdsk_add_fact;
    }
  }
  else {
    ctrl_id = Proc % pio_info->dsk_list_cnt;
    pio_info->rdisk = pio_info->dsk_list[ctrl_id];
  }

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void gen_par_filename(char *scalar_fname, char *par_fname,
                      PARIO_INFO_PTR pio_info, int proc_for, int nprocs)
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

  int i1, iTemp1;
  int iMaxDigit=0, iMyDigit=0;
  char cTemp[FILENAME_MAX];

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
  if (pio_info->num_dsk_ctrlrs > 0) {
    if(pio_info->zeros) {
      if(pio_info->rdisk <= 9) {
        sprintf(par_fname, "%s%d%d/%s%s", pio_info->pdsk_root,0,
                pio_info->rdisk, pio_info->pdsk_subdir, cTemp);
      }
      else {
        sprintf(par_fname, "%s%d/%s%s", pio_info->pdsk_root,
                pio_info->rdisk, pio_info->pdsk_subdir, cTemp);
      }
    }
    else {
      sprintf(par_fname, "%s%d/%s%s", pio_info->pdsk_root, pio_info->rdisk,
              pio_info->pdsk_subdir, cTemp);
    }
  }
  else
    sprintf(par_fname, "%s/%s", pio_info->pdsk_root, cTemp);

  return;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
