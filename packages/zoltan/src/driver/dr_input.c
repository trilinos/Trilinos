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
#include <math.h>
#include <mpi.h>
#include <ctype.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"
#include "ch_init_dist_const.h"


#define SKIPW    "%*[,\t ]"         /* eat up white space amd comma */
#define SKIPEQ   " = "                /* eat up white space and req'd = sign */
#define BIGSKIP  "%*[={()},\t ]"    /* eat up list starts, stops, white space */
#define NEXTARG  "%*[,\t ]%[^=\t ]"  /* argument w/o comma, whitespace */
#define LASTARG  SKIPEQ "%[^,=\t\n ]" /* last arg w/o comma, whitespace */
#define NEXTLIST BIGSKIP "%[^,)}=\t\n ]"
  

/* Purpose: Determine file types for command files and read in the parallel
 * ExodusII command file. Taken from nemesis utilites nem_spread and nem_join. 
 */

int read_cmd_file (
 char *filename,            /* The name of the command file. */
 PROB_INFO_PTR prob,
 PARIO_INFO_PTR pio_info)   /* pio_info - parallel I/O information. */
{
  FILE *file_cmd;
  char  line[MAX_INPUT_STR_LN + 1], *pline, *pmax;
  char  original_line[MAX_INPUT_STR_LN + 1];  /* preserved upper/lower cases */
  char  string[100], value[100], dummy[100];
  int   i, n;

  
  /* Open the file */
  if ((file_cmd = fopen (filename, "r")) == NULL)
    return 0;

  /* Begin parsing the input file */
  prob->num_params = 0;
  prob->params = (Parameter_Pair*) malloc (sizeof(Parameter_Pair));
  while (fgets (line, MAX_INPUT_STR_LN, file_cmd)) {
    strcpy(original_line, line);  /* used when upper/lower case matters */
    pmax = line + strlen (line);
    for (pline = line; pline < pmax; pline++)
      *pline = tolower (*pline);
        
    if (line[0] == '#' || line[0] == '\n')
      continue;                      /* skip blank lines and comment lines */
        
    else if (sscanf(line, " chaco input assignment inverse" SKIPEQ "%d",
      &Chaco_In_Assign_Inv) == 1)       /* ??? Is this used anymore?  */
      continue;

    else if (sscanf(line," decomposition method" SKIPEQ "%s", prob->method)== 1)
      continue;

    else if (sscanf(line, " file name" SKIPEQ "%s", pio_info->pexo_fname) == 1)
      sscanf(original_line, "%*[^=]= %s", pio_info->pexo_fname); /*save case*/

    else if (sscanf(line, " file type" LASTARG "%n", value, &n) == 1) {
      if (!strcmp(value, "chaco"))  {
        pio_info->file_type       = CHACO_FILE;
        pio_info->init_dist_type  = INITIAL_LINEAR;
        pio_info->init_dist_procs = -1;
        pline = line;

        while (pline+n < pmax)  {
          i = sscanf(pline += n, SKIPW "initial" NEXTARG LASTARG "%n", string,
           value, &n);
          if (i != 2)
            break;

          if (!strcmp(string, "distribution")) {
            if      (!strcmp(value, "linear"))  i = INITIAL_LINEAR;
            else if (!strcmp(value, "block"))   i = INITIAL_LINEAR;
            else if (!strcmp(value, "cyclic"))  i = INITIAL_CYCLIC;
            else if (!strcmp(value, "file"))    i = INITIAL_FILE;
            else  {
              Gen_Error(0, "fatal: bad initial distribution argument");
              return 0;
            }
            pio_info->init_dist_type = i;
          }
          else if (!strcmp(string, "procs"))  {
            if (sscanf(value, " %d%n", &pio_info->init_dist_procs, &n) != 1) {
              Gen_Error(0, "fatal: initial procs value must be integal");
              return 0;
            }
          }
          else {
            Gen_Error(0, "fatal: unrecognizable file type arguments");
            return 0;
          }
        }
      }
      else if (strcmp(value, "hypergraph") == 0)  {
        pio_info->file_type       = HYPERGRAPH_FILE;
        pio_info->init_dist_type  = INITIAL_LINEAR;
        pio_info->init_dist_procs = -1;
      }
      else if (strcmp(value, "nemesisi") == 0)  {
        pio_info->file_type      = NEMESIS_FILE;
        pio_info->init_dist_type = INITIAL_FILE;
      }
      else if (strcmp(value, "random") == 0)  {
        /* No input file; generate random coordinates. */
        pio_info->file_type       = NO_FILE;
        pio_info->init_dist_type  = INITIAL_LINEAR;
        pio_info->init_dist_procs = -1;

        pio_info->init_size     = 100;       /* default */
        pio_info->init_dim      = 3;         /* default */
        pio_info->init_vwgt_dim = 1;     /* default */

        strcpy(pio_info->pexo_fname, "random");
        while (pline+n < pmax && 
               sscanf(pline += n, NEXTARG LASTARG "%n", string, value, &n)==2) {
          if (!strcmp(string, "dimension")
              && sscanf(value, "%d%n", &pio_info->init_dim, &n) == 1)
            continue;
          else if (!strcmp(string, "obj_weight_dim")
                   && sscanf(value, "%d%n", &pio_info->init_vwgt_dim, &n) == 1)
            continue;
          else if (!strcmp(string, "size")
                   && sscanf(value, "%d%n", &pio_info->init_size, &n) == 1)
            continue;
          else  {
            Gen_Error(0, "fatal: bad file type = random file parameters");
            return 0;
          }
        }
      }
      else  {
        Gen_Error(0, "fatal: bad file type parameter");
        return 0;
      }
    }

    else if (sscanf(line, " gnuplot output" SKIPEQ "%d", &Output.Gnuplot) == 1)
      continue;                                 /* Generate GNUplot output */

    else if (sscanf(line, " nemesis output" SKIPEQ "%d", &Output.Nemesis) == 1)
      continue;                                 /* Generate Nemesis output */

    else if (sscanf(line, " number of iterations" SKIPEQ "%d",
      &Number_Iterations) == 1)
      continue;      /* The Number of iterations of the balancer to perform */

    else if (sscanf(line, " parallel disk info %[= \t\n]%n", dummy, &n) == 1) {
      pline = line;
      if (sscanf(pline += n, " number = %d%n", &i, &n) != 1)  {
        Gen_Error(0, "fatal: First sup-option for disk info must be number");
        return 0;
      }
      if (i < 0)  {
        Gen_Error(0, "fatal: Invalid value for # of raid controllers.");
        return 0;
      }
      pio_info->num_dsk_ctrlrs = i;

      /* if number_dsk_ctrlrs = 0, then the input file(s) is in the root */
      /* directory given by the parallel disk info, or it is in the same */
      /* directory as the executable if nothing is given for the root    */
      /* infomation. So, no other options can be given when = 0          */
      if (pio_info->num_dsk_ctrlrs == 0)
        continue;

      while (pline+n < pmax) {
        i = sscanf(pline += n, NEXTLIST "%n", string, &n) ;
        if (i != 1)
          break;

        if (strcmp(string, "zeros") == 0)
          pio_info->zeros = 1;
        else if (strcmp(string, "offset") == 0) {
          if (sscanf(pline += n, " =%d%n", &pio_info->pdsk_add_fact, &n) != 1
              || pio_info->pdsk_add_fact < 0)  {
            Gen_Error(0, "fatal: Invalid value for offset.");
            return 0;
          }
        }
        else if (strcmp(string, "list") == 0)  {
          pio_info->dsk_list_cnt = pio_info->num_dsk_ctrlrs;
          pio_info->num_dsk_ctrlrs = -1;

          /* allocate memory for to hold the values */
          pio_info->dsk_list=(int*)malloc(pio_info->dsk_list_cnt*sizeof(int));

          for (i = 0; i < pio_info->dsk_list_cnt; i++)
            if (pline+n < pmax && sscanf(pline += n,  BIGSKIP "%d%n",
             &pio_info->dsk_list[i], &n) == 1)
              continue;
            else {
              Gen_Error(0, "Unknown parameter for parallel disk information");
              return 0;
            }
        }
      }
    }

    else if (sscanf(line," parallel file location %[=]%n", dummy, &n)==1) {
      pline = line;
      while (pline+n < pmax)  {
        i = sscanf(pline += n, NEXTARG LASTARG "%n", string, value, &n);
        if (i != 2)
          break;
        sscanf(original_line + (pline-line), NEXTARG LASTARG, dummy, value); /* reread value from orig line to preserve case */
                       
        if (strcmp(string, "root") == 0)
          strcpy(pio_info->pdsk_root, value);
        if (strcmp(string, "subdir") == 0)   {
          strcpy(pio_info->pdsk_subdir, value);
          if (value [strlen(value)-1] != '/')  {
            pio_info->pdsk_subdir [strlen(value)] = '/';
            pio_info->pdsk_subdir [strlen(value)+1] = 0;
          }
        }
      }
    }

    else if (sscanf(line, " plot partitions" SKIPEQ "%d%n",
                    &Output.Plot_Partitions, &n) == 1)
      continue;                    /* Plot processor numbers or partitions? */
    else if (sscanf(line, " plot partition" SKIPEQ "%d%n",
                    &Output.Plot_Partitions, &n) == 1)
      continue;                    /* Plot processor numbers or partitions? */
      
    else if (sscanf(line, " print mesh info file" SKIPEQ "%d%n",
                    &Output.Mesh_Info_File, &n) == 1)
      continue;                                /* Generate ASCII mesh file? */

    else if (sscanf(line, " test ddirectory" SKIPEQ "%d%n", 
                    &Test.DDirectory, &n) == 1)
      continue;                                  /* DDirectory testing flag */

    else if (sscanf(line, " test drops" SKIPEQ "%d%n", &Test.Drops, &n) == 1)
      continue;                         /* Box- and Point-drop testing flag */

    else if (sscanf(line, " test generate files" SKIPEQ "%d%n",
                    &Test.Gen_Files, &n) == 1)
      continue;                             /* file generation testing flag */
    else if (sscanf(line, " test generate file" SKIPEQ "%d%n",
                    &Test.Gen_Files, &n) == 1)
      continue;                             /* file generation testing flag */      

    else if (sscanf(line, " test local partitions" SKIPEQ "%d%n",
                    &Test.Local_Partitions, &n) == 1)
      continue;                /* Unusual Partition generation testing flag */
    else if (sscanf(line, " test local partition" SKIPEQ "%d%n",
                    &Test.Local_Partitions, &n) == 1)
      continue;                /* Unusual Partition generation testing flag */

    else if (sscanf(line, " test multi callbacks" SKIPEQ "%d%n",
                    &Test.Multi_Callbacks, &n) == 1)
      continue;             /* List-based (MULTI) callback function testing */
    else if (sscanf(line, " test multi callback" SKIPEQ "%d%n",
                    &Test.Multi_Callbacks, &n) == 1)
      continue;             /* List-based (MULTI) callback function testing */
      
    else if (sscanf(line, " test null export lists" SKIPEQ "%d%n", &i,&n)==1) {
      if (i == 1)              /* Null export lists to Help_Migrate testing */
        Test.Null_Lists = EXPORT_LISTS;
      }
    else if (sscanf(line, " test null export list" SKIPEQ "%d%n", &i,&n) == 1) {
      if (i == 1)              /* Null export lists to Help_Migrate testing */
        Test.Null_Lists = EXPORT_LISTS;        
      }

    else if (sscanf(line, " test null import lists" SKIPEQ "%d%n", &i,&n) == 1){
      if (i == 1)              /* Null import lists to Help_Migrate testing */
        Test.Null_Lists = IMPORT_LISTS;
      }
    else if (sscanf(line, " test null import list" SKIPEQ "%d%n", &i,&n) == 1) {
      if (i == 1)              /* Null import lists to Help_Migrate testing */
        Test.Null_Lists = IMPORT_LISTS;
    }

    else if (sscanf(line," zdrive action" SKIPEQ "%d%n",&Driver_Action,&n) == 1)
      continue;            /* zdrive action: Do load-balancing or ordering? */

    else if (sscanf(line," zdrive debug level" SKIPEQ "%d%n",
                    &Debug_Driver,&n)==1)
      continue;                                /* The Debug reporting level */

    else if (sscanf(line, " zoltan parameter %*[s\t ]%[=]%n", dummy, &n)==1) {
      pline = line;
      while (pline+n < pmax && sscanf(pline += n, NEXTARG LASTARG "%n",
                                   prob->params[prob->num_params].Name,
                                   prob->params[prob->num_params].Val, &n)==2) {
        prob->params[prob->num_params++].Index = -1;
        prob->params = (Parameter_Pair*) realloc(prob->params,
                       (prob->num_params+1) * sizeof(Parameter_Pair));
        if (prob->params == NULL)  {
          Gen_Error(0, "fatal: realloc failed for Zoltan Parameters");
          return 0;
        }
      }
    }

    else if (sscanf(line," zoltan vector parameters" LASTARG "%n",string,&n)==1
       ||    sscanf(line," zoltan vector parameter"  LASTARG "%n",string,&n)==1
       ||    sscanf(line," zoltan parameter vectors" LASTARG "%n",string,&n)==1
       ||    sscanf(line," zoltan parameter vector"  LASTARG "%n",string,&n)==1)
    {
      pline = line;
      i = 0;
      while (pline+n < pmax && sscanf(pline += n, BIGSKIP "%[^,\t\n) ]%n",
                               prob->params[prob->num_params].Val, &n)==1) {
        prob->params[prob->num_params].Index = i++;
        strcpy(prob->params[prob->num_params++].Name, string);
          
        prob->params = (Parameter_Pair*) realloc(prob->params,
                       (prob->num_params+1) * sizeof(Parameter_Pair));
        if (prob->params == NULL)  {
          Gen_Error(0, "fatal, realloc failed for Zoltan Parameters");
          return 0;
        }
      }             
    }

    else {
      char buffer[200];
      sprintf (buffer, "fatal error, unrecognized command line: %s\n", line);
      Gen_Error(0, buffer);
      return 0;
    }
  }

  if (prob->num_params == 0) {
    free (prob->params);
    prob->params = NULL;
  }
     
  fclose (file_cmd);
  return 1;
}



int check_inp (
  PROB_INFO_PTR prob, 
  PARIO_INFO_PTR pio_info
)
{
  /* check for the parallel Nemesis file for proc 0 */
  if (strlen(pio_info->pexo_fname) <= 0) {
    Gen_Error (0, "fatal: must specify file base name");
    return 0;
  }

  /* default file type is nemesis */
  if (pio_info->file_type < 0)
    pio_info->file_type = NEMESIS_FILE;

#ifndef ZOLTAN_NEMESIS
  /* if not compiling with the ZOLTAN_NEMESIS flag (i.e., not linking with */
  /* Nemesis library), can't use NEMESIS_FILE file type. */
  if (pio_info->file_type == NEMESIS_FILE) {
    Gen_Error(0, "fatal: must compile for and link with Nemesis libraries for"
                 "Nemesis file types");
    return 0;
  }
#endif

  /* check that there is a list of disks, or a number of raids */
  if ((pio_info->dsk_list_cnt <= 0) && (pio_info->num_dsk_ctrlrs < 0))
    pio_info->num_dsk_ctrlrs = 0;    /* default to single directory */

  /* default is not to have preceeding 0's in the disk names */
  if (pio_info->zeros < 0)
    pio_info->zeros = 0;

  /* most systems that we deal with start their files systems with 1 not 0 */
  if (pio_info->pdsk_add_fact < 0)
    pio_info->pdsk_add_fact = 1;

  /* if there are parallel disks, then specify the root and subdir locations */
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

  /* Check Zoltan specs to insure a load-balancing method is selected */
  if (strlen(prob->method) == 0) {
    Gen_Error(0, "fatal: load balance method must be specified");
    return 0;
  }

  return 1;
}



void brdcst_cmd_info (
  int Proc, 
  PROB_INFO_PTR prob, 
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh
)
{
  int ctrl_id;
  int size;
  int int_params[13];  /* Make sure this array is large enough */

  int j = 0;
  int_params[j++] = Debug_Driver;
  int_params[j++] = Test.DDirectory;
  int_params[j++] = Test.Local_Partitions;
  int_params[j++] = Test.Multi_Callbacks;
  int_params[j++] = Test.Null_Lists;
  int_params[j++] = Output.Gnuplot;
  int_params[j++] = Output.Nemesis;
  int_params[j++] = Output.Plot_Partitions;
  int_params[j++] = Output.Mesh_Info_File;
  int_params[j++] = Number_Iterations;
  int_params[j++] = Driver_Action;
  int_params[j++] = Test.Drops;
  int_params[j++] = Test.Gen_Files;

  MPI_Bcast (int_params, j, MPI_INT, 0, MPI_COMM_WORLD);

  j = 0;
  Debug_Driver           = int_params[j++];
  Test.DDirectory        = int_params[j++];
  Test.Local_Partitions  = int_params[j++];
  Test.Multi_Callbacks   = int_params[j++];
  Test.Null_Lists        = int_params[j++];
  Output.Gnuplot         = int_params[j++];
  Output.Nemesis         = int_params[j++];
  Output.Plot_Partitions = int_params[j++];
  Output.Mesh_Info_File  = int_params[j++];
  Number_Iterations      = int_params[j++];
  Driver_Action          = int_params[j++];
  Test.Drops             = int_params[j++];
  Test.Gen_Files         = int_params[j++];

  MPI_Bcast (pio_info, sizeof(PARIO_INFO), MPI_BYTE, 0, MPI_COMM_WORLD);

  switch (pio_info->file_type) {
  case CHACO_FILE:
  case NO_FILE:
    mesh->data_type = GRAPH;
    break;
  case NEMESIS_FILE:
    mesh->data_type = MESH;
    break;
  case HYPERGRAPH_FILE:
    mesh->data_type = HYPERGRAPH;
    break;
  }

  if (pio_info->dsk_list_cnt > 0) {
    if (Proc != 0)
      pio_info->dsk_list = (int*) malloc (pio_info->dsk_list_cnt*sizeof(int));
    MPI_Bcast (pio_info->dsk_list, pio_info->dsk_list_cnt, MPI_INT, 0,
    MPI_COMM_WORLD);
  }

  /* and broadcast the problem specifications */
  MPI_Bcast (prob, sizeof(PROB_INFO), MPI_BYTE, 0, MPI_COMM_WORLD);
  if (prob->num_params > 0) {
    size = prob->num_params * sizeof(Parameter_Pair);
    if (Proc != 0)
      prob->params = (Parameter_Pair*) malloc(size);
    MPI_Bcast (prob->params, size, MPI_CHAR, 0, MPI_COMM_WORLD);
  }

  /* now calculate where the file for this processor is */
  if (pio_info->dsk_list_cnt <= 0) {
    if (pio_info->num_dsk_ctrlrs > 0) {
      ctrl_id = (Proc % pio_info->num_dsk_ctrlrs);
      pio_info->rdisk = ctrl_id + pio_info->pdsk_add_fact;
    }
  }
  else {
    ctrl_id = Proc % pio_info->dsk_list_cnt;
    pio_info->rdisk = pio_info->dsk_list[ctrl_id];
  }
}



/*----------------------------------------------------------------------------
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
 *
 *----------------------------------------------------------------------------
 *  Determine the number of digits needed to specify the processor ID. 
 *  This allows
 *  numbers like 01-99, i.e., prepending zeros to the name to preserve proper
 *  alphabetic sorting of the files. Comments by Gary Hennigan (1421) */

void gen_par_filename (
  char *scalar_fname, 
  char *par_fname,
  PARIO_INFO_PTR pio_info, 
  int myproc, 
  int nprocs
)
{
  if (pio_info->num_dsk_ctrlrs <= 0)
    sprintf(par_fname, "%s/%s.%d.%0*d", pio_info->pdsk_root, scalar_fname,
            nprocs, 1+(int)log10(nprocs), myproc);
  else if (pio_info->zeros)
    sprintf(par_fname, "%s%0*d/%s%s.%d.%0*d", pio_info->pdsk_root,
            pio_info->num_dsk_ctrlrs<9 ? 2 
                                       : 1+(int)log10(pio_info->num_dsk_ctrlrs),
            pio_info->rdisk, pio_info->pdsk_subdir, scalar_fname, nprocs,
            1 + (int) log10(nprocs), myproc);
  else
    sprintf(par_fname, "%s%d/%s/%s.%d.%0*d", pio_info->pdsk_root,
            pio_info->rdisk, pio_info->pdsk_subdir, scalar_fname, nprocs,
            1 + (int) log10(nprocs), myproc);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

