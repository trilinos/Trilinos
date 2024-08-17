// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "dr_const.h"
#include "dr_externs.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"
#include "ch_init_dist_const.h"
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


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
 const char *filename,            /* The name of the command file. */
 PROB_INFO_PTR prob,
 PARIO_INFO_PTR pio_info,   /* pio_info - parallel I/O information. */
 UNDEFINED_INFO_PTR undef)  /* optional list for unrecognized commands */
{
  FILE *file_cmd;
  char  line[MAX_INPUT_STR_LN + 1], *pline, *pmax;
  char  original_line[MAX_INPUT_STR_LN + 1];  /* preserved upper/lower cases */
  char  string[MAX_INPUT_STR_LN], value[MAX_INPUT_STR_LN];
  char  dummy[MAX_INPUT_STR_LN];
  int   i, n, nv, nread, keepreading;


  /* Open the file */
  if ((file_cmd = fopen (filename, "r")) == NULL) {
      return 0;
  }

  /* Begin parsing the input file */
  strcpy(prob->zoltanParams_file, "");

  if (undef){
    undef->list_size = 0;
  }

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

    else if (sscanf(line, " compression" SKIPEQ "%s", value) == 1) {
      if (!strcmp(value, "uncompressed")) {
	pio_info->file_comp = STANDARD;
	continue;
      }
      if (!strcmp(value, "gzip")) {
	pio_info->file_comp = GZIP;
	continue;
      }
    }
    else if (sscanf(line, " zoltan memory debug level" SKIPEQ "%s", value) == 1) {
      if (strcmp(value, "1")) {
	Zoltan_Memory_Debug(1);
      }
      else if (strcmp(value, "2")) {
	Zoltan_Memory_Debug(2);
      }
      else {
	Zoltan_Memory_Debug(3);  /* highest level */
      }
      continue;
    }


    else if (sscanf(line, " file type" LASTARG "%n", value, &n) == 1) {
      if ((!strcmp(value, "chaco")) || (!strcmp(value, "graph")))  {
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
	    else if (!strcmp(value, "owner"))   i = INITIAL_OWNER;
	    else  {
	      Gen_Error(0, "fatal: bad initial distribution argument");
	      return 0;
	    }
	    pio_info->init_dist_type = i;
	  }
	  else if (!strcmp(string, "procs"))  {
	    if (sscanf(value, " %d%n", &pio_info->init_dist_procs, &nv) != 1) {
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
      else if ((strcmp(value, "hypergraph") == 0)
	       || (strcmp(value, "matrixmarket") == 0)
	       || (strcmp(value, "matrixmarket+") == 0)) {
	if (strcmp(value, "hypergraph") == 0){
	  pio_info->file_type       = HYPERGRAPH_FILE;
	  pio_info->init_dist_type  = INITIAL_LINEAR;
	  pio_info->init_dist_pins  = INITIAL_ROW;
	}
	else if (strcmp(value, "matrixmarket") == 0){
	  pio_info->file_type       = MATRIXMARKET_FILE;
	  pio_info->init_dist_type  = INITIAL_LINEAR;
	  pio_info->init_dist_pins  = INITIAL_ROW;
	}
	else if (strcmp(value, "matrixmarket+") == 0){
	  pio_info->file_type       = MATRIXMARKET_PLUS_FILE;
	  pio_info->init_dist_type  = INITIAL_LINEAR;
	  pio_info->init_dist_pins  = INITIAL_ROW;
	}

	pio_info->init_dist_procs = -1;
	pio_info->matrix_obj = COLUMNS;

	pline = line;
	keepreading = 1; /* dummy value to enter loop */

	while ((pline+n < pmax) && keepreading)  {
	  /* Check for options starting with "initial" */
	  keepreading = 0; /* only keep reading if we found new options */
	  pline += n; /* update pline based on previous token */
	  n = 0;      /* no. of chars read. set to zero in case sscanf fails */
	  nread = sscanf(pline, SKIPW "initial" NEXTARG LASTARG "%n", string,
	   value, &n);
	  keepreading += nread;
	  if (nread == 2){
	    if (!strcmp(string, "read")) {
	      if ((value[0] == 'c') && (value[1] == 'h') &&
		  (value[2] == 'u') &&(value[3] == 'n') &&(value[4] == 'k')){
		pio_info->chunk_reader = 1;
	      }
	    }
	    else if (!strcmp(string, "distribution")) {
	      if      (!strcmp(value, "linear"))  i = INITIAL_LINEAR;
	      else if (!strcmp(value, "block"))   i = INITIAL_LINEAR;
	      else if (!strcmp(value, "cyclic"))  i = INITIAL_CYCLIC;
	      else if (!strcmp(value, "owner"))   i = INITIAL_OWNER;
	      else if (!strcmp(value, "file"))    i = INITIAL_FILE;
	      else  {
		Gen_Error(0, "fatal: bad initial distribution argument");
		return 0;
	      }
	      pio_info->init_dist_type = i;
	    }
	    else if (!strcmp(string, "pins")) {
	      if      (!strcmp(value, "linear"))  i = INITIAL_LINEAR;
	      else if (!strcmp(value, "block"))   i = INITIAL_LINEAR;
	      else if (!strcmp(value, "cyclic"))  i = INITIAL_CYCLIC;
	      else if (!strcmp(value, "file"))    i = INITIAL_FILE;
	      else if (!strcmp(value, "zero"))    i = INITIAL_ZERO;
	      else if (!strcmp(value, "row"))     i = INITIAL_ROW;
	      else if (!strcmp(value, "rows"))    i = INITIAL_ROW;
	      else if (!strcmp(value, "column"))  i = INITIAL_COL;
	      else if (!strcmp(value, "columns"))  i = INITIAL_COL;
	      else if (!strcmp(value, "col"))     i = INITIAL_COL;
	      else if (!strcmp(value, "cols"))    i = INITIAL_COL;
	      else  {
		Gen_Error(0, "fatal: bad initial pins argument");
		return 0;
	      }
	      pio_info->init_dist_pins = i;
	    }
	    else if (!strcmp(string, "procs"))  {
	      if (sscanf(value, " %d%n", &pio_info->init_dist_procs, &nv) != 1){
		Gen_Error(0, "fatal: initial procs value must be integal");
		return 0;
	      }
	    }
	    else {
	      Gen_Error(0, "fatal: unrecognizable file type arguments");
	      return 0;
	    }
	  }

	  /* Check for options without "initial" */
	  pline += n; /* update pline based on previous token */
	  n = 0;      /* set to zero in case sscanf fails */
	  i = 0;
	  if (pline < pmax){
	    nread = sscanf(pline, NEXTARG LASTARG "%n", string,
	     value, &n);
	    keepreading += nread;
	    if (nread == 2){
	      if (!strcmp(string, "objects")) {
		if      (!strcmp(value, "rows"))     i = ROWS;
		else if (!strcmp(value, "columns"))  i = COLUMNS;
		else if (!strcmp(value, "nonzeros")) i = NONZEROS;
		else  {
		  Gen_Error(0, "fatal: bad objects type argument");
		  return 0;
		}
		pio_info->matrix_obj = i;
		if (pio_info->matrix_obj != COLUMNS &&
		    pio_info->file_type == MATRIXMARKET_PLUS_FILE) {
		  Gen_Error(0, "fatal: objects != COLUMNS not supported "
			       "for MATRIXMARKET_PLUS");
		  return 0;
		}
	      }
	      else {
		Gen_Error(0, "fatal: unrecognizable file type arguments");
		return 0;
	      }
	    }
	  }
	}
      }
      else if (strcmp(value, "create-a-graph") == 0){
        pio_info->file_type = NO_FILE_GRAPH;         /* zdrive creates a graph */
	strcpy(pio_info->pexo_fname, "zdrive-created-graph");
	pio_info->init_dist_type = INITIAL_NO_DIST;  /* it's already distributed */

	pio_info->init_size     = 10000;       /* default */
	pio_info->init_dim      = 3;           /* default */
	pio_info->init_vwgt_dim = 1;           /* default */

	pline = line;
	while (pline+n < pmax &&
	       sscanf(pline += n, NEXTARG LASTARG "%n", string, value, &n)==2) {
	  if (!strcmp(string, "dimension")
	      && sscanf(value, "%d%n", &pio_info->init_dim, &nv) == 1)
	    continue;
	  else if (!strcmp(string, "obj_weight_dim")
		   && sscanf(value, "%d%n", &pio_info->init_vwgt_dim, &nv) == 1)
	    continue;
	  else if (!strcmp(string, "size")
		   && sscanf(value, ZOLTAN_ID_SPEC "%n", &pio_info->init_size, &nv) == 1)
	    continue;
	  else  {
	    Gen_Error(0, "fatal: bad create-a-graph file parameters");
	    return 0;
	  }
	}
      }
      else if (strcmp(value, "nemesisi") == 0)  {
	pio_info->file_type      = NEMESIS_FILE;
	pio_info->init_dist_type = INITIAL_FILE;
      }
      else if ((strcmp(value, "random-triangles") == 0) || (strcmp(value, "random") == 0))  {
	/* No input file; generate random coordinates. */
	if (strcmp(value, "random-triangles") == 0){
	  pio_info->file_type = NO_FILE_TRIANGLES;
	  strcpy(pio_info->pexo_fname, "random-triangles");
	}
	else{
	  pio_info->file_type       = NO_FILE_POINTS;
	  strcpy(pio_info->pexo_fname, "random");
	}
	pio_info->init_dist_type  = INITIAL_LINEAR;
	pio_info->init_dist_procs = -1;

	pio_info->init_size     = 100;       /* default */
	pio_info->init_dim      = 3;         /* default */
	pio_info->init_vwgt_dim = 1;     /* default */

	pline = line;
	while (pline+n < pmax &&
	       sscanf(pline += n, NEXTARG LASTARG "%n", string, value, &n)==2) {
	  if (!strcmp(string, "dimension")
	      && sscanf(value, "%d%n", &pio_info->init_dim, &nv) == 1)
	    continue;
	  else if (!strcmp(string, "obj_weight_dim")
		   && sscanf(value, "%d%n", &pio_info->init_vwgt_dim, &nv) == 1)
	    continue;
	  else if (!strcmp(string, "size")
		   && sscanf(value,  ZOLTAN_ID_SPEC "%n", &pio_info->init_size, &nv) == 1)
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

    else if (sscanf(line, " text output" SKIPEQ "%d", &Output.Text) == 1)
      continue;                                 /* Generate text output */

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
		    &Output.Plot_Partition, &n) == 1)
      continue;                    /* Plot processor numbers or partitions? */
    else if (sscanf(line, " plot partition" SKIPEQ "%d%n",
		    &Output.Plot_Partition, &n) == 1)
      continue;                    /* Plot processor numbers or partitions? */

    else if (sscanf(line, " print mesh info file" SKIPEQ "%d%n",
		    &Output.Mesh_Info_File, &n) == 1)
      continue;                                /* Generate ASCII mesh file? */

    else if (sscanf(line, " test ddirectory" SKIPEQ "%d%n",
		    &Test.DDirectory, &n) == 1)
      continue;                                  /* DDirectory testing flag */

    else if (sscanf(line, " test drops" SKIPEQ "%d%n", &Test.Drops, &n) == 1)
      continue;                         /* Box- and Point-drop testing flag */

    else if (sscanf (line, " test rcb box" SKIPEQ "%d", &Test.RCB_Box) == 1)
      continue;                         /* Zoltan_RCB_Box testing flag */

    else if (sscanf(line, " test generate files" SKIPEQ "%d%n",
		    &Test.Gen_Files, &n) == 1)
      continue;                             /* file generation testing flag */
    else if (sscanf(line, " test generate file" SKIPEQ "%d%n",
		    &Test.Gen_Files, &n) == 1)
      continue;                             /* file generation testing flag */

    else if (sscanf(line, " test local partitions" SKIPEQ "%d%n",
		    &Test.Local_Parts, &n) == 1)
      continue;                /* Unusual Partition generation testing flag */
    else if (sscanf(line, " test local partition" SKIPEQ "%d%n",
		    &Test.Local_Parts, &n) == 1)
      continue;                /* Unusual Partition generation testing flag */

    else if (sscanf(line, " test fixed objects" SKIPEQ "%d%n",
		    &Test.Fixed_Objects, &n) == 1)
      continue;                /* Fixed objects test flag */
    else if (sscanf(line, " test fixed object" SKIPEQ "%d%n",
		    &Test.Fixed_Objects, &n) == 1)
      continue;                /* Fixed objects test flag */
    else if (sscanf(line, " test dynamic weights" SKIPEQ "%f%n",
		    &Test.Dynamic_Weights, &n) == 1)
      continue;                /* Dynamic weights; changes between iter. */
    else if (sscanf(line, " test dynamic graph" SKIPEQ "%f%n",
		    &Test.Dynamic_Graph, &n) == 1)
      continue;                /* Dynamic graph; edges/verts change between iter. */
    else if (sscanf(line, " test vertex increment" SKIPEQ "%i%n",
                    &Test.Vtx_Inc, &n) == 1)
      continue;                /* Add more vertices in each iteration. */
    else if (sscanf(line, " test multi callbacks" SKIPEQ "%d%n",
		    &Test.Multi_Callbacks, &n) == 1)
      continue;             /* List-based (MULTI) callback function testing */
    else if (sscanf(line, " test multi callback" SKIPEQ "%d%n",
		    &Test.Multi_Callbacks, &n) == 1)
      continue;             /* List-based (MULTI) callback function testing */

    else if (sscanf(line, " test graph callbacks" SKIPEQ "%d%n",
		    &Test.Graph_Callbacks, &n) == 1)
      continue;             /* Graph-based callback function testing */
    else if (sscanf(line, " test graph callback" SKIPEQ "%d%n",
		    &Test.Graph_Callbacks, &n) == 1)
      continue;             /* Graph-based callback function testing */
    else if (sscanf(line, " test hypergraph callbacks" SKIPEQ "%d%n",
		    &Test.Hypergraph_Callbacks, &n) == 1)
      continue;             /* HyperGraph-based callback function testing */
    else if (sscanf(line, " test hypergraph callback" SKIPEQ "%d%n",
		    &Test.Hypergraph_Callbacks, &n) == 1)
      continue;             /* HyperGraph-based callback function testing */
    else if (sscanf(line, " test no global objects" SKIPEQ "%d%n",
		    &Test.No_Global_Objects, &n) == 1)
      continue;             /* HyperGraph-based callback function testing */
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
      continue;            /* zdrive action: Do coloring, load-balancing or ordering? */

    else if (sscanf(line," zdrive debug level" SKIPEQ "%d%n",
		    &Debug_Driver,&n)==1)
      continue;                                /* The Debug reporting level */

    else if (sscanf(line, " zoltan parameter%*[^=]%[=]%n", string, &n) == 1)  {
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

    else if (sscanf(line, " zoltanparams file" SKIPEQ "%s", string) == 1) {
      /* specify zoltanParams-format filename to be read later --
	 included to support hierarchical balancing tests */
      strcpy(prob->zoltanParams_file, string);
    }
    else {
      if (undef){
	if (undef->list_size < UNDEFINED_LIST_MAX){
	  strncpy((char *)(undef->line + undef->list_size), line,
	    UNDEFINED_LENGTH_MAX-1);
	  undef->list_size++;
	}
	else{
	  char buffer[200];
	  sprintf (buffer,
	    "fatal error, too many unrecognized commands: %s\n", line);
	  Gen_Error(0, buffer);
	  return 0;
	}
      }
      else{
	printf("Warning: ignoring unrecognized line in input file:\n\t%s\n",
	  line);
      }
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



/* Broadcast commands to all processes.
   NOTE: This requires manual updating when a new field 
   is added to a struct, like Test!
 */
void brdcst_cmd_info (
  int Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh
)
{
  int ctrl_id, j, k;
  int size;
  int int_params[20];  /* Make sure this array is large enough */
  float float_params[2];  /* Make sure this array is large enough */

  k = 0;
  float_params[k++] = Test.Dynamic_Weights;
  float_params[k++] = Test.Dynamic_Graph;

  MPI_Bcast (float_params, k, MPI_FLOAT, 0, zoltan_get_global_comm());

  k = 0;
  Test.Dynamic_Weights = float_params[k++];
  Test.Dynamic_Graph = float_params[k++];

  j = 0;
  int_params[j++] = Debug_Driver;
  int_params[j++] = Test.DDirectory;
  int_params[j++] = Test.Local_Parts;
  int_params[j++] = Test.Fixed_Objects;
  int_params[j++] = Test.Multi_Callbacks;
  int_params[j++] = Test.Graph_Callbacks;
  int_params[j++] = Test.Hypergraph_Callbacks;
  int_params[j++] = Test.No_Global_Objects;
  int_params[j++] = Test.Null_Lists;
  int_params[j++] = Output.Text;
  int_params[j++] = Output.Gnuplot;
  int_params[j++] = Output.Nemesis;
  int_params[j++] = Output.Plot_Partition;
  int_params[j++] = Output.Mesh_Info_File;
  int_params[j++] = Number_Iterations;
  int_params[j++] = Driver_Action;
  int_params[j++] = Test.Drops;
  int_params[j++] = Test.RCB_Box;
  int_params[j++] = Test.Gen_Files;
  int_params[j++] = Test.Vtx_Inc;

  MPI_Bcast (int_params, j, MPI_INT, 0, zoltan_get_global_comm());

  j = 0;
  Debug_Driver           = int_params[j++];
  Test.DDirectory        = int_params[j++];
  Test.Local_Parts  = int_params[j++];
  Test.Fixed_Objects     = int_params[j++];
  Test.Multi_Callbacks   = int_params[j++];
  Test.Graph_Callbacks   = int_params[j++];
  Test.Hypergraph_Callbacks   = int_params[j++];
  Test.No_Global_Objects = int_params[j++];
  Test.Null_Lists        = int_params[j++];
  Output.Text            = int_params[j++];
  Output.Gnuplot         = int_params[j++];
  Output.Nemesis         = int_params[j++];
  Output.Plot_Partition  = int_params[j++];
  Output.Mesh_Info_File  = int_params[j++];
  Number_Iterations      = int_params[j++];
  Driver_Action          = int_params[j++];
  Test.Drops             = int_params[j++];
  Test.RCB_Box           = int_params[j++];
  Test.Gen_Files         = int_params[j++];
  Test.Vtx_Inc           = int_params[j++];

  MPI_Bcast (pio_info, sizeof(PARIO_INFO), MPI_BYTE, 0, zoltan_get_global_comm());

  switch (pio_info->file_type) {
  case CHACO_FILE:
  case NO_FILE_POINTS:
  case NO_FILE_TRIANGLES:
    mesh->data_type = ZOLTAN_GRAPH;
    break;
  case NEMESIS_FILE:
    mesh->data_type = MESH;
    break;
  case HYPERGRAPH_FILE:
  case MATRIXMARKET_FILE:
  case MATRIXMARKET_PLUS_FILE:
    mesh->data_type = ZOLTAN_HYPERGRAPH;
    break;
  }

  if (pio_info->dsk_list_cnt > 0) {
    if (Proc != 0)
      pio_info->dsk_list = (int*) malloc (pio_info->dsk_list_cnt*sizeof(int));
    MPI_Bcast (pio_info->dsk_list, pio_info->dsk_list_cnt, MPI_INT, 0,
    zoltan_get_global_comm());
  }

  /* and broadcast the problem specifications */
  MPI_Bcast (prob, sizeof(PROB_INFO), MPI_BYTE, 0, zoltan_get_global_comm());
  if (prob->num_params > 0) {
    size = prob->num_params * sizeof(Parameter_Pair);
    if (Proc != 0)
      prob->params = (Parameter_Pair*) malloc(size);
    MPI_Bcast (prob->params, size, MPI_CHAR, 0, zoltan_get_global_comm());
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
  double dnprocs = (double)nprocs;
  double dndisks = (double)pio_info->num_dsk_ctrlrs;

  if (pio_info->num_dsk_ctrlrs <= 0)
    sprintf(par_fname, "%s/%s.%d.%0*d", pio_info->pdsk_root, scalar_fname,
	    nprocs, 1+(int)log10(dnprocs), myproc);
  else if (pio_info->zeros)
    sprintf(par_fname, "%s%0*d/%s%s.%d.%0*d", pio_info->pdsk_root,
	    pio_info->num_dsk_ctrlrs<9 ? 2 : 1+(int)log10(dndisks),
	    pio_info->rdisk, pio_info->pdsk_subdir, scalar_fname, nprocs,
	    1 + (int) log10(dnprocs), myproc);
  else
    sprintf(par_fname, "%s%d/%s/%s.%d.%0*d", pio_info->pdsk_root,
	    pio_info->rdisk, pio_info->pdsk_subdir, scalar_fname, nprocs,
	    1 + (int) log10(dnprocs), myproc);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
