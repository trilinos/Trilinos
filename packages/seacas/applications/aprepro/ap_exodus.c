/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
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

#if !defined(NO_EXODUSII)
#include "my_aprepro.h"
#include "exodusII.h"
#include "ne_nemesisI.h"
#include "netcdf.h"

#include "y.tab.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

extern char *get_temp_filename(void);
void LowerCaseTrim(char *name);

int open_exodus_file(char *filename)
{
  int   cpu = sizeof(double);
  int   io  = 0;
  int   exo;
  float version; 

 exo=ex_open(filename,EX_READ,&cpu,&io,&version);
 if (exo < 0) {
   yyerror("Error opening exodusII file.");
 } else {
  symrec *ptr;
  ptr = putsym("ex_version", VAR, 0);
  ptr->value.var = version;
 }
 return exo;
}

char *do_exodus_info(char *filename)
{
  int exoid;
  int size = 0;
  int count;
  char **info;
  int i;
  char *lines;
  char *ret_string = NULL; 
  
  /*
   * Open the specified exodusII file, read the info records
   * then parse them as input to aprepro.
   */
  exoid = open_exodus_file(filename);
  if (exoid < 0) return "";

  ex_inquire(exoid, EX_INQ_INFO, &count, (float *) NULL, (char *) NULL);
  
  if (count > 0) {
    info = (char**)malloc(count * sizeof(char*));
    for (i=0; i < count; i++) {
      info[i] = (char*)malloc((MAX_LINE_LENGTH+1) * sizeof(char));
      memset(info[i], '\0', MAX_LINE_LENGTH+1);
    }
    
    ex_get_info(exoid, info);
    
    /* Count total size of info records.. */
    for (i=0; i<count; i++) {
      size += strlen(info[i]);  
    }
    size += count-1; /* newlines */
    lines = malloc(size * sizeof(char) + 1);
    lines[0] = '\0';

    for (i=0; i<count; i++) {
      strcat(lines, info[i]);
      strcat(lines, "\n");
    }

    NEWSTR(lines, ret_string);
    if (lines) free(lines);
    
    if (count > 0) {
      for (i=0; i < count; i++) {
	free(info[i]);
      }
      free(info);
    }
    return ret_string;
  } else {
    return "";
  }
  ex_close(exoid);
}

char *do_exodus_meta(char *filename)
{
  int exoid;
  int ndim, nnodes, nelems, nblks, nnsets, nssets;
  char *title;
  symrec *ptr;

  int *ids = NULL;
  
  /*
   * Open the specified exodusII file, read the metadata and set
   * variables for each item.
   * Examples include "node_count", "element_count", ...
   */
  exoid = open_exodus_file(filename);
  if (exoid < 0) return "";

  /* read database paramters */
  title = (char *)calloc ((MAX_LINE_LENGTH+1),sizeof(char *));
  ex_get_init(exoid,title,&ndim,&nnodes,&nelems,&nblks,&nnsets,&nssets);
  
  ptr = putsym("ex_title", SVAR, 0);
  ptr->value.svar = title;

  ptr = putsym("ex_dimension", VAR, 0);
  ptr->value.var = ndim;

  ptr = putsym("ex_node_count", VAR, 0);
  ptr->value.var = nnodes;

  ptr = putsym("ex_element_count", VAR, 0);
  ptr->value.var = nelems;

  ptr = putsym("ex_block_count", VAR, 0);
  ptr->value.var = nblks;

  ptr = putsym("ex_nset_count", VAR, 0);
  ptr->value.var = nnsets;

  ptr = putsym("ex_sset_count", VAR, 0);
  ptr->value.var = nssets;
  
  { /* Nemesis Information */
    int proc_count;
    int proc_in_file;
    char file_type[MAX_STR_LENGTH+1];

    int global_nodes;
    int global_elements;
    int global_blocks;
    int global_nsets;
    int global_ssets;
    int error;
    
    error = ne_get_init_info(exoid, &proc_count, &proc_in_file, file_type);

    if (error >= 0) {
      ptr = putsym("ex_processor_count", VAR, 0);
      ptr->value.var = proc_count;
      
      ne_get_init_global(exoid, &global_nodes, &global_elements,
			 &global_blocks, &global_nsets, &global_ssets);
      
      ptr = putsym("ex_node_count_global", VAR, 0);
      ptr->value.var = global_nodes;
      
      ptr = putsym("ex_element_count_global", VAR, 0);
      ptr->value.var = global_elements;
    }
  }
    
  /*
   * Read The Element Blocks And Set Variables For Those Also.
   * The Scheme Is:
   * -- 'Ex_Block_Ids' Is A List Of Ids.  Due To Aprepro Limitations,
   *     This List Is A String, Not An Integer List...
   * -- Each Block Is Named 'Ex_Block_X' Where X Is Replaced By The
   *    Blocks Position In The List. For Example, The First Block Will
   *    Be Named 'Ex_Block_1'
   *
   * -- Each Block Will Have The Following Symbols:
   *    -- Ex_Block_X_Id = Id Of This Element Block
   *    -- Ex_Block_X_Name = Composed Name "Block_" + Id
   *    -- Ex_Block_X_Element_Count = Number Of Elements In Block
   *    -- Ex_Block_X_Nodes_Per_Element = Number Of Nodes Per Element
   *    -- Ex_Block_X_Topology = Type Of Elements In Block
   *      (Lowercased)
   *    -- Ex_Block_X_Attribute_Count = Number Of Attributes.
   */

  ids = malloc(nblks * sizeof(int));
  ex_get_elem_blk_ids (exoid, ids);

  {
    int i;
    char *buffer = NULL;
    char cid[33];     /* arbitrary size, large enough for INT_MAX */
    int size = 2048;
    char *tmp = NULL;
    
    buffer = calloc(size, sizeof(char));
    for (i=0; i < nblks; i++) {
      sprintf(cid, "%d ", ids[i]);
      if (strlen(buffer) + strlen(cid) +1 > size) {
	if (realloc(buffer, size *=2) == NULL) {
	  yyerror("Error allocating memory.");
	}
	memset(&buffer[size/2], 0, size/2);
      }
      strcat(buffer, cid);
    }
    NEWSTR(buffer, tmp);
    ptr = putsym("ex_block_ids", SVAR, 0);
    ptr->value.svar = tmp;

    if (buffer != NULL)
      free(buffer);
  }
    
  {
    int i;
    char var[128];
    char type[MAX_STR_LENGTH+1];
    char *tmp = NULL;
    int nel;
    int nnel;
    int natr;
    
    for (i=0; i < nblks; i++) {
      ex_get_elem_block(exoid, ids[i], type, &nel, &nnel, &natr);

      sprintf(var, "ex_block_seq_%d_id", i+1);
      ptr = putsym(var, VAR, 0);
      ptr->value.var = ids[i];

      sprintf(var, "ex_block_%d_name", ids[i]);
      ptr = putsym(var, SVAR, 0);
      sprintf(var, "block_%d", ids[i]);
      NEWSTR(var, tmp);
      ptr->value.svar = tmp;

      sprintf(var, "ex_block_%d_element_count", ids[i]);
      ptr = putsym(var, VAR, 0);
      ptr->value.var = nel;

      sprintf(var, "ex_block_%d_nodes_per_element", ids[i]);
      ptr = putsym(var, VAR, 0);
      ptr->value.var = nnel;

      sprintf(var, "ex_block_%d_topology", ids[i]);
      ptr = putsym(var, SVAR, 0);
      NEWSTR(type, tmp);

      /* lowercase the string */
      LowerCaseTrim(tmp);
      ptr->value.svar = tmp;

      sprintf(var, "ex_block_%d_attribute_count", ids[i]);
      ptr = putsym(var, VAR, 0);
      ptr->value.var = natr;
    }
  }
  if (ids != NULL) free(ids);

  {
    /* Get timestep count */
    int ts_count;
    ex_inquire(exoid, EX_INQ_TIME, &ts_count, (float *) NULL, (char *) NULL);
    ptr = putsym("ex_timestep_count", VAR, 0);
    ptr->value.var = ts_count;
    
    if (ts_count > 0) {
      int i;
      symrec *format = getsym("_FORMAT");
      char *buffer = NULL;
      char cid[33];     /* arbitrary size, large enough for double... */
      int size = 2048;
      char *tmp = NULL;
      double *timesteps = malloc(ts_count * sizeof(double));

      ex_get_all_times(exoid, timesteps);

      buffer = calloc(size, sizeof(char));

      for (i=0; i < ts_count; i++) {
	sprintf(cid, format->value.svar, timesteps[i]);
	if (strlen(buffer) + strlen(cid) +2 > size) {
	  if (realloc(buffer, size *=2) == NULL) {
	    yyerror("Error allocating memory.");
	  }
	  memset(&buffer[size/2], 0, size/2);
	}
	strcat(buffer, cid);
	strcat(buffer, " ");
      }
      NEWSTR(buffer, tmp);
      ptr = putsym("ex_timestep_times", SVAR, 0);
      ptr->value.svar = tmp;

      if (buffer != NULL)
	free(buffer);
    }
  }
  
  ex_close(exoid);
  return "";
}

void LowerCaseTrim(char *name)
{
  /*
   * Convert all characters to lowercase. Strip leading whitespace and
   * trim string at first 'whitespace' character (after the leading)
   */

  char *p = name; char *o = name;
  while(*p != '\0' && isspace(*p))
    ++p;
 
  while (*p != '\0') {
    if (isspace(*p)) {
      *o = '\0';
      return;
    } else {
      *o=tolower(*p);
    }
    o++;
    p++;
  }
  *o = '\0';
}
#endif
