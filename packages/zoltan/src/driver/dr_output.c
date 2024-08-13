// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "dr_const.h"
#include "dr_externs.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_par_util_const.h"
#include "dr_err_const.h"
#include "dr_output_const.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Prototypes */
static void echo_cmd_file(FILE *fp, const char *cmd_file);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_distributed_mesh(
     int Proc, 
     int Num_Proc, 
     MESH_INFO_PTR mesh)
{
int i, j, k;
ZOLTAN_ID_TYPE elem;
int offset;
ELEM_INFO_PTR current_elem;


/*
 * Print the distributed mesh description for each processor.  This routine
 * is useful for debugging the input meshes (Nemesis or Chaco).  It is 
 * serial, so it should not be used for production runs.
 */

  print_sync_start(Proc, 1);

  printf ("############## Mesh information for processor %d ##############\n",
          Proc);
  printf ("Number Dimensions:     %d\n", mesh->num_dims);
  printf ("Number Nodes:          %d\n", mesh->num_nodes);
  printf ("Number Elements:       %d\n", mesh->num_elems);
  printf ("Number Element Blocks: %d\n", mesh->num_el_blks);
  printf ("Number Node Sets:      %d\n", mesh->num_node_sets);
  printf ("Number Side Sets:      %d\n", mesh->num_side_sets);

  for (i = 0; i < mesh->num_el_blks; i++) {
    printf("\nElement block #%d\n", (i+1));
    printf("\tID:                 %d\n", mesh->eb_ids[i]);
    printf("\tElement Type:       %s\n", mesh->eb_names[i]);
    printf("\tElement Count:      " ZOLTAN_ID_SPEC "\n", mesh->eb_cnts[i]);
    printf("\tNodes Per Element:  %d\n", mesh->eb_nnodes[i]);
    printf("\tAttrib Per Element: %d\n", mesh->eb_nattrs[i]);
  }

  printf("\nElement connect table, partition, weights and coordinates:\n");
  for (i = 0; i < mesh->elem_array_len; i++) {
    current_elem = &(mesh->elements[i]);
    if (current_elem->globalID == ZOLTAN_ID_INVALID) continue;

    printf(ZOLTAN_ID_SPEC " in part %d (%f):\n", current_elem->globalID, 
           current_elem->my_part, current_elem->cpu_wgt[0]);
    for (j = 0; j < mesh->eb_nnodes[current_elem->elem_blk]; j++) {
      printf("\t" ZOLTAN_ID_SPEC " |", current_elem->connect[j]);
      for (k = 0; k < mesh->num_dims; k++) {
        printf(" %f", current_elem->coord[j][k]);
      }
      printf("\n");
    }
  }

  /* now print the adjacencies */
  printf("\nElement adjacencies:\n");
  printf("elem\tnadj(adj_len)\tadj,proc\n");
  for (i = 0; i < mesh->elem_array_len; i++) {
    current_elem = &(mesh->elements[i]);
    if (current_elem->globalID == ZOLTAN_ID_INVALID) continue;

    printf(ZOLTAN_ID_SPEC "\t", current_elem->globalID);
    printf("%d(%d)\t", current_elem->nadj, current_elem->adj_len);
    for (j = 0; j < current_elem->adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (current_elem->adj[j] == ZOLTAN_ID_INVALID) continue;

      if (current_elem->adj_proc[j] == Proc)
        elem = mesh->elements[current_elem->adj[j]].globalID;
      else
        elem = current_elem->adj[j];
      printf(ZOLTAN_ID_SPEC ",%d ", elem, current_elem->adj_proc[j]);
    }
    printf("\n");
  }

  printf("\nCommunication maps\n");
  printf("Number of maps: %d\n", mesh->necmap);
  printf("Map Ids(and Counts):");
  for (i = 0; i < mesh->necmap; i++)
    printf(" %d(%d)", mesh->ecmap_id[i], mesh->ecmap_cnt[i]);
  printf("\n");
  offset = 0;
  for (i = 0; i < mesh->necmap; i++) {
    printf("Map %d:\n", mesh->ecmap_id[i]);
    printf("    elem   side   globalID  neighID\n");
    for (j = 0; j < mesh->ecmap_cnt[i]; j++) {
      k = j + offset;
      printf("    %d     %d     " ZOLTAN_ID_SPEC "    " ZOLTAN_ID_SPEC "\n", mesh->ecmap_elemids[k],
           mesh->ecmap_sideids[k], 
           mesh->elements[mesh->ecmap_elemids[k]].globalID,
           mesh->ecmap_neighids[k]);
    }
    offset += mesh->ecmap_cnt[i];
  }

  printf("\nHyperedges\n");
  printf("Number of global hyperedges:  " ZOLTAN_ID_SPEC "\n", mesh->gnhedges);
  if (mesh->format == ZOLTAN_COMPRESSED_EDGE){
    printf("Number of rows (edges):   %d\n", mesh->nhedges);
  }
  else{
    printf("Number of columns (vertices):   %d\n", mesh->nhedges);
  }
  for (i = 0; i < mesh->nhedges; i++) {
    printf("  " ZOLTAN_ID_SPEC " (%d):  (", mesh->hgid[i], i);
    for (j = mesh->hindex[i]; j < mesh->hindex[i+1]; j++)
      printf(ZOLTAN_ID_SPEC, mesh->hvertex[j]);
    printf(")\n");
  }
  if (mesh->hewgt_dim && (mesh->heNumWgts > 0)){
    printf("\nHyperedge Weights\n");
    for (i=0; i<mesh->heNumWgts; i++){ 
      if (mesh->heWgtId){
        printf("Hyperedge " ZOLTAN_ID_SPEC " (%d):  (", mesh->heWgtId[i], i);
      }
      else{
        printf("Hyperedge " ZOLTAN_ID_SPEC " (%d):  (", mesh->hgid[i], i);
      }
      for (j = 0; j < mesh->hewgt_dim; j++) 
        printf("%f ", mesh->hewgts[i*mesh->hewgt_dim + j]);
      printf(")\n");
    }
  }
  print_sync_end(Proc, Num_Proc, 1);
}


/*--------------------------------------------------------------------------*/
/* Purpose: Output the new element assignments.                             */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
int output_results(const char *cmd_file,
                   const char *tag,
                   int Proc,
                   int Num_Proc,
                   PROB_INFO_PTR prob,
                   PARIO_INFO_PTR pio_info,
                   MESH_INFO_PTR mesh)
/*
 * For the first swipe at this, don't try to create a new
 * exodus/nemesis file or anything. Just get the global ids,
 * sort them, and print them to a new ascii file.
 */
{
  /* Local declarations. */
  const char  *yo = "output_results";
  char   par_out_fname[FILENAME_MAX+1], ctemp[FILENAME_MAX+1];
  char cmsg[256];

  ZOLTAN_ID_TYPE   *global_ids = NULL;
  ZOLTAN_ID_TYPE   *parts = NULL;
  ZOLTAN_ID_TYPE   *perm = NULL;
  ZOLTAN_ID_TYPE   *invperm = NULL;
  int              *index = NULL;
  int    i, j;

  FILE  *fp;
/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (mesh->num_elems) {
     global_ids = (ZOLTAN_ID_TYPE*) malloc(4 * mesh->num_elems * sizeof(ZOLTAN_ID_TYPE) + mesh->num_elems * sizeof(int));
     if (!global_ids) {
       Gen_Error(0, "fatal: insufficient memory");
       return 0;
     }
     parts = global_ids + mesh->num_elems;
     perm = parts + mesh->num_elems;
     invperm = perm + mesh->num_elems;
     index = (int*)(invperm + mesh->num_elems);
  }

  for (i = j = 0; i < mesh->elem_array_len; i++) {
    if (mesh->elements[i].globalID != ZOLTAN_ID_INVALID) {
      global_ids[j] = (ZOLTAN_ID_TYPE) mesh->elements[i].globalID;
      parts[j] =(ZOLTAN_ID_TYPE) mesh->elements[i].my_part;
      perm[j] = (ZOLTAN_ID_TYPE)mesh->elements[i].perm_value;
      invperm[j] = (ZOLTAN_ID_TYPE)mesh->elements[i].invperm_value;
      index[j] = j;
      j++;
    }
  }

  quicksort_pointer_inc_id_id(index, global_ids, NULL,
                              0, mesh->num_elems-1);

  /* generate the parallel filename for this processor */
  strcpy(ctemp, pio_info->pexo_fname);
  strcat(ctemp, ".");
  strcat(ctemp, tag);
  gen_par_filename(ctemp, par_out_fname, pio_info, Proc, Num_Proc);

  fp = fopen(par_out_fname, "w");

  if (fp == NULL){
    sprintf(cmsg, "Error in %s; %s can not be opened for writing.", yo, par_out_fname);
    Gen_Error(0, cmsg);
    return 0;
  }

  if (Proc == 0) 
    echo_cmd_file(fp, cmd_file);

  fprintf(fp, "Global element ids assigned to processor %d\n", Proc);
  fprintf(fp, "GID\tPart\tPerm\tIPerm\n");
  for (i = 0; i < mesh->num_elems; i++) {
    j = index[i];
    fprintf(fp, ZOLTAN_ID_SPEC "\t%d\t%d\t%d\n", 
              global_ids[j], (int)parts[j], (int)perm[j], (int)invperm[j]);
  }

  fclose(fp);
  free(global_ids);

  if (Output.Mesh_Info_File) {

    ELEM_INFO_PTR current_element;
    int total_nodes = 0;
    float *x, *y, *z;
    int k;
    ZOLTAN_ID_TYPE prev_id;

    for (i = 0; i < mesh->num_elems; i++) {
      total_nodes += mesh->eb_nnodes[mesh->elements[i].elem_blk];
    }
    global_ids = (ZOLTAN_ID_TYPE *) malloc(total_nodes * (sizeof(ZOLTAN_ID_TYPE)+sizeof(int)));
    index = (int*)(global_ids + total_nodes);
    x = (float *) calloc(3 * total_nodes,  sizeof(float));
    y = x + total_nodes;
    z = y + total_nodes;

    for (k = 0, i = 0; i < mesh->num_elems; i++) {
      current_element = &(mesh->elements[i]);
      for (j = 0; j < mesh->eb_nnodes[current_element->elem_blk]; j++) {
        global_ids[k] = (ZOLTAN_ID_TYPE)(current_element->connect[j]);
        x[k] = current_element->coord[j][0];
        if (mesh->num_dims > 1) 
          y[k] = current_element->coord[j][1];
        if (mesh->num_dims > 2)
          z[k] = current_element->coord[j][2];
        index[k] = k;
        k++;
      }
    }

    quicksort_pointer_inc_id_id(index, global_ids, NULL, 
                                0, total_nodes-1);

    strcat(par_out_fname, ".mesh");
    fp = fopen(par_out_fname, "w");
    fprintf(fp, "Vertex IDs and coordinates\n");
    prev_id = ZOLTAN_ID_INVALID;
    for (k = 0; k < total_nodes; k++) {
      j = (int)index[k];
      if (global_ids[j] == prev_id)
        continue;
      prev_id = global_ids[j];
      fprintf(fp, "  " ZOLTAN_ID_SPEC "  (%e, %e, %e)\n", global_ids[j], x[j], y[j], z[j]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Element connectivity:\n");
    for (i = 0; i < mesh->num_elems; i++) {
      current_element = &(mesh->elements[i]);
      fprintf(fp, "  " ZOLTAN_ID_SPEC "  (", current_element->globalID);
      for (j = 0; j < mesh->eb_nnodes[current_element->elem_blk]; j++) {
        fprintf(fp, ZOLTAN_ID_SPEC "  ", current_element->connect[j]);
      }
      fprintf(fp, ")\n");
    }
    
    fclose(fp);
    free(global_ids);
    free(x);
  }

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void echo_cmd_file(FILE *fp, const char *cmd_file)
{
/* Routine to echo the input file into the output results (so that
 * we know what conditions were used to produce a given result).
 */
const char *yo = "echo_cmd_file";
char cmsg[256];
char inp_line[MAX_INPUT_STR_LN + 1];
FILE *cmd_fp;

  /* Open the file */
  if((cmd_fp=fopen(cmd_file, "r")) == NULL) {
    sprintf(cmsg, "Error in %s; input file %s does not exist.", yo, cmd_file);
    Gen_Error(0, cmsg);
    return;
  }

  while(fgets(inp_line, MAX_INPUT_STR_LN, cmd_fp)) {
    /* skip any line that is a comment */
    if((inp_line[0] != '#') && (inp_line[0] != '\n')) 
      fprintf(fp, "%s", inp_line);
  }

  fclose(cmd_fp);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
