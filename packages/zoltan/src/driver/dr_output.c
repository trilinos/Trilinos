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

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_par_util_const.h"
#include "dr_err_const.h"
#include "dr_output_const.h"

/* Prototypes */
static void echo_cmd_file(FILE *fp, char *cmd_file);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_distributed_mesh(
     int Proc, 
     int Num_Proc, 
     MESH_INFO_PTR mesh)
{
int i, j, k;
int elem;
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
    printf("\tElement Count:      %d\n", mesh->eb_cnts[i]);
    printf("\tNodes Per Element:  %d\n", mesh->eb_nnodes[i]);
    printf("\tAttrib Per Element: %d\n", mesh->eb_nattrs[i]);
  }

  printf("\nElement connect table, partition, weights and coordinates:\n");
  for (i = 0; i < mesh->num_elems; i++) {
    current_elem = &(mesh->elements[i]);
    printf("%d in part %d (%f):\n", current_elem->globalID, 
           current_elem->my_part, current_elem->cpu_wgt[0]);
    for (j = 0; j < mesh->eb_nnodes[current_elem->elem_blk]; j++) {
      printf("\t%d |", current_elem->connect[j]);
      for (k = 0; k < mesh->num_dims; k++) {
        printf(" %f", current_elem->coord[j][k]);
      }
      printf("\n");
    }
  }

  /* now print the adjacencies */
  printf("\nElement adjacencies:\n");
  printf("elem\tnadj(adj_len)\tadj,proc\n");
  for (i = 0; i < mesh->num_elems; i++) {
    current_elem = &(mesh->elements[i]);
    printf("%d\t", current_elem->globalID);
    printf("%d(%d)\t", current_elem->nadj, current_elem->adj_len);
    for (j = 0; j < current_elem->adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (current_elem->adj[j] == -1) continue;

      if (current_elem->adj_proc[j] == Proc)
        elem = mesh->elements[current_elem->adj[j]].globalID;
      else
        elem = current_elem->adj[j];
      printf("%d,%d ", elem, current_elem->adj_proc[j]);
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
      printf("    %d     %d     %d    %d\n", mesh->ecmap_elemids[k],
           mesh->ecmap_sideids[k], 
           mesh->elements[mesh->ecmap_elemids[k]].globalID,
           mesh->ecmap_neighids[k]);
    }
    offset += mesh->ecmap_cnt[i];
  }

  printf("\nHyperedges\n");
  printf("Number of global hyperedges:  %d\n", mesh->gnhedges);
  printf("Number of local hyperedges:   %d\n", mesh->nhedges);
  for (i = 0; i < mesh->nhedges; i++) {
    printf("Local Hyperedge %d:  (", i);
    for (j = mesh->hindex[i]; j < mesh->hindex[i+1]; j++)
      printf("%d ", mesh->hvertex[j]);
    printf(")\n");
    if (mesh->hewgt_dim) {
      printf("Edge Weights:  ");
      for (j = 0; j < mesh->hewgt_dim; j++) 
        printf("%f ", mesh->hewgts[i*mesh->hewgt_dim + j]);
      printf("\n");
    }
  }
  print_sync_end(Proc, Num_Proc, 1);
}


/*--------------------------------------------------------------------------*/
/* Purpose: Output the new element assignments.                             */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
int output_results(char *cmd_file,
                   char *tag,
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
  char  *yo = "output_results";
  char   par_out_fname[FILENAME_MAX+1], ctemp[FILENAME_MAX+1];

  int   *global_ids = NULL;
  int   *parts = NULL;
  int   *perm = NULL;
  int   *invperm = NULL;
  int   *index = NULL;
  int    i, j;

  FILE  *fp;
/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (mesh->num_elems) {
     global_ids = (int *) malloc(5 * mesh->num_elems * sizeof(int));
     if (!global_ids) {
       Gen_Error(0, "fatal: insufficient memory");
       return 0;
     }
     parts = global_ids + mesh->num_elems;
     perm = parts + mesh->num_elems;
     invperm = perm + mesh->num_elems;
     index = invperm + mesh->num_elems;
  }

  for (i = j = 0; i < mesh->elem_array_len; i++) {
    if (mesh->elements[i].globalID >= 0) {
      global_ids[j] = mesh->elements[i].globalID;
      parts[j] = mesh->elements[i].my_part;
      perm[j] = mesh->elements[i].perm_value;
      invperm[j] = mesh->elements[i].invperm_value;
      index[j] = j;
      j++;
    }
  }

  sort_index(mesh->num_elems, global_ids, index);

  /* generate the parallel filename for this processor */
  strcpy(ctemp, pio_info->pexo_fname);
  strcat(ctemp, ".");
  strcat(ctemp, tag);
  gen_par_filename(ctemp, par_out_fname, pio_info, Proc, Num_Proc);

  fp = fopen(par_out_fname, "w");
  if (Proc == 0) 
    echo_cmd_file(fp, cmd_file);

  fprintf(fp, "Global element ids assigned to processor %d\n", Proc);
  fprintf(fp, "GID\tPart\tPerm\tIPerm\n");
  for (i = 0; i < mesh->num_elems; i++) {
    j = index[i];
    fprintf(fp, "%d\t%d\t%d\t%d\n", global_ids[j], parts[j], perm[j], invperm[j]);
  }

  fclose(fp);
  free(global_ids);

  if (Output.Mesh_Info_File) {

    ELEM_INFO_PTR current_element;
    int total_nodes = 0;
    float *x, *y, *z;
    int k;
    int prev_id;

    for (i = 0; i < mesh->num_elems; i++) {
      total_nodes += mesh->eb_nnodes[mesh->elements[i].elem_blk];
    }
    global_ids = (int *) malloc(2 * total_nodes * sizeof(int));
    index = global_ids + total_nodes;
    x = (float *) calloc(3 * total_nodes,  sizeof(float));
    y = x + total_nodes;
    z = y + total_nodes;

    for (k = 0, i = 0; i < mesh->num_elems; i++) {
      current_element = &(mesh->elements[i]);
      for (j = 0; j < mesh->eb_nnodes[current_element->elem_blk]; j++) {
        global_ids[k] = current_element->connect[j];
        x[k] = current_element->coord[j][0];
        if (mesh->num_dims > 1) 
          y[k] = current_element->coord[j][1];
        if (mesh->num_dims > 2)
          z[k] = current_element->coord[j][2];
        index[k] = k;
        k++;
      }
    }

    sort_index(total_nodes, global_ids, index);

    strcat(par_out_fname, ".mesh");
    fp = fopen(par_out_fname, "w");
    fprintf(fp, "Vertex IDs and coordinates\n");
    prev_id = -1;
    for (k = 0; k < total_nodes; k++) {
      j = index[k];
      if (global_ids[j] == prev_id)
        continue;
      prev_id = global_ids[j];
      fprintf(fp, "  %d  (%e, %e, %e)\n", global_ids[j], x[j], y[j], z[j]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Element connectivity:\n");
    for (i = 0; i < mesh->num_elems; i++) {
      current_element = &(mesh->elements[i]);
      fprintf(fp, "  %d  (", current_element->globalID);
      for (j = 0; j < mesh->eb_nnodes[current_element->elem_blk]; j++) {
        fprintf(fp, "%d  ", current_element->connect[j]);
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

static void echo_cmd_file(FILE *fp, char *cmd_file)
{
/* Routine to echo the input file into the output results (so that
 * we know what conditions were used to produce a given result).
 */
char *yo = "echo_cmd_file";
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
