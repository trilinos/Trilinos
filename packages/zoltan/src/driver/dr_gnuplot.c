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

#include "dr_const.h"
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


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Purpose: Output the element assignments in gnuplot format.               */
/*--------------------------------------------------------------------------*/
int output_gnu(const char *cmd_file,
               const char *tag,
               int Proc,
               int Num_Proc,
               PROB_INFO_PTR prob,
               PARIO_INFO_PTR pio_info,
               MESH_INFO_PTR mesh)
/*
 * For 2D problems, output files that can be read by gnuplot for looking at
 * results.
 * We'll do 3D problems later.
 *
 * One gnuplot file is written for each partition.  
 * When number of partitions == number of processors, there is one file per
 * processor.
 *
 * For Chaco input files, the file written contains coordinates of owned
 * nodes and all nodes in that partition connected to the owned nodes. When
 * drawn "with linespoints", the subdomains are drawn, but lines connecting the
 * subdomains are not drawn.
 *
 * For Nemesis input files, the file written contains the coordinates of
 * each node of owned elements.  When drawn "with lines", the element outlines
 * for each owned element are drawn.
 *
 * In addition, processor 0 writes a gnuplot command file telling gnuplot how
 * to process the individual coordinate files written.  This file can be used
 * with the gnuplot "load" command to simplify generation of the gnuplot.
 */
{
  /* Local declarations. */
  const char  *yo = "output_gnu";
  char   par_out_fname[FILENAME_MAX+1], ctemp[FILENAME_MAX+1];
  ELEM_INFO *current_elem, *nbor_elem;
  int    nbor, num_nodes;
  const char  *datastyle = NULL;
  int    i, j, nelems;
  int    prev_part = -1;
  int    max_part = -1;
  float    locMaxX = INT_MIN;
  float    locMinX = INT_MAX;
  float    locMaxY = INT_MIN;
  float    locMinY = INT_MAX;
  float    globMaxX = INT_MIN;
  float    globMinX = INT_MAX;
  float    globMaxY = INT_MIN;
  float    globMinY = INT_MAX;
  int    gmax_part = Num_Proc-1;
  int    gnum_part = Num_Proc;
  int   *parts = NULL;
  int   *index = NULL;
  int   *elem_index = NULL;
  FILE  *fp = NULL;
/***************************** BEGIN EXECUTION ******************************/

  if(Output.Gnuplot < 0)
  {
    Gen_Error(0,"warning: 'gnuplot output' parameter set to invalid negative value.");
    return 0;
  }

  DEBUG_TRACE_START(Proc, yo);

  if (mesh->num_dims > 2) {
    Gen_Error(0, "warning: cannot generate gnuplot data for 3D problems.");
    DEBUG_TRACE_END(Proc, yo);
    return 0;
  }

  if (mesh->eb_nnodes[0] == 0) {
    /* No coordinate information is available.  */
    Gen_Error(0, "warning: cannot generate gnuplot data when no coordinate"
                 " input is given.");
    DEBUG_TRACE_END(Proc, yo);
    return 0;
  }

  /* 
   * Build arrays of partition number to sort by.  Index and elem_index arrays 
   * will be used even when plotting by processor numbers (for generality), 
   * so build it regardless. 
   */
  nelems = mesh->num_elems - mesh->blank_count;

  if (nelems > 0) {
    parts = (int *) malloc(3 * nelems * sizeof(int));
    index = parts + nelems;
    elem_index = index + nelems;
    for (j = 0, i = 0; i < mesh->elem_array_len; i++) {
      current_elem = &(mesh->elements[i]);
      if (current_elem->globalID >= 0) {

        if (mesh->blank_count && (mesh->blank[i] == 1)) continue;
        
        if (current_elem->my_part > max_part) max_part = current_elem->my_part;
        parts[j] = (Output.Plot_Partition ? current_elem->my_part : Proc);
        index[j] = j;
        elem_index[j] = i;
        j++;
      }
    }
  }
  if (Output.Plot_Partition) {
    /* Sort by partition numbers.  Assumes # parts >= # proc. */
    if (nelems > 0) 
      sort_index(nelems, parts, index);
    MPI_Allreduce(&max_part, &gmax_part, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    gnum_part = gmax_part + 1;
  }

  /* generate the parallel filename for this processor */
  strcpy(ctemp, pio_info->pexo_fname);
  strcat(ctemp, ".");
  strcat(ctemp, tag);
  strcat(ctemp, ".gnu");


  if (pio_info->file_type == CHACO_FILE ||
      pio_info->file_type == NO_FILE_POINTS ||
      pio_info->file_type == NO_FILE_TRIANGLES ||
      pio_info->file_type == HYPERGRAPH_FILE) {
    /* 
     * For each node of Chaco graph, print the coordinates of the node.
     * Then, for each neighboring node on the processor, print the neighbor's
     * coordinates.
     */
    datastyle = "linespoints";
    for (i = 0; i < nelems; i++) {
      current_elem = &(mesh->elements[elem_index[index[i]]]);
      if (parts[index[i]] != prev_part) {
        if (fp != NULL) fclose(fp);
        gen_par_filename(ctemp, par_out_fname, pio_info, 
                         parts[index[i]], Num_Proc);
        fp = fopen(par_out_fname, "w");
        prev_part = parts[index[i]];
      }
    
      /* Include the point itself, so that even if there are no edges,
       * the point will appear.  */
      fprintf(fp, "\n%e %e\n", 
              current_elem->coord[0][0], current_elem->coord[0][1]);

      /* save max and min x/y coords */
      if(current_elem->coord[0][0] < locMinX)
      {
        locMinX = current_elem->coord[0][0];
      }
      if(current_elem->coord[0][0] > locMaxX)
      {
        locMaxX = current_elem->coord[0][0];
      }
      if(current_elem->coord[0][1] < locMinY)
      {
        locMinY = current_elem->coord[0][1];
      }
      if(current_elem->coord[0][1] > locMaxY)
      {
        locMaxY = current_elem->coord[0][1];
      }

      if (Output.Gnuplot>1)
      {

        for (j = 0; j < current_elem->nadj; j++) {
          if (current_elem->adj_proc[j] == Proc) {  /* Nbor is on same proc */
            if (mesh->blank_count && (mesh->blank[current_elem->adj[j]] == 1))
              continue;
            if (!Output.Plot_Partition || 
                mesh->elements[current_elem->adj[j]].my_part == 
                             current_elem->my_part) {  
              /* Not plotting partitions, or nbor is in same partition */
              /* Plot the edge.  Need to include current point and nbor point
               * for each edge. */
              fprintf(fp, "\n%e %e\n", 
                  current_elem->coord[0][0], current_elem->coord[0][1]);
              nbor = current_elem->adj[j];
              nbor_elem = &(mesh->elements[nbor]);
              fprintf(fp, "%e %e\n",
                      nbor_elem->coord[0][0], nbor_elem->coord[0][1]);
            }
          }
        }

      }


    }

    MPI_Reduce(&locMinX,&globMinX,1,MPI_FLOAT,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(&locMinY,&globMinY,1,MPI_FLOAT,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(&locMaxX,&globMaxX,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&locMaxY,&globMaxY,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

  }
  else if (pio_info->file_type == NEMESIS_FILE) { /* Nemesis input file */
    /* 
     *  For each element of Nemesis input file, print the coordinates of its
     *  nodes.  No need to follow neighbors, as decomposition is by elements.
     */
    double sum[2];
    datastyle = "lines";
    for (i = 0; i < nelems; i++) {
      current_elem = &(mesh->elements[elem_index[index[i]]]);
      if (parts[index[i]] != prev_part) {
        if (fp != NULL) fclose(fp);
        gen_par_filename(ctemp, par_out_fname, pio_info, 
                         parts[index[i]], Num_Proc);
        fp = fopen(par_out_fname, "w");
        prev_part = parts[index[i]];
      }
      num_nodes = mesh->eb_nnodes[current_elem->elem_blk];
      sum[0] = sum[1] = 0.0;
      for (j = 0; j < num_nodes; j++) {
        fprintf(fp, "%e %e\n", 
                current_elem->coord[j][0], current_elem->coord[j][1]);
        sum[0] += current_elem->coord[j][0];
        sum[1] += current_elem->coord[j][1];
      }
      fprintf(fp, "%e %e\n", current_elem->coord[0][0], 
                             current_elem->coord[0][1]);
      fprintf(fp, "\n");
      /* Print small + in center of element */
      sum[0] /= num_nodes;
      sum[1] /= num_nodes;
      fprintf(fp, "%e %e\n",   sum[0] - 0.001, sum[1]);
      fprintf(fp, "%e %e\n\n", sum[0] + 0.001, sum[1]);
      fprintf(fp, "%e %e\n",   sum[0], sum[1] - 0.001);
      fprintf(fp, "%e %e\n\n", sum[0], sum[1] + 0.001);
    }
  }
  
  if (nelems == 0 && !Output.Plot_Partition) { 
    /* Open a file just so one exists; satisfies the gnuload file. */
    gen_par_filename(ctemp, par_out_fname, pio_info, Proc, Num_Proc);
    fp = fopen(par_out_fname, "w");
  }
    
  if (fp != NULL) fclose(fp);
  safe_free((void **)(void *) &parts);

  if (Proc == 0) {
    /* Write gnu master file with gnu commands for plotting */
    strcpy(ctemp, pio_info->pexo_fname);
    strcat(ctemp, ".");
    strcat(ctemp, tag);
    strcat(ctemp, ".gnuload");
    fp = fopen(ctemp, "w");
    fprintf(fp, "set nokey\n");
    fprintf(fp, "set nolabel\n");
    fprintf(fp, "set noxzeroaxis\n");
    fprintf(fp, "set noyzeroaxis\n");
    fprintf(fp, "set noxtics\n");
    fprintf(fp, "set noytics\n");
    fprintf(fp, "set data style %s\n", datastyle);

    /* resize range so that there is a 5% border around data */
    fprintf(fp, "set xrange [%f:%f] \n ",globMinX-(globMaxX-globMinX)/20
	                            ,globMaxX+(globMaxX-globMinX)/20);
    fprintf(fp, "set yrange [%f:%f] \n ",globMinY-(globMaxY-globMinY)/20
	                            ,globMaxY+(globMaxY-globMinY)/20);


    fprintf(fp, "plot ");
    strcpy(ctemp, pio_info->pexo_fname);
    strcat(ctemp, ".");
    strcat(ctemp, tag);
    strcat(ctemp, ".gnu");
    for (i = 0; i < gnum_part; i++) {
      gen_par_filename(ctemp, par_out_fname, pio_info, i, Num_Proc);
      fprintf(fp, "\"%s\"", par_out_fname);
      if (i != gnum_part-1) {
        fprintf(fp, ",\\\n");
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
