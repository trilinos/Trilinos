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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elb_const.h"
#include "elb_err_const.h"

/*****************************************************************************/
static void scandescriptor(const char *d, int *blkids, int n, int nblks,
                           PROB_INFO_PTR prob);
static void chgrp(int grp_id, int blk, int *blkids, int nblks,
                  PROB_INFO_PTR prob);
extern int ilog2i (unsigned int n);

/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function parse_groups() begins:
 *----------------------------------------------------------------------------
 * This function will parse the group designation string and set up which
 * element blocks will be in which groups.
 *
 * the group designator prob->groups follows these rules:
 *  - Blocks are grouped using the slash "/" character.
 *  - Ids are separated with white space, comma, or by the hyphen "-" character. *  - Any blocks not included in the list, are added to a separate group.
 *  - Duplicates in the list are permitted, but the last group to which a
 *    block is placed is where the block will go.
 *  - Block IDs not in the exodus file are quietly ignored.
 *
 * Examples.
 *   Assume block IDs= 1-20 31-45
 *
 *     descriptor              group1          group2     group3
 *   - "1-20"                  1-20             31-45
 *   - "30-45 3/ 10-12"        3, 30-45         10,11,12   1,2,4-20
 *   - "1-20/40-45/5-10 21-41   1-4,11-20       42,43,45   5-10 31-41
 *
 *****************************************************************************/
int parse_groups(int *el_blk_ids,
                 int *el_blk_cnts,
                 MESH_INFO_PTR mesh,
                 PROB_INFO_PTR prob
  )
{
  char *id;
  int   i, last, found;

/*---------------------------Execution Begins--------------------------------*/

  /* allocate memory for the groups */
  prob->group_no = (int *) malloc (mesh->num_el_blks * sizeof(int));
  mesh->eb_cnts = (int *) malloc (mesh->num_el_blks * sizeof(int));
  if (!(prob->group_no) || !(mesh->eb_cnts))
  {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /* prepare the group number array, and copy the element block counts */
  for (i = 0; i < mesh->num_el_blks; i++) {
    prob->group_no[i] = -1;
    mesh->eb_cnts[i] = el_blk_cnts[i];
  }

  /* convert any comma's to blank spaces in the designator string */
  for (i = 0; i < strlen(prob->groups); i++)
    if (prob->groups[i] == ',') prob->groups[i] = ' ';

  /* fill in the group identifier for each block */
  id = prob->groups;
  i = 0;
  do {
    if (*id == '/') id++;
    scandescriptor(id, el_blk_ids, i, mesh->num_el_blks, prob);
    id = strchr(id, '/');
    i++;
  } while (id != NULL);
  last = i;

  /* set any remaining blocks to new group */
  found = 0;
  for (i = 0; i < mesh->num_el_blks; i++)
    if (prob->group_no[i] < 0) {
      prob->group_no[i] = last;
      found = 1;
     }

  if (found) last++;

  prob->num_groups = last;

  {
    int first_el = 0;
    printf("\nNumber of blocks: %d\n", mesh->num_el_blks);
    printf("Block ID and associated groups:\n");
    printf("   block   #elems  group   type\n");
    for (i = 0; i < mesh->num_el_blks; i++) {
      printf("%8d%8d%8d%8s\n", el_blk_ids[i], mesh->eb_cnts[i], prob->group_no[i], elem_names[mesh->elem_type[first_el]]);
      first_el += mesh->eb_cnts[i];
    }
    printf("There are %d groups of blocks\n", prob->num_groups);
  }

  /* finnished with the group designator string */
  free (prob->groups);

  return 1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_group_info() begins:
 *----------------------------------------------------------------------------
 * This function will get the information that is needed to break up the
 * groups before passing them to Chaco.
 * It generates the following:
 *   - an array of with the group number for each element
 *   - an array of elements per group
 *   - an array of processors to be used per group
 *   - the max number of vertecies and adjacencies for all of the groups
 *     this value is needed so that arrays used to pass information
 *     to Chaco can be allocated
 *****************************************************************************/
int get_group_info(MACHINE_PTR machine,
                   PROB_INFO_PTR prob,
                   MESH_INFO_PTR mesh,
                   GRAPH_INFO_PTR graph,
                   int elem2grp[],
                   int nprocg[],
                   int nelemg[],
                   int *max_vtx,
                   int *max_adj
  )
{
  int nproc, sum, iblk;
  int i, j;
  int *nadj_per_grp;

/*---------------------------Execution Begins--------------------------------*/

  /* allocate array to hold adjacency counts, if necessary */
  if (prob->alloc_graph == ELB_TRUE) {
    nadj_per_grp = (int *) malloc(prob->num_groups * sizeof(int));
    if (!nadj_per_grp) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
  }

  /* initialize the group counts arrays */
  for (i = 0; i < prob->num_groups; i++) {
    nelemg[i] = 0;
    if (prob->alloc_graph == ELB_TRUE) nadj_per_grp[i] = 0;
  }

  /*
   * fill the vertex2proc array with the group number of the individual
   * elements, calculate how many elements in each group and determine
   * how many adjacencies are in each group (if necessary).
   */
  iblk = sum = 0;
  for (i = 0; i < prob->num_vertices; i++) {

    /* figure out which element block this is */
    if (sum == mesh->eb_cnts[iblk]) {
      sum = 0;
      iblk++;
    }
    sum++;

    /*
     * use negative numbers to specify the groups, in order to
     * avoid having a problem with group 0, add 1 to the group
     * number
     */
    elem2grp[i] =  -(prob->group_no[iblk] + 1);

    nelemg[prob->group_no[iblk]]++;

    if (prob->alloc_graph == ELB_TRUE)
      nadj_per_grp[prob->group_no[iblk]] += graph->start[i+1]
                                            - graph->start[i];
  }

  /*
   * calculate how many processors to use for each group
   *   using method from the materials group, haven't really checked it
   */
  if (machine->type == MESH)
    nproc = machine->procs_per_box;
  else if (machine->type == HCUBE)
    nproc = ilog2i((unsigned int) machine->procs_per_box);
  for (i = 0; i < prob->num_groups; i++) {
    nprocg[i] = (nproc * (nelemg[i] + .5)) / (float) prob->num_vertices;
    if (nelemg[i] && !nprocg[i]) nprocg[i] = 1;
  }

  /*
   * check to see if correct number of processors have been allocated
   * and get the maximum number of vertices
   */
  sum = 0;
  j = 0;
  *max_vtx = 0;
  *max_adj = 0;
  for (i = 0; i < prob->num_groups; i++) {
    sum += nprocg[i];
    if (nprocg[i] > nprocg[j]) {
      j = i;
      *max_vtx = nelemg[j];  /* most processors implies most elements */
    }

    /* determine how large to make temporary arrays */
    if (nelemg[i] > *max_vtx) *max_vtx = nelemg[i];
    if(prob->alloc_graph == ELB_TRUE)
      if (nadj_per_grp[i] > *max_adj) *max_adj = nadj_per_grp[i];
  }

  if (sum != nproc) {
    /* correct group with most processors (j determined above) */
    nprocg[j] -= (sum - nproc);
    if (nprocg[j] <= 0) {
      Gen_Error(0,"Unable to balance # processors in get_group_info().");
      return 0;
    }
  }

  printf("Load balance information\n");
  for (i = 0; i < prob->num_groups; i++)
    printf("group[%d]  #elements=%-10d  #proc=%d\n",i,nelemg[i],nprocg[i]);

  if(prob->alloc_graph == ELB_TRUE)
    free (nadj_per_grp);

  return 1;

}

/**********************************************************************/
static void scandescriptor(const char *d, int *blkids, int n, int nblks,
                           PROB_INFO_PTR prob)
{
  /* reads the descriptor up to the next "/" character. interprets the
     descriptor. The ranges specified in the descriptor are then stored
     in the grp array. */

  const char* p=d;
  int i;        /* integer read from string */
  int last=0;   /* last integer read */
  int stop;     /* stop value in a string range */
  int q;        /* number of ints read */
  int qn;       /* number of bytes read */
  int c;        /* integer index when spanning a range */

  while (*p != '/' && *p != 0) {
    q = sscanf(p, "%d%n", &i, &qn);
    if (q == 0 || i < 0){
      if (p[qn-1] == '/' || *p == 0 ) return;
      else if (i < 0){
        stop = -i;
        for (c = last; c <= stop; c++)
          chgrp(n, c, blkids, nblks, prob);
      }
      else if (p[qn-1] == '-'){
        p += qn;
        sscanf(p, "%d%n", &stop, &qn);
        for (c = last; c <= stop; c++)
          chgrp(n, c, blkids, nblks, prob);
      }
      else {
        /* check for minus sign */
        for (c = 0; c < qn; c++)
          if (p[c] == '-') break;
        if (c < qn){
          p += qn;
          sscanf(p, "%d%n", &stop, &qn);
          for (c = last; c <= stop; c++)
            chgrp(n, c, blkids, nblks, prob);
        }
        else {
          printf("Error reading descriptor '%s'\n", d);
          printf("                          ");
          for (c = 0;c < qn; c++)
            printf(" ");
          printf("^\n");
          return;
        }
      }
    }
    else
      last = i;
    chgrp(n, i, blkids, nblks, prob);
    p += qn;
  }
}


/**********************************************************************/
/* changes the grp[] entry corresponding to block=blk to the value grp_id
 * It uses the extern variables "nblks", "grp" and "blkid". only grp
 * is altered.
 *
 * Not very efficient. To find the blk, it loops through all blks.
 */
static void chgrp(int grp_id, int blk, int *blkids, int nblks,
                  PROB_INFO_PTR prob)
{
  int j;

  for (j = 0; j < nblks; j++)
    if (blkids[j] == blk){
      prob->group_no[j] = grp_id;
      return;
    }
}


/**********************************************************************/
/* builds a map of the element ID for a group. The length of the map is
 * the total number of elements. Elements that are within the group are
 * numbered from 1:n. Elements outside the group are assigned -1.
 *
 * Note that adjacency lists are 1:n, not 0:n-1.
 *
 */
void build_grp_map(MESH_INFO_PTR mesh, PROB_INFO_PTR prob, int* map, int g)
{
  int blk;      /* block index */
  int id = 1;   /* element id */
  int elem = 0; /* element index */
  int k;        /* element index within a block */

  for (blk = 0; blk < mesh->num_el_blks; blk++){
    for (k = 0; k < mesh->eb_cnts[blk]; k++, elem++) {
      if (prob->group_no[blk] == g )
        map[elem] = id++;
      else
        map[elem] = -1;
    }
  }
}

