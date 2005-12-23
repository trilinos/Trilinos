/*
** $Id$
**
** Functions to support writing Zoltan parallel hypergraph examples.
**
**   Defines a hypergraph and divide it among processes for
**   testing purposes.  The division can be controlled with
**   settings described in the header file exzoltan.h.
**
**   Also defines query functions for Zoltan.
*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "exzoltan.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NROWS 28
#define NCOLS 72
static int nprocs, rank;

struct _division{
 int row0, row1, col0, col1;   /* my pins */
 int ew0, ew1;                 /* my edge weights */
 int format;    /* compressed rows or compressed columns */
};

static void divide_interval(int from, int to, 
  int *l, int *r, int nprocs, int rank);

/* Very serious hypergraph, a non-space is a pin */

static char *hg[NROWS]={
"                        b                         .$                    ",
"                        4r                       .$$                    ",
"                        4$r                     z$$$                    ",
"                        $$$c                   d. $$                    ",
"                        $$d$$c.eeeeePeee..   z-4-'P$                    ",
"                       J$$$$$dP3---*$*$$ee$$$$$ - $$F                   ",
"                      z$*.P*$ -- .=-         '|  $$$$                   ",
"                    .P    'F- *e-               ||$$$                   ",
"                    -      4  ^                  ^ $$b                  ",
"                           4$%                    ^$$$r                 ",
"                            F                      4$$$                 ",
"                           4                        $$$b                ",
"                .    4$b   P           .ec          $$$$b               ",
"                     $-               4$P--         $$$$$$c             ",
"               ^    4$$$r 4-          4$eec        4$$$$P-*$e           ",
"     |         .     $$$   4          ^$$$F        $$$$$E.zd$$$-        ",
"      ^|       F            r                     d$$P-  *$$$P .^       ",
"         -.    .        %   ^.                   J$P      $$F.      |   ",
"              -$|    .z.      -                 $P- ^    .P      .^     ",
"                  --3$$$$$$e    -.           .d-     ^   $F  --         ",
"                ^   $$$$$$$$L       ^--==-----      |   ^               ",
"                 -  3$ee= -$-                      |  .-                ",
"                  ^. ^****-                      z-  ^                  ",
"                                               z-  |                    ",
"                      -.                   .e$-  =                      ",
"                        ^ ec..... ....ee$$*-  .^                        ",
"                         -.  ---**---     .=^                           ",
"                              -  === ---                                "
};

/*
** Set options on how the hypergraph is divided up initially.
*/

void *set_initial_division(
  int initial_partitioning,
  int who_has_edge_weights,
  int edge_weights_overlap,
  int format)
{
  int row0, row1, col0, col1;   /* my pins */
  int ew0, ew1;                 /* my edge weights */
  int a, b;
  struct _division *div = NULL;

  if ((initial_partitioning != COMPLETE_ROWS) &&
      (initial_partitioning != COMPLETE_COLS) &&
      (initial_partitioning != BLOCKS)){
    printf("bad initial partitioning value\n");
    return NULL;
  }
  if ((who_has_edge_weights != OMIT_EDGE_WEIGHTS) &&
      (who_has_edge_weights != ALL_HAVE_EDGE_WEIGHTS) &&
      (who_has_edge_weights != SOME_HAVE_EDGE_WEIGHTS) &&
      (who_has_edge_weights != ONE_HAS_EDGE_WEIGHTS)){
    printf("bad who has edge weights\n");
    return NULL;
  } 
  if ((format != ZOLTAN_COMPRESSED_ROWS) &&
      (format != ZOLTAN_COMPRESSED_COLS)){
    printf("bad format\n");
    return NULL;
  }

  div = (struct _division *)malloc(sizeof(struct _division));
  div->format = format;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((nprocs > NROWS) || (nprocs > NCOLS)){
    printf("sorry, we don't handle so many processes\n");
    return NULL;
  }

  if ((nprocs < 3) && (initial_partitioning == BLOCKS)){
    initial_partitioning = COMPLETE_ROWS;
  }

  if (initial_partitioning == COMPLETE_ROWS){
    col0 = 0; col1 = NCOLS-1;
    divide_interval(0, NROWS, &row0, &row1, nprocs, rank);
  }
  else if (initial_partitioning == COMPLETE_COLS){
    row0 = 0; row1 = NROWS-1;
    divide_interval(0, NCOLS, &col0, &col1, nprocs, rank);
  }
  else{
    a = nprocs / 2;
    b = NROWS / 2;
    if (rank < a){
      row0 = 0; row1 = b-1;
      divide_interval(0, NCOLS, &col0, &col1, a, rank);
    }
    else if (rank < 2*a){
      row0 = b; row1 = NROWS-1;
      divide_interval(0, NCOLS, &col0, &col1, a, rank-a);
    }
    else{  /* odd number of processes, last one doesn't get a block */
      col0 = row0 = -1;
      col1 = row1 = -1;
    }
  }
  ew0 = ew1 = -1;
  if (who_has_edge_weights == ALL_HAVE_EDGE_WEIGHTS){
    /* we have edge weights, but not necc. for the same rows of our pins */
    divide_interval(0, NROWS, &ew0, &ew1, nprocs, nprocs-rank);
  }
  else if (who_has_edge_weights == ONE_HAS_EDGE_WEIGHTS){
    if (rank == 0){
      ew0 = 0;  ew1 = NROWS-1;
    }
  }
  else if (who_has_edge_weights == SOME_HAVE_EDGE_WEIGHTS){
    a = nprocs / 2;
    if (rank < a){
      divide_interval(0, NROWS, &ew0, &ew1, a, rank);
    }
  }

  if (edge_weights_overlap){
    if ((ew1 >= 0) && (ew1 < NROWS-1)){
      ew1++;
    }
    if ((ew0 > 1)){
      ew0--;
    }
  }

  div->row0 = row0;
  div->row1 = row1;
  div->col0 = col0;
  div->col1 = col1;
  div->ew0 = ew0;
  div->ew1 = ew1;

  return div;
}

/*
** The query functions
*/
int get_hg_size_and_format(void *data, 
  int *num_lists, int *num_pins, int *format)
{
struct _division *div;
int nedges=0;
int nverts=0;
int r, c;

  div = (struct _division *)data;

  *format = div->format;
  *num_pins = 0;

  if (div->row0 >= 0){
    nedges = div->row1 - div->row0 + 1;
    nverts = div->col1 - div->col0 + 1;

    for (r=div->row0; r<=div->row1; r++){
      for (c=div->col0; c<=div->col1; c++){
        if (!isspace(hg[r][c])){
          *num_pins++;
        }
      }
    }
  }

  if (div->format == ZOLTAN_COMPRESSED_ROWS){
    *num_lists = nedges;
  }
  else{
    *num_lists = nverts;
  }
  return ZOLTAN_OK; 
}

int get_hg(void *data,  int num_gid_entries,
  int nrowcol, int npins, int format,
  ZOLTAN_ID_PTR rowcol_GID, int *rowcol_ptr, ZOLTAN_ID_PTR pin_GID)
{
  int i, j, e, v;
  struct _division *div;

  div = (struct _division *)data;

  if (npins > 0){
    if (format == ZOLTAN_COMPRESSED_ROWS){
      for (i=0, e=div->row0, j=0; e <= div->row1; i++, e++){
        rowcol_ptr[i] = j;
        rowcol_GID[i] = e;

        for (v=div->col0; v <= div->col1; v++){
          if (!isspace(hg[e][v])){
            pin_GID[j++] = v;            
          }
        }
      }
    }
  }

  return ZOLTAN_OK;
}

int get_hg_edge_weight_size(void *data, int *num_edges)
{
  struct _division *div;

  div = (struct _division *)data;

  if (div->ew0 >= 0){
    *num_edges = div->ew1 - div->ew0 + 1;
  }
  return ZOLTAN_OK;
}

int get_hg_edge_weights(void *data,  int num_gid_entries,
  int num_lid_entries, int nedges, int edge_weight_dim,
  ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float *edge_weight)
{
  struct _division *div;
  int i, e;

  div = (struct _division *)data;

  if (div->ew0 >= 0){

    for (i=0, e=div->ew0; e <= div->ew1; i++, e++){
      edge_GID[i] = e;
      edge_LID[i] = i;
      edge_weight[i] = 1.0;
    }
  }

  return ZOLTAN_OK;
}

static void divide_interval(int from, int to, 
             int *l, int *r, int nprocs, int rank)
{
  int a, b;
  int len = to - from + 1;

  a = len / nprocs;
  b = len % nprocs;

  if (rank < b){
    *l = (a+1) * rank;
    *r = *l + a;
  }
  else{
    *l = (a+1)*b + (rank-b)*a;
    *r = *l + a - 1;
  }

  *l += from;
  *r += from;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

