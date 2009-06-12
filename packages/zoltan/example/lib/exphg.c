/*
** $Id$
**
** Functions to support writing Zoltan parallel hypergraph examples.
**
**   Defines a hypergraph and divides it among processes for
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

/************** sample hypergraphs *******************************************/
#if 0
/* felix: non-space is a pin, has vertices (cols) with no hyperedges, 
 * has dense and non-dense edges, ascii art helps in spotting errors
 */
#define NROWS 28
#define NCOLS 72
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
#endif
/* From Zoltan's test directory: the graph ch_simple where hyperedges are
 * created from each vertex by including the vertex and all it's graph 
 * neighbors.
 */
#define NROWS 25
#define NCOLS 25
static char *hg[NROWS]={
"12   6                   ",
"123   7                  ",
" 234   8                 ",
"  345   9                ",
"   45    0               ",
"1    67   1              ",
" 2   678   2             ",
"  3   789   3            ",
"   4   890   4           ",
"    5   90    5          ",
"     6    12   6         ",
"      7   123   7        ", 
"       8   234   8       ",
"        9   345   9      ",
"         0   45    0     ",
"          1    67   1    ",
"           2   678   2   ",
"            3   789   3  ",
"             4   890   4 ",
"              5   90    5",
"               6    12   ", 
"                7   123  ",
"                 8   234 ",
"                  9   345", 
"                   0   45" 
};
/************** sample hypergraphs *******************************************/

#define DIFFERENT_PROCS_DIFFERENT_EDGE_WEIGHTS 1
#define EDGE_WEIGHT_EQUALS_GID 0
#define VTX_WEIGHT_EQUALS_GID 0

static int nprocs, rank;

struct _division{
 int row0, row1, col0, col1;   /* my pins */
 int ew0, ew1;                 /* my edge weights */
 int format;    /* compressed rows or compressed columns */
 int num_obj;   /* my number of vertices */
 int obj_gid[NCOLS];  /* global ID of my vertices */
};

static void divide_interval(int from, int to, 
  int *l, int *r, int nprocs, int rank);
static void show_partitions(int nrows, int ncols, int *parts);

/*
** Set options specifying how the objects (vertices) and
** hypergraph pins are divided up initially.
*/

void *exSetHGDivisions(
  int initial_partitioning,
  int who_has_edge_weights,
  int edge_weights_overlap,
  int format)
{
  int row0, row1, col0, col1;   /* my pins */
  int ew0, ew1;                 /* my edge weights */
  int a, b, l, r, i;
  int *id;
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
  if ((format != ZOLTAN_COMPRESSED_EDGE) &&
      (format != ZOLTAN_COMPRESSED_VERTEX)){
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

  /*
   * Initial partitioning of vertices (objects to be partitioned) 
   * among processes
   */
  divide_interval(0, NCOLS-1, &l, &r, nprocs, rank);

  div->num_obj = r - l + 1;
  id = div->obj_gid;
  for (i=l; i<=r; i++){
    *id++ = i;
  }

  /*
   * Division of hypergraph pins
   */

  if ((nprocs < 3) && (initial_partitioning == BLOCKS)){
    initial_partitioning = COMPLETE_ROWS;
  }

  if (initial_partitioning == COMPLETE_ROWS){
    col0 = 0; col1 = NCOLS-1;
    divide_interval(0, NROWS-1, &row0, &row1, nprocs, rank);
  }
  else if (initial_partitioning == COMPLETE_COLS){
    row0 = 0; row1 = NROWS-1;
    divide_interval(0, NCOLS-1, &col0, &col1, nprocs, rank);
  }
  else{
    a = nprocs / 2;
    b = NROWS / 2;
    if (rank < a){
      row0 = 0; row1 = b-1;
      divide_interval(0, NCOLS-1, &col0, &col1, a, rank);
    }
    else if (rank < 2*a){
      row0 = b; row1 = NROWS-1;
      divide_interval(0, NCOLS-1, &col0, &col1, a, rank-a);
    }
    else{  /* odd number of processes, last one doesn't get a block */
      col0 = row0 = -1;
      col1 = row1 = -1;
    }
  }
  ew0 = ew1 = -1;
  if (who_has_edge_weights == ALL_HAVE_EDGE_WEIGHTS){
    /* we have edge weights, but not necc. for the same rows of our pins */
    divide_interval(0, NROWS-1, &ew0, &ew1, nprocs, nprocs-rank-1);
  }
  else if (who_has_edge_weights == ONE_HAS_EDGE_WEIGHTS){
    if (rank == 0){
      ew0 = 0;  ew1 = NROWS-1;
    }
  }
  else if (who_has_edge_weights == SOME_HAVE_EDGE_WEIGHTS){
    a = nprocs / 2;
    if (rank < a){
      divide_interval(0, NROWS-1, &ew0, &ew1, a, rank);
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

  return (void *)div;
}
/*
** After re-partitioning, update my list of objects.
*/
void exUpdateDivisions(void *data, int nexport, int nimport, 
                      ZOLTAN_ID_PTR exportGIDs, ZOLTAN_ID_PTR importGIDs)
{
int n_old_obj, n_new_obj = 0;
int i;
char gid_buf[NCOLS];
struct _division *div;

  div = (struct _division *)data;
  n_old_obj = div->num_obj;

  for (i=0; i<NCOLS; i++){
    gid_buf[i] = 0;
  }
  for (i=0; i<nimport; i++){
    gid_buf[importGIDs[i]] = 1;
  }
  for (i=0; i<n_old_obj; i++){
    gid_buf[div->obj_gid[i]] = 1;
  }
  for (i=0; i<nexport; i++){
    gid_buf[exportGIDs[i]] = 0;
  }
  for (i=0; i<NCOLS; i++){
    if (gid_buf[i]){
      div->obj_gid[n_new_obj++] = i;
    }
  }
  div->num_obj = n_new_obj;
}

/*
** Display the matrix representing the hypergraph, with partition numbers
*/
void exShowPartitions(void *data)
{
  struct _division *div = (struct _division *)data;
  int i, j, sendsize;
  int buf[NCOLS];
  int partitions[NCOLS];
  int tag=0x0101;
  MPI_Status status;

  if (rank == 0){
     
    for (j=0; j < NCOLS; j++){
      partitions[j] = -1;
    }

    for (i=0; i<div->num_obj; i++){
      partitions[div->obj_gid[i]] = 0;
    }

    for (i=1; i < nprocs; i++){
      MPI_Recv(buf, NCOLS, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

      for (j=0; j < NCOLS; j++){
        if (buf[j] < 0) break;
        partitions[buf[j]] = i;
      }
    }

    show_partitions(NROWS, NCOLS, partitions);
  }
  else{
    if (div->num_obj < NCOLS){
      div->obj_gid[div->num_obj] = -1;
      sendsize = div->num_obj+1;
    }
    else{
      sendsize = div->num_obj;
    }
    MPI_Send(div->obj_gid, sendsize, MPI_INT, 0, tag, MPI_COMM_WORLD);
  }
}

/*
** Release memory used by library
*/
void exFreeDivisions(void *data)
{
  free(data);
}

/*
** The query functions
*/
void exGetHgSizeAndFormat(void *data, 
   int *num_lists, int *num_pins, int *format, int *ierr)
{
struct _division *div;
int nedges=0;
int nverts=0;
int r, c;

  /*
   * Supply this query function to Zoltan with Zoltan_Set_HG_Size_CS_Fn().
   * It tells Zoltan the number of rows or columns to be supplied by
   * the process, the number of pins (non-zeroes) in the rows or columns, 
   * and whether the pins are provided in compressed row format or
   * compressed column format.
   */

  div = (struct _division *)data;

  *format = div->format;
  *num_pins = 0;

  if (div->row0 >= 0){
    nedges = div->row1 - div->row0 + 1;
    nverts = div->col1 - div->col0 + 1;

    for (r=div->row0; r<=div->row1; r++){
      for (c=div->col0; c<=div->col1; c++){
        if (!isspace(hg[r][c])){
          (*num_pins)++;
        }
      }
    }
  }

  if (div->format == ZOLTAN_COMPRESSED_EDGE){
    *num_lists = nedges;
  }
  else{
    *num_lists = nverts;
  }
  *ierr =ZOLTAN_OK; 
}

void exGetHg(void *data,  int num_gid_entries,
  int nrowcol, int npins, int format,
  ZOLTAN_ID_PTR rowcol_GID, int *rowcol_ptr, ZOLTAN_ID_PTR pin_GID, int *ierr)
{
  int i, j, e, v;
  struct _division *div;

  /*
   * Supply this query function to Zoltan with Zoltan_Set_HG_CS_Fn().
   * It supplies hypergraph pins to Zoltan.
   */

  div = (struct _division *)data;

  if (npins > 0){
    if (format == ZOLTAN_COMPRESSED_EDGE){
      for (i=0, e=div->row0, j=0; e <= div->row1; i++, e++){
        if (i == nrowcol) break;
        rowcol_ptr[i] = j;
        rowcol_GID[i] = e;

        for (v=div->col0; v <= div->col1; v++){
          if (!isspace(hg[e][v])){
            pin_GID[j++] = v;            
            if (j == npins) break;
          }
        }
      }
    }
    else{
      for (i=0, v=div->col0, j=0; v <= div->col1; i++, v++){
        if (i == nrowcol) break;
        rowcol_ptr[i] = j;
        rowcol_GID[i] = v;

        for (e=div->row0; e <= div->row1; e++){
          if (!isspace(hg[e][v])){
            pin_GID[j++] = e;            
            if (j == npins) break;
          }
        }
      }
    }
  }

  *ierr = ZOLTAN_OK;
}

void exGetHgEdgeWeightSize(void *data, int *num_edges, int *ierr)
{
  struct _division *div;

  /*
   * Supply this query function to Zoltan with Zoltan_Set_HG_Size_Edge_Weights_Fn().
   * It tells Zoltan the number edges for which this process will supply
   * edge weights.  The number of weights per edge was specified with the
   * parameter EDGE_WEIGHT_DIM.  If more than one process provides a weight
   * for the same edge, the multiple weights are resolved according the
   * value for the PHG_EDGE_WEIGHT_OPERATION parameter.
   */

  div = (struct _division *)data;

  if (div->ew0 >= 0){
    *num_edges = div->ew1 - div->ew0 + 1;
  }
  else{
    *num_edges = 0;
  }
  *ierr = ZOLTAN_OK;
}

void exGetHgEdgeWeights(void *data,  int num_gid_entries,
  int num_lid_entries, int nedges, int edge_weight_dim,
  ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float *edge_weight, int *ierr)
{
  struct _division *div;
  int i, e;

  /*
   * Supply this query function to Zoltan with Zoltan_Set_HG_Edge_Weights_Fn().
   * It tells Zoltan the weights for some subset of edges.
   */

  div = (struct _division *)data;

  if (div->ew0 >= 0){

    for (i=0, e=div->ew0; e <= div->ew1; i++, e++){
      edge_GID[i] = e;
      edge_LID[i] = i;
#if DIFFERENT_PROCS_DIFFERENT_EDGE_WEIGHTS
      edge_weight[i] = rank+1;
#else
  #if EDGE_WEIGHT_EQUALS_GID
      edge_weight[i] = e+1;
  #else
      edge_weight[i] = 1.0;
  #endif
#endif
    }
  }

  *ierr = ZOLTAN_OK;
}
int exGetHgNumVertices(void *data, int *ierr)
{
  struct _division *div;

  /*   
   * Supply this query function to Zoltan with Zoltan_Set_Num_Obj_Fn().
   * It returns the number of vertices that this process owns.
   *
   * The parallel hypergraph method requires that vertex global IDs and 
   * weights are returned by the application with query functions.  The
   * global IDs should be unique, and no two processes should
   * return the same vertex.
   */
  div = (struct _division *)data;

  return div->num_obj;
}
void exGetHgVerticesAndWeights(void *data, int num_gid_entries,
  int num_lid_entries, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
  int wgt_dim, float *obj_weights, int *ierr)
{
  int i, nv, j;
  float *w;
  struct _division *div;

  div = (struct _division *)data;

  /*
   * Supply this query function to Zoltan with Zoltan_Set_Obj_List_Fn().
   *
   * It supplies vertex global IDs, local IDs,
   * and weights to the Zoltan library.  The application has previously
   * indicated the number of weights per vertex by setting the
   * OBJ_WEIGHT_DIM parameter.
   */

  nv = div->num_obj;
  w = obj_weights;

  if (nv > NCOLS){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0; i<nv; i++){
    gids[i] = div->obj_gid[i];
    lids[i] = i;

    for (j=0; j<wgt_dim; j++){
#if VTX_WEIGHT_EQUALS_GID
      *w++ = (j ? 1.0 : gids[i]);
#else
      *w++ = 1.0;
#endif
    }
  }

  *ierr = ZOLTAN_OK;
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

static void show_partitions(int nrows, int ncols, int *parts)
{
  char p[128];
  int r, c, part, ncuts, maxpart, nspaces;
 
  maxpart = -1;
  for (c=0; c<ncols; c++){
    if (parts[c] > maxpart) maxpart = parts[c];
  }
 
  for (r=0; r<nrows; r++){
     memset(p, 0, 128);
     ncuts = -1;
     for (c=0; c<ncols; c++){
       if (!isspace(hg[r][c])){
         part = parts[c];
         if (maxpart > 9) printf("%02d", part);
         else             printf("%d", part);
         if (p[part] == 0){
           ncuts++;
           p[part] = 1;
         }
       }
       else{
         printf(" ");
         if (maxpart > 9) printf(" ");
       }
     }
     printf(" (%d cuts)\n",ncuts);
  }
}
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
