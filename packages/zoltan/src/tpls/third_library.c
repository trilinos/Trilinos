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


#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "order_const.h"
#include "third_library.h"

/**********  parameters structure used by PHG, ParMetis and Jostle **********/
static PARAM_VARS Graph_Package_params[] = {
        { "GRAPH_PACKAGE", NULL, "STRING", 0 },
        { "ORDER_TYPE", NULL, "STRING", 0 },
        { NULL, NULL, NULL, 0 } };


/**********************************************************/
/* Interface routine for Graph methods.                   */
/**********************************************************/

int Zoltan_Graph(
  ZZ *zz,               /* Zoltan structure */
  float *part_sizes,    /* Input:  Array of size zz->Num_Global_Parts
                           containing the percentage of work to be
                           assigned to each partition.               */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int **imp_to_part,    /* list of partitions to which imported objects are
                           assigned.  */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,      /* list of processors to export to */
  int **exp_to_part     /* list of partitions to which exported objects are
                           assigned. */
)
{
  static char* yo = "Zoltan_Graph";

char *defaultMethod= "PHG";
char package[MAX_PARAM_STRING_LEN];
int rc;

  strcpy(package, defaultMethod);
  Zoltan_Bind_Param(Graph_Package_params, "GRAPH_PACKAGE", package);
  Zoltan_Assign_Param_Vals(zz->Params, Graph_Package_params, zz->Debug_Level,
          zz->Proc, zz->Debug_Proc);

  if (!strcasecmp(package, "PARMETIS")){
#ifdef ZOLTAN_PARMETIS
    rc = Zoltan_ParMetis(zz, part_sizes, num_imp, imp_gids, imp_lids,
                         imp_procs, imp_to_part,
                         num_exp, exp_gids, exp_lids, exp_procs, exp_to_part);
#else
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "ParMETIS partitioning was requested but "
                       "Zoltan was compiled without ParMETIS.\n");
    rc = ZOLTAN_FATAL;
#endif /* ZOLTAN_PARMETIS */
  }
  else if (!strcasecmp(package, "SCOTCH")){
#ifdef ZOLTAN_SCOTCH
    rc = Zoltan_Scotch(zz, part_sizes, num_imp, imp_gids, imp_lids,
                         imp_procs, imp_to_part,
                         num_exp, exp_gids, exp_lids, exp_procs, exp_to_part);
#else
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Scotch partitioning was requested but "
                       "Zoltan was compiled without Scotch.\n");
    rc = ZOLTAN_FATAL;
#endif /* ZOLTAN_SCOTCH */
  }
  else if (!strcasecmp(package, "ZOLTAN") ||
           !strcasecmp(package, "PHG")) {

    rc = Zoltan_PHG(zz, part_sizes, num_imp, imp_gids, imp_lids,
                         imp_procs, imp_to_part,
                         num_exp, exp_gids, exp_lids, exp_procs, exp_to_part);
  }
  else{
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Invalid value for GRAPH_PACKAGE parameter\n");
    rc = ZOLTAN_FATAL;
  }

  return rc;
}


/*********************************************************************/
/* Graph_Package parameter routine                                        */
/*********************************************************************/

int Zoltan_Graph_Package_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int status, i;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */
  char *valid_methods[] = {
    "PARMETIS", "PHG", "ZOLTAN", "SCOTCH",
    NULL };

  status = Zoltan_Check_Param(name, val, Graph_Package_params,
			      &result, &index);

  if (status == 0){
    /* OK so far, do sanity check of parameter values */

    if (strcmp(name, "GRAPH_PACKAGE") == 0){
      status = 2;
      for (i=0; valid_methods[i] != NULL; i++){
	if (strcmp(val, valid_methods[i]) == 0){
	  status = 0;
#ifndef ZOLTAN_PARMETIS
          if (strcmp(val, "PARMETIS") == 0){
            status = 2;
          }
#endif
#ifndef ZOLTAN_SCOTCH
          if (strcmp(val, "SCOTCH") == 0){
            status = 2;
          }
#endif
	  break;
	}
      }
    }
  }

  return(status);
}


#define MEMFREE(ptr) do { if (ptr) ZOLTAN_FREE(&(ptr)); } while (0);


void Zoltan_Third_Exit(ZOLTAN_Third_Graph *gr, ZOLTAN_Third_Geom *geo,
		       ZOLTAN_Third_Part *prt, ZOLTAN_Third_Vsize *vsp,
		       ZOLTAN_Output_Part *part, ZOLTAN_Output_Order *ord)
{
  if (gr) {

    Zoltan_Matrix2d_Free(&gr->graph.mtx);

    MEMFREE(gr->vwgt);
    MEMFREE(gr->vtxdist);
    MEMFREE(gr->xadj);
    MEMFREE(gr->adjncy);
    MEMFREE(gr->ewgts);
    MEMFREE(gr->float_ewgts);
    MEMFREE(gr->adjproc);
    
    Zoltan_ZG_Free(&gr->graph);
  }

  if (geo) {
    MEMFREE(geo->xyz);
  }

  if (prt) {
    MEMFREE(prt->part);
    MEMFREE(prt->input_part);
    MEMFREE(prt->part_orig);
    if (prt->part_sizes != prt->input_part_sizes)
      MEMFREE(prt->part_sizes);
    if (sizeof(realtype) != sizeof(float))
      MEMFREE(prt->input_part_sizes);
  }

  if (vsp) {
    if (!vsp->vsize_malloc) {
      MEMFREE(vsp->vsize);
    }
    MEMFREE(vsp->vsizeBACKUP);
  }

  if (ord) {
    MEMFREE(ord->sep_sizes);
    MEMFREE(ord->rank);
    MEMFREE(ord->iperm);
  }
}

int Zoltan_Third_Init(ZOLTAN_Third_Graph *gr, ZOLTAN_Third_Part  *prt, ZOLTAN_Third_Vsize *vsp, ZOLTAN_Output_Part *part,
		      ZOLTAN_ID_PTR *imp_gids, ZOLTAN_ID_PTR *imp_lids, int **imp_procs, int **imp_to_part,
		      ZOLTAN_ID_PTR *exp_gids, ZOLTAN_ID_PTR *exp_lids, int **exp_procs, int **exp_to_part)
{

  memset (gr, 0, sizeof(ZOLTAN_Third_Graph));
  memset (prt, 0, sizeof(ZOLTAN_Third_Part));
  memset (vsp, 0, sizeof(ZOLTAN_Third_Vsize));
  memset (part, 0, sizeof(ZOLTAN_Output_Part));

  /* Initialize return-argument arrays to return arguments so that F90 works. */
  part->imp_gids = imp_gids;
  part->imp_lids = imp_lids;
  part->imp_procs = imp_procs;
  part->imp_part = imp_to_part;

  part->exp_gids = exp_gids;
  part->exp_lids = exp_lids;
  part->exp_procs = exp_procs;
  part->exp_part = exp_to_part;

  part->num_imp = part->num_exp = -1;

  /* Most ParMetis methods use only graph data */
  gr->get_data = 1;

  return (ZOLTAN_OK);
}

/* export to user variables */
int Zoltan_Third_Export_User(ZOLTAN_Output_Part *part,
			     int *num_imp, ZOLTAN_ID_PTR *imp_gids, ZOLTAN_ID_PTR *imp_lids, int **imp_procs, int **imp_to_part,
			     int *num_exp, ZOLTAN_ID_PTR *exp_gids, ZOLTAN_ID_PTR *exp_lids, int **exp_procs, int **exp_to_part)
{
  /* Write results in user variables */
  *num_imp = part->num_imp;
  *imp_gids = *(part->imp_gids);
  *imp_lids = *(part->imp_lids);
  *imp_procs = *(part->imp_procs);
  *imp_to_part = *(part->imp_part);
  *num_exp = part->num_exp;
  *exp_gids = *(part->exp_gids);
  *exp_lids = *(part->exp_lids);
  *exp_procs = *(part->exp_procs);
  *exp_to_part = *(part->exp_part);

  return (ZOLTAN_OK);
}

int Zoltan_matrix_Print(Zoltan_matrix *m, char *s)
{
int i, j, k;
float *wgts;

  if (s) fprintf(stderr,"Zoltan_matrix, %s\n",s);
  fprintf(stderr,"\nOptions: enforceSquare %d, pinwgtop %s, randomize %d, pinwgt %d\n",
     m->opts.enforceSquare,
     ((m->opts.pinwgtop == 0) ? "add weight" : ((m->opts.pinwgtop == 1) ? "max weight" : "cmp weight")),
     m->opts.randomize, m->opts.pinwgt);

  fprintf(stderr,"Options: local %d, final_output %d, symmetrize %d keep_distribution %d speed %s\n",
    m->opts.local, m->opts.final_output, m->opts.symmetrize, m->opts.keep_distribution,
     ((m->opts.speed == 0) ? "full dd" : ((m->opts.speed == 1) ? "fast" : "no redist")));
   
  fprintf(stderr,"redist %d, completed %d, bipartite %d\n", m->redist, m->completed, m->bipartite);

  fprintf(stderr,"globalX " ZOLTAN_GNO_SPEC ", globalY " ZOLTAN_GNO_SPEC ", nY %d, nY_ori %d, ywgtdim %d, nPins %d\n",
                m->globalX, m->globalY, m->nY, m->nY_ori, m->ywgtdim, m->nPins);
  fprintf(stderr,"Edges and non-zeroes:\n");
  wgts = m->pinwgt;
  if (m->yGNO && m->pinGNO){
    for (i=0; i < m->nY; i++){
      fprintf(stderr, ZOLTAN_GNO_SPEC ": ",m->yGNO[i]);
      for (j=m->ystart[i]; j < m->ystart[i+1]; j++){
        fprintf(stderr, ZOLTAN_GNO_SPEC " ", m->pinGNO[j]);
        if (wgts && (m->pinwgtdim > 0)){
          fprintf(stderr,"("); 
          for (k=0; k < m->pinwgtdim; k++){
            fprintf(stderr,"%f ",*wgts++);
          }
          fprintf(stderr,") "); 
        }
      }
      fprintf(stderr,"\n");
    }
  }
  else{
    fprintf(stderr,"not set");
  }
  fprintf(stderr,"\n");
  fflush(stderr);
  return ZOLTAN_OK;
}

int Zoltan_ZG_Print(ZZ *zz, ZG *gr, char *s)
{
int i, me, proc;
Zoltan_matrix_2d *m2d = &gr->mtx;
Zoltan_matrix *m = &m2d->mtx;
int nproc_x = m2d->comm->nProc_x;
int nproc_y = m2d->comm->nProc_y;

  me = zz->Proc;
  for (proc=0; proc < zz->Num_Proc; proc++){
    if (proc == me){
      if (proc == 0) fprintf(stderr,"\n%s\n",s);
      fprintf(stderr,"Process: %d) flags: bipartite %d fixObj %d, fixed vertices buffer %p:\n",
               zz->Proc, gr->bipartite, gr->fixObj, (void *)gr->fixed_vertices);
      fprintf(stderr,"GNO distribution in x direction: ");
      if (m2d->dist_x){
        for (i=0; i <= nproc_x; i++){
          fprintf(stderr, ZOLTAN_GNO_SPEC " ",m2d->dist_x[i]);
        }
      }
      else{
        fprintf(stderr,"not set");
      }
      fprintf(stderr,"\nGNO distribution in y direction: ");
      if (m2d->dist_y){
        for (i=0; i <= nproc_y; i++){
          fprintf(stderr, ZOLTAN_GNO_SPEC " ",m2d->dist_y[i]);
        }
      }
      else{
        fprintf(stderr,"not set");
      }

      Zoltan_matrix_Print(m, NULL);
  
      fflush(stderr);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return ZOLTAN_OK;
}
int Zoltan_Third_Graph_Print(ZZ *zz, ZOLTAN_Third_Graph *gr, char *s)
{
int i, me, proc;
indextype numvtx, offset, j;
me = zz->Proc;

  for (proc=0; proc < zz->Num_Proc; proc++){
    if (proc == me){
      if (proc == 0) fprintf(stderr,"\n%s\n",s);
      fprintf(stderr,"Process: %d) graph type %d, check graph %d, final output %d, showMoveVol %d, scatter %d\n",
      me, gr->graph_type, gr->check_graph, gr->final_output, gr->showMoveVol, gr->scatter);
      fprintf(stderr,"scatter min %d, get data %d, obj wgt dim %d, edge wgt dim %d\n",
      gr->scatter_min, gr->get_data, gr->obj_wgt_dim, gr->edge_wgt_dim);
      fprintf(stderr,"num obj %d, num obj orig %d, num edges %d\n",
      gr->num_obj, gr->num_obj_orig, gr->num_edges);
  
      if (gr->vtxdist){
        numvtx = gr->vtxdist[proc+1] - gr->vtxdist[proc];
        offset = gr->vtxdist[proc];
        fprintf(stderr,"Num vertices: " TPL_IDX_SPEC "\n",numvtx);
        if (gr->xadj){
          fprintf(stderr,"Num edges: " TPL_IDX_SPEC "\n",gr->xadj[numvtx]);
        }

        for (i=0; i < numvtx; i++){
          fprintf(stderr,TPL_IDX_SPEC ": ",i+offset);
          if (gr->xadj){
            for (j=gr->xadj[i];j < gr->xadj[i+1]; j++){
              if (gr->adjncy){
                fprintf(stderr,"gid " TPL_IDX_SPEC,gr->adjncy[j]);
              }
              if (gr->adjproc){
                fprintf(stderr," proc %d ",gr->adjproc[j]);
              }
            }
          }
          else{
            fprintf(stderr,"adjacency info is null");
          }
          fprintf(stderr,"\n");
        }
      }
  
      fflush(stderr);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return ZOLTAN_OK;
}

/*
** Copy of Zoltan_Order_Init_Tree in order_tools.c which uses TPL data types
*/
#define CHECK_ZOLTAN_FREE(ptr) do { if ((ptr) != NULL) ZOLTAN_FREE(&(ptr)); } while (0)

int Zoltan_TPL_Order_Init_Tree (struct Zoltan_TPL_Order_Struct *order, indextype blocknbr, indextype leavesnbr)
{   
  Zoltan_TPL_Order_Free_Struct(order);

  order->ancestor = (indextype *) ZOLTAN_MALLOC(blocknbr*sizeof(indextype));
  order->start = (indextype *) ZOLTAN_MALLOC((blocknbr+1)*sizeof(indextype));
  order->leaves = (indextype *) ZOLTAN_MALLOC((leavesnbr+1)*sizeof(indextype));
  
  if ((order->ancestor == NULL) || (order->start == NULL) || (order->leaves == NULL)) {
    Zoltan_TPL_Order_Free_Struct(order);
    return (ZOLTAN_MEMERR);
  }
  order->needfree = 1;
  return (ZOLTAN_OK);
} 

void Zoltan_TPL_Order_Free_Struct(struct Zoltan_TPL_Order_Struct *order)
{
  if (order->needfree == 0)
    return;

  CHECK_ZOLTAN_FREE(order->start);
  CHECK_ZOLTAN_FREE(order->ancestor);
  CHECK_ZOLTAN_FREE(order->leaves);

  order->needfree = 0;
}

    


#ifdef __cplusplus
}
#endif
