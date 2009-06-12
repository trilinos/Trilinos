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
    MEMFREE(gr->vtxdist);
    MEMFREE(gr->xadj);
    MEMFREE(gr->adjncy);
    MEMFREE(gr->vwgt);
    MEMFREE(gr->ewgts);
    MEMFREE(gr->float_ewgts);
    MEMFREE(gr->adjproc);
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


#ifdef __cplusplus
}
#endif
