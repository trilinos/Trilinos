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
#include "matrix.h"
#include "third_library.h"
#include "ccolamd_interface.h"

  /**********  parameters structure for parmetis methods **********/
static PARAM_VARS CColAMD_params[] = {
  { NULL, NULL, NULL, 0 } };


/***************************************************************************
 * External function to compute CCOLAMD ordering,
 * used by Zoltan_order.
 **************************************************************************/

int Zoltan_CColAMD_Order(
  ZZ *zz,               /* Zoltan structure */
  int num_obj,		/* Number of (local) objects to order. */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
  /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
/* The application must allocate enough space */
  int *rank,		/* rank[i] is the rank of gids[i] */
  int *iperm,
  ZOOS *order_opt 	/* Ordering options, parsed by Zoltan_Order */
)
{
  static char *yo = "Zoltan_CColAMD_Order";
  int n, ierr;

  return (ZOLTAN_OK);
}


/***************************************************************************
 * External function to compute CCOLAMD ordering,
 * used by Zoltan_order.
 * Compute by /slices/ of the matrix, according to the partition given
 * in parameter.
 **************************************************************************/

int Zoltan_CColAMD(
  ZZ *zz,               /* Zoltan structure */
  struct Zoltan_DD_Struct *dd_constraint,
  int *num_obj,
  ZOLTAN_ID_PTR *gids,
  int **rank
)
{
  static char *yo = "Zoltan_CColAMD";
  int n, ierr = ZOLTAN_OK;
  Zoltan_matrix_options opt;
  double knobs [CCOLAMD_KNOBS];
  int stats [CCOLAMD_STATS];
  size_t Alen;
  int *pins = NULL;         /* Ccolamd needs a copy of the non-zeros */
  int *cmember = NULL;      /* constraints */
  void *partdata = NULL;
  int *ystart = NULL;
  int n_col, n_row, n_nnz;
  Zoltan_matrix_2d mtx;
  int i;
  int offset = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);
  memset (&opt, 0, sizeof(Zoltan_matrix_options));
  opt.speed = MATRIX_NO_REDIST;

  mtx.comm = (PHGComm*)ZOLTAN_MALLOC (sizeof(PHGComm));
  if (mtx.comm == NULL) MEMORY_ERROR;
  Zoltan_PHGComm_Init (mtx.comm);

  /* Construct a CSC matrix */
  /* TODO: take this in parameter instead */
  ierr = Zoltan_Matrix_Build(zz, &opt, &mtx.mtx);
  CHECK_IERR;

  /* Set up the correct distribution function */
  /* TODO: do here to avoid a first distribution */

  ierr = Zoltan_Distribute_LinearY(zz, mtx.comm);
  CHECK_IERR;
  ierr = Zoltan_Matrix2d_Distribute (zz, mtx.mtx, &mtx, 0);
  CHECK_IERR;
  ierr = Zoltan_Matrix_Complete(zz, &mtx.mtx);
  CHECK_IERR;

  n_col = mtx.mtx.nY;

  ierr = Zoltan_DD_Find (dd_constraint, mtx.mtx.yGID, (ZOLTAN_ID_PTR)cmember, NULL, NULL,
			 mtx.mtx.nY, NULL);
  CHECK_IERR;
  partdata = Zoltan_Distribute_Partition_Register(zz, n_col, mtx.mtx.yGNO, cmember);
  ZOLTAN_FREE(&cmember);
  Zoltan_Distribute_Set(&mtx, Zoltan_Distribute_Partition, partdata);
  ierr = Zoltan_Matrix2d_Distribute (zz, mtx.mtx, &mtx, 0);
  CHECK_IERR;
  ierr = Zoltan_Matrix_Complete(zz, &mtx.mtx);
  CHECK_IERR;
  Zoltan_Distribute_Partition_Free(&partdata);


  (*num_obj) = n_col = mtx.mtx.nY;
  n_row = mtx.mtx.globalX;
  n_nnz = mtx.mtx.nPins;


  /* Prepare call to CCOLAMD */
  ccolamd_set_defaults (knobs);

  cmember = (int*) ZOLTAN_MALLOC(n_col * sizeof(int));
  if (n_col > 0 && cmember == NULL) MEMORY_ERROR;

  (*gids) = ZOLTAN_MALLOC_GID_ARRAY(zz , n_col);
  if (n_col > 0 && (*gids) == NULL) MEMORY_ERROR;
  memcpy ((*gids), mtx.mtx.yGID, n_col*sizeof(int)*zz->Num_GID);

  Alen = ccolamd_recommended (n_nnz, n_row, n_col);
  pins = (int*) ZOLTAN_MALLOC(Alen * sizeof(int));
  if (Alen >0 && pins == NULL) MEMORY_ERROR;
  memcpy (pins, mtx.mtx.pinGNO, mtx.mtx.nPins*sizeof(int));

  ystart = (int*) ZOLTAN_MALLOC((n_col + 1) * sizeof(int));
  if (ystart == NULL) MEMORY_ERROR;
  memcpy (ystart, mtx.mtx.ystart, (n_col + 1) * sizeof(int));

  Zoltan_Matrix2d_Free(&mtx);
  ZOLTAN_FREE(&mtx.comm);

  /* Compute ordering */
  /* Upon return, ystart is the invert permutation ... */
  ierr = ccolamd (n_row, n_col, Alen, pins, ystart,
	   knobs, stats, cmember);

  ZOLTAN_FREE(&pins);
  ZOLTAN_FREE(&cmember);

  (*rank) = (int*) ZOLTAN_MALLOC(n_col * sizeof(int));
  if (n_col > 0 && (*rank) == NULL) MEMORY_ERROR;


  /* Compute offset in the global graph */
  MPI_Scan(&n_col, &offset, 1, MPI_INT, MPI_SUM, zz->Communicator);
  offset -= n_col;
  /* Compute direct permutation */
  for (i = 0 ; i < n_col ; ++i) {
    (*rank)[ystart[i]] = i + offset;
  }

/*   memcpy ((*rank), ystart, n_col * sizeof(int)); */

 End:
  ZOLTAN_FREE(&pins);
  ZOLTAN_FREE(&ystart);
  ZOLTAN_FREE(&cmember);


  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}



/*********************************************************************/
/* ParMetis parameter routine                                        */
/*********************************************************************/

int Zoltan_CColAMD_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status, i;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    char *valid_methods[] = {
         NULL };

    status = Zoltan_Check_Param(name, val, CColAMD_params, &result, &index);

    if (status == 0){
      /* OK so far, do sanity check of parameter values */

      if (strcmp(name, "PARMETIS_METHOD") == 0){
        status = 2;
        for (i=0; valid_methods[i] != NULL; i++){
          if (strcmp(val, valid_methods[i]) == 0){
            status = 0;
            break;
          }
        }
      }
    }
    return(status);
}

#ifdef __cplusplus
}
#endif

