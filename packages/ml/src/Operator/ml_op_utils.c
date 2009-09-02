/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_struct.h"
#include "ml_op_utils.h"
#include "ml_agg_genP.h"
#include "ml_memory.h"
#include "ml_agg_Zoltan.h"
#if defined(HAVE_ML_ZOLTAN) && defined(HAVE_MPI)
#include "zoltan.h"
#include "ml_viz_stats.h"
#endif

/* ******************************************************************** */
/* Blow away any inter-mixing between boundary and interior points in   */
/* the matrix ml_handle->Pmat[level2].                                  */
/* -------------------------------------------------------------------- */

int oldML_Mdfy_Prolongator_DirBdry(ML *ml_handle, int level2, 
                                   double *boundary, double *cboundary)
{
   int *cols, *row_ptr, i, j;
   double *vals;
   struct ML_CSR_MSRdata *temp;

   int size;

   if (ml_handle->Pmat[level2].getrow->func_ptr != CSR_getrow)
     perror("ML_Mdfy_Prolongator_DirBdry can only be used with CSR matrices\n");

   temp    = (struct ML_CSR_MSRdata *) ml_handle->Pmat[level2].data;
   cols    = temp->columns;
   vals    = temp->values;
   row_ptr = temp->rowptr;
   size    = ml_handle->Pmat[level2].outvec_leng;

   /* zero out any inter-mixing between boundary  */
   /* and interior in the interpolation operator. */

   for (i = 0; i < size; i++) {
      if (boundary[i] == 1.) {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( cboundary[cols[j]] == 0.) vals[j] = 0.0;
      }
      else {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( cboundary[cols[j]] == 1.) vals[j] = 0.0;
      }
   }
   return(1);
}

/* ******************************************************************** */
/* Blow away any inter-mixing between boundary and interior points in   */
/* the matrix ml_handle->Pmat[level2].                                  */
/* -------------------------------------------------------------------- */

int ML_Mdfy_Prolongator_DirBdry(ML *ml_handle, int level2, int size,
     int fine_size )
{
   int fboundary_length,cboundary_length, *fboundary_list, *cboundary_list;
   char *f_bdry;
   int *cols, *row_ptr, i, j, jk;
   double *vals;
   struct ML_CSR_MSRdata *temp;
   double *dtemp;
   ML_CommInfoOP *comm_info;

   comm_info        = ml_handle->Pmat[level2].getrow->pre_comm;
   fboundary_length = ml_handle->Pmat[level2].to->BCs->Dirichlet_grid_length;
   fboundary_list   = ml_handle->Pmat[level2].to->BCs->Dirichlet_grid_list;
   cboundary_length = ml_handle->Pmat[level2].from->BCs->Dirichlet_grid_length;
   cboundary_list   = ml_handle->Pmat[level2].from->BCs->Dirichlet_grid_list;
   dtemp       = (double *) ML_allocate((size+1)*sizeof(double));
   f_bdry      = (char *) ML_allocate((fine_size+1)*sizeof(char));
   if (f_bdry == NULL) {
        printf("No space to compute coarse boundary\n");
        exit(1);
   }
   for (jk = 0; jk < fine_size; jk++)       f_bdry[jk] = 'i';
   for (jk = 0; jk < fboundary_length; jk++) f_bdry[fboundary_list[jk]] = 'b';
   for (jk = 0; jk < size; jk++)       dtemp[jk] = 0.0;
   for (jk = 0; jk < cboundary_length; jk++) dtemp[cboundary_list[jk]] = 1.;
   if ( comm_info != NULL)
      ML_exchange_bdry(dtemp, comm_info, size, ml_handle->comm,
                       ML_OVERWRITE,NULL);



   if (ml_handle->Pmat[level2].getrow->func_ptr != CSR_getrow)
     perror("ML_Mdfy_Prolongator_DirBdry can only be used with CSR matrices\n");

   temp    = (struct ML_CSR_MSRdata *) ml_handle->Pmat[level2].data;
   cols    = temp->columns;
   vals    = temp->values;
   row_ptr = temp->rowptr;

   /* zero out any inter-mixing between boundary  */
   /* and interior in the interpolation operator. */

   for (i = 0; i < fine_size; i++) {
      if (f_bdry[i] == 'b') {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( dtemp[cols[j]] == 0.0) vals[j] = 0.0;
      }
      else {
         for (j = row_ptr[i] ; j < row_ptr[i+1]; j++)
            if ( dtemp[cols[j]] == 1.0) vals[j] = 0.0;
      }
   }
   ML_free(dtemp);
   ML_free(f_bdry);

   return(1);
}

/* ******************************************************************** */
/* -------------------------------------------------------------------- */

int ML_Compute_Coarse_Bdry(ML *ml_handle, int level, int size, int fine_size)
{
   int     *cols, *row_ptr, jk, Ncoarse;
   struct  ML_CSR_MSRdata *temp;
   int     *boundary_list, boundary_length, *cboundary_list, cboundary_length = 0;
   char    *f_bdry, *c_bdry;

   Ncoarse     = ml_handle->Pmat[level].invec_leng;
   c_bdry      = (char *) ML_allocate((size+1)*sizeof(char));
   f_bdry      = (char *) ML_allocate((fine_size+1)*sizeof(char));
   if (f_bdry == NULL) {
      printf("No space to compute coarse boundary\n");
      exit(1);
   }
   boundary_length = ml_handle->Pmat[level].to->BCs->Dirichlet_grid_length;
   boundary_list   = ml_handle->Pmat[level].to->BCs->Dirichlet_grid_list;
   for (jk = 0; jk < fine_size; jk++) f_bdry[jk] = 'i';
   for (jk = 0; jk < boundary_length; jk++) f_bdry[boundary_list[jk]] = 'b';

   /* Mark the coarse grid boundary */

   temp        = (struct ML_CSR_MSRdata*) ml_handle->Pmat[level].data;
   cols        = temp->columns;
   row_ptr     = temp->rowptr;

   for (jk = 0; jk < size; jk++) c_bdry[jk] = 'i';
   for (jk = 0; jk < fine_size; jk++) {
      if ((row_ptr[jk+1] - row_ptr[jk] == 1) && ( f_bdry[jk] == 'b')) {
         c_bdry[ cols[row_ptr[jk]]] = 'b';
      }
   }

   /* stuff the coarse boundary information into ML */

   for (jk = 0; jk < Ncoarse; jk++) {
      if (c_bdry[jk] == 'b') cboundary_length++;
   }
   cboundary_list = (int *) ML_allocate((cboundary_length+1)*sizeof(int));
   if (cboundary_list == NULL) {
      printf("No space to compute coarse boundary\n");
      exit(1);
   }
   cboundary_length = 0;
   for (jk = 0; jk < Ncoarse; jk++) {
      if (c_bdry[jk] == 'b') cboundary_list[cboundary_length++] = jk;
   }
   ML_Set_BoundaryTypes(ml_handle, ml_handle->Pmat[level].from->levelnum,
                        ML_BDRY_DIRICHLET,cboundary_length,cboundary_list);

   ML_free(c_bdry);
   ML_free(f_bdry);
   ML_free(cboundary_list);
   return(1);
}

/* ******************************************************************** */
/*
 * Take the 'matrix' defined by getrow() and vec_comm() and make a CSR
 * copy of it which will be stored in the matrix ml_handle->Pmat[level2].
 *
 *
 * Parameters
 * ==========
 *
 *    isize      On input, number of local columns in matrix to be copied.
 *
 *    osize      On input, number of local rows in matrix to be copied.
 *
 *    getrow()   On input, user's function to get rows of the matrix.
 *
 *    vec_comm() On input, user's function to update a vector via communication
 *               so that a local matrix-vector product can occur.
 *
 * -------------------------------------------------------------------- */

int ML_Gen_Prolongator_Getrow(ML *ml_handle, int level2, int level, int isize,
   int osize, int (*getrow)(void* , int , int *, int , int *, double *, int *),
   int (*vec_comm)(double *, void*), void *data, int Nghost)
{
   int *cols, *row_ptr, space, flag, nz_ptr, i, length;
   double *vals;
   double dsize, di;
   struct ML_CSR_MSRdata *temp;

   row_ptr = (int *) ML_allocate(sizeof(int)*(osize+1));

   space = osize*5+30;

   flag = 0;

   while (flag == 0) {
      cols    = (int    *) ML_allocate(sizeof(int)*space);
      vals    = (double *) ML_allocate(sizeof(double)*space);

      nz_ptr = 0;
      row_ptr[0] = nz_ptr;
      for (i = 0; i < osize; i++) {
         flag = getrow(data, 1, &i, space-nz_ptr, &(cols[nz_ptr]),
                       &(vals[nz_ptr]), &length);
         if (flag == 0) break;
         nz_ptr += length;
         row_ptr[i+1] = nz_ptr;
      }
      if (flag == 0) {
         dsize = (double) osize;
         di    = (double) (i+1);
         dsize = 1.2*dsize/di;
         space = (int) ( ((double) space)*dsize);
         space++;
         ML_free(vals);
         ML_free(cols);
      }
   }

   /* store the matrix into ML */

   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
   temp->columns = cols;
   temp->values  = vals;
   temp->rowptr = row_ptr;

   ml_handle->Pmat[level2].data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Init_Prolongator(ml_handle, level2, level, isize, osize, (void *) temp);
   ML_Operator_Set_ApplyFunc(&(ml_handle->Pmat[level2]),CSR_matvec);
/*
printf("we've changed the data pointer ? ....\n");
*/
   if (vec_comm != NULL) {
      ML_CommInfoOP_Generate(&(ml_handle->Pmat[level2].getrow->pre_comm), 
			  vec_comm, data, ml_handle->comm, 
			  ml_handle->Pmat[level2].invec_leng, Nghost);
   }
   else ml_handle->Pmat[level2].getrow->pre_comm = NULL;

   ML_Operator_Set_Getrow(&(ml_handle->Pmat[level2]), 
			  ml_handle->Pmat[level2].outvec_leng, CSR_getrow);

/*
   ML_CommInfoOP_Generate(&(ml_handle->Pmat[level2].getrow->pre_comm),
                   vec_comm, data, ml_handle->comm, 
		   ml_handle->Pmat[level2].invec_leng, Nghost);
*/

   return(1);
}

/* ******************************************************************** */
/* Take the Prolongator given by ml_handle->Pmat[level2], transpose
 * it and store it into the matrix given by ml_handle->Rmat[level].
 *
 * Added on Jan-06, MS.
 * If the operator to transpose is available outside the ml_handle
 * thing, then one can pass it using Ptent; otherwise set Ptent to NULL
 * to transpose Pmat = &(ml_handle->Pmat[level2]).
 * -------------------------------------------------------------------- */

int ML_Gen_Restrictor_TransP(ML *ml_handle, int level, int level2,
                             ML_Operator* Ptent)
{

   ML_Operator *Pmat, *Rmat;
   int  *row_ptr, *colbuf, *cols;
   int isize, osize, i, j, N_nzs, flag, length, sum, new_sum;
   int Nneighbors, *neigh_list, *send_list, *rcv_list, Nsend, Nrcv;
   void *data = NULL;
   double *valbuf, *vals;
   int (*getrow)(ML_Operator* , int , int *, int , int *, double *, int *) = NULL;
   struct ML_CSR_MSRdata *temp;
   int Nghost = 0, Nghost2 = 0;
   int *remap, remap_leng;
   ML_CommInfoOP *c_info, **c2_info;


   /* pull out things from ml_handle */

   if (Ptent == NULL) /* this is the default behavior */
     Pmat  = &(ml_handle->Pmat[level2]);
   else
     Pmat  = Ptent;

   temp = (struct ML_CSR_MSRdata *) Pmat->data;
   Rmat  = &(ml_handle->Rmat[level]);
   isize = Pmat->outvec_leng;
   osize = Pmat->invec_leng;
   data   = (void *) Pmat;
   getrow = Pmat->getrow->func_ptr;

   /* transpose Pmat's communication list. This means that PRE communication */
   /* is replaced by POST, ML_OVERWRITE is replaced by ML_ADD, and the send  */
   /* send and receive lists are swapped.                                    */

   c_info     = Pmat->getrow->pre_comm;
   Nneighbors = ML_CommInfoOP_Get_Nneighbors(c_info);
   neigh_list = ML_CommInfoOP_Get_neighbors(c_info);
   remap_leng = osize;
   Nrcv = 0;
   Nsend = 0;
   for (i = 0; i < Nneighbors; i++) {
      Nrcv  += ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nsend += ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
   }
   remap_leng = osize + Nrcv + Nsend;
   remap = (int *) ML_allocate( remap_leng*sizeof(int));
   for (i = 0; i < osize; i++) remap[i] = i;
   for (i = osize; i < osize+Nrcv+Nsend; i++) 
      remap[i] = -1;
 
   c2_info     = &(Rmat->getrow->post_comm);
   ML_CommInfoOP_Set_neighbors(c2_info, Nneighbors,
 			      neigh_list,ML_ADD,remap,remap_leng);
   ML_free(remap);
   for (i = 0; i < Nneighbors; i++) {
      Nsend      = ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
      send_list  = ML_CommInfoOP_Get_sendlist (c_info, neigh_list[i]);
      Nrcv       = ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nghost    += Nrcv;
      rcv_list   = ML_CommInfoOP_Get_rcvlist(c_info, neigh_list[i]);
      /* handle empty rows ... i.e. ghost variables not used */
      if (rcv_list != NULL) {
         for (j = 0; j < Nrcv; j++) {
            if (rcv_list[j] > Nghost2 + osize - 1)
               Nghost2 = rcv_list[j] - osize + 1;
         }
      }
 
      ML_CommInfoOP_Set_exch_info(*c2_info, neigh_list[i], Nsend, send_list,
 				 Nrcv,rcv_list);
      if (send_list != NULL) ML_free(send_list);
      if ( rcv_list != NULL) ML_free( rcv_list);
   }
   if (Nghost2 > Nghost) Nghost = Nghost2;
   if (neigh_list != NULL) ML_free(neigh_list);

   row_ptr = (int    *) ML_allocate(sizeof(int)*(Nghost+osize+1));
   colbuf  = (int    *) ML_allocate(sizeof(int)*(Nghost+osize+1));
   valbuf  = (double *) ML_allocate(sizeof(double)*(Nghost+osize+1));

   /* count the total number of nonzeros and compute */
   /* the length of each row in the transpose.       */
 
   for (i = 0; i < Nghost+osize; i++) row_ptr[i] = 0;

   N_nzs = 0;
   for (i = 0; i < isize; i++) {
      flag = getrow((ML_Operator *) data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      if (flag == 0) pr_error("ML_Transpose_Prolongator: sizes don't work\n");
      N_nzs += length;
      for (j = 0; j < length; j++)
         row_ptr[  colbuf[j] ]++;
   }

   cols    = (int    *) ML_allocate(sizeof(int   )*(N_nzs+1));
   vals    = (double *) ML_allocate(sizeof(double)*(N_nzs+1));
   if (vals == NULL) 
      pr_error("ML_Gen_Restrictor_TransP: Out of space\n");

   /* set 'row_ptr' so it points to the beginning of each row */

   sum = 0;
   for (i = 0; i < Nghost+osize; i++) {
      new_sum = sum + row_ptr[i];
      row_ptr[i] = sum;
      sum = new_sum;
   }
   row_ptr[osize+Nghost] = sum;

   /* read in the prolongator matrix and store transpose in Rmat */

   for (i = 0; i < isize; i++) {
      getrow((ML_Operator *) data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      for (j = 0; j < length; j++) {
         cols[ row_ptr[ colbuf[j] ]   ] = i;
         vals[ row_ptr[ colbuf[j] ]++ ] = valbuf[j];
      }
   }

   /* row_ptr[i] now points to the i+1th row.    */
   /* Reset it so that it points to the ith row. */

   for (i = Nghost+osize; i > 0; i--)
      row_ptr[i] = row_ptr[i-1];
   row_ptr[0] = 0;

   ML_free(valbuf);
   ML_free(colbuf);

   /* store the matrix into ML */

   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
   temp->columns = cols;
   temp->values  = vals;
   temp->rowptr  = row_ptr;
 
   ml_handle->Rmat[level].data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Init_Restrictor(ml_handle, level, level2, isize, osize, (void *) temp);
   ML_Operator_Set_ApplyFunc(Rmat, CSR_matvec);
   ML_Operator_Set_Getrow(&(ml_handle->Rmat[level]), 
                                 Nghost+osize, CSR_getrow);
   Rmat->N_nonzeros = N_nzs;
  return(1);
}
/*  SYMB_GRID::partitionBlocksNodes ****************************************
 *
 *   - do local partition
 *
 *   INPUT:
 *    - nblk: number of blocks
 *    - pnode_part[nLocalNd]: local nodes to block map (out)
 *
 */
#ifdef HAVE_ML_METIS
#ifdef __cplusplus
extern "C" {
#endif
#include "metis.h"
#ifdef __cplusplus
}
#endif
#endif

#ifdef HAVE_ML_JOSTLE
#include "jostle.h"
#endif

#ifdef HAVE_ML_PARMETIS
#include "parmetis.h"
#endif

int ML_Operator_BlockPartition(ML_Operator *matrix, int n, int *nblks,
                         int *pnode_part, ML_Partitioner which_partitioner, 
			 double *x_coord, double *y_coord, double *z_coord,int num_PDE_eqns)
{
#ifndef METIS
#define idxtype int
#endif

  idxtype *vtxdist = NULL, *tpwts = NULL, *adjncy = NULL, *xadj = NULL, *node_wt = NULL; 
  int *map = NULL, *bindx = NULL, *blks = NULL, nprocs, myid, j, ii, jj;
  int allocated = 0, row_length, itemp1, itemp2, Nrows;
  double *val = NULL; 
  int    offset = -1;

#if defined(HAVE_ML_METIS) || defined(HAVE_ML_PARMETIS) || defined(HAVE_ML_JOSTLE) 
  int     Cstyle = 0, dummy = -1;
#endif
#ifdef HAVE_ML_PARMETIS
  idxtype ncon = 1;
  float  itr = 1000., ubvec = 1.05;
#endif
#ifdef HAVE_ML_JOSTLE
  static int jostle_called = 0;
  int    output_level, nnodes, *itmp = NULL, *proc_ids;
  int    *total_nz = NULL, *total_rows = NULL, *temp_part, zz;
  int    msgtype;
  USR_REQ *request = NULL; 
  char str[80];
#endif
#if defined(HAVE_ML_METIS) || defined(HAVE_ML_PARMETIS)
  int     options[5]={0,3,1,1,0};
  int     weightflag = 0;
#endif
#ifdef HAVE_ML_ZOLTAN
  float ZoltanVersion;
  ML_Aggregate_Viz_Stats * grid_info;
  grid_info = (ML_Aggregate_Viz_Stats *) matrix->to->Grid->Grid;
#endif
  nprocs = matrix->comm->ML_nprocs;
  myid   = matrix->comm->ML_mypid;

  /* Trivial cases */

  if (*nblks == 1) {
    for( ii = 0 ; ii < n ; ii++ ) pnode_part[ii] = 0;
    return 0;
  }
  if ( (n < 1) && (which_partitioner == ML_USEMETIS) ) return 0;

  /* check what is requested versus what is compiled */

#ifndef HAVE_ML_METIS
  if ( which_partitioner == ML_USEMETIS) {
    printf("ML_partitionBlocksNodes: Metis not linked\n");
    for (ii = 0; ii < n; ii++) pnode_part[ii] = ii;
    return 1;
  }
#endif
#ifndef HAVE_ML_ZOLTAN
  if ( which_partitioner == ML_USEZOLTAN) {
    if (myid == 0) 
      printf("ML_partitionBlocksNodes: Zoltan not linked\n");
    for (ii = 0; ii < n; ii++) pnode_part[ii] = myid;
    return 1;
  }
#endif
#if !defined(HAVE_ML_PARMETIS)
  if ( which_partitioner == ML_USEPARMETIS ) {
    if (myid == 0) 
      printf("ML_partitionBlocksNodes: Parmetis is not linked\n");
    for (ii = 0; ii < n; ii++) pnode_part[ii] = myid;
    return 1;
  }
#endif
#if !(defined(HAVE_ML_JOSTLE))
  if ( which_partitioner==ML_USEJOSTLE) {
    if (myid == 0) 
      printf("ML_partitionBlocksNodes: Jostle is not linked\n");
    for (ii = 0; ii < n; ii++) pnode_part[ii] = myid;
    return 1;
  }
#endif

  /* Build adjacency graph. If METIS is used, this corresponds to the     */
  /* local matrix. If ParMetis or Jostle is used the adjacency data       */
  /* is my part of the global matrix. Zoltan does not need this structure.*/
  /* Note: matrix does not include self connections (diag entries)        */
 
  Nrows = n;
  if (which_partitioner != ML_USEZOLTAN) {
    /* count number of nonzeros */

    ii = 0;
    for( jj = 0 ; jj < n ; jj++ ) {
      ML_get_matrix_row(matrix,1,&jj,&allocated,&bindx,&val,&row_length,0);
      ii += row_length;
    }

    /* Allow additional space for matrix if we use Jostle and some */
    /* processors are idle on the next level. This space will hold */
    /* the fine grid graph from "idle" processors that will be     */
    /*  shipped to the active processors.                          */
    
#ifdef HAVE_ML_JOSTLE
    if ( (which_partitioner == ML_USEJOSTLE) && (*nblks != nprocs)) {
      total_nz  = (int *) ML_allocate(sizeof(int)*nprocs);
      total_rows= (int *) ML_allocate(sizeof(int)*nprocs);
      itmp      = (int *) ML_allocate(sizeof(int)*nprocs);
      for (jj = 0; jj < nprocs; jj++) total_nz[jj]   = 0;
      for (jj = 0; jj < nprocs; jj++) total_rows[jj] = 0;
      total_nz[myid]   = ii;
      total_rows[myid] = n;
      ML_gsum_vec_int(&total_nz,   &itmp, nprocs, matrix->comm);
      ML_gsum_vec_int(&total_rows, &itmp, nprocs, matrix->comm);
      /* in order to use Jostle contiguous storage, we must assign    */
      /* global column id's so that after moving data proc 0 has the  */
      /* lowest numbered rows followed by proc 1, etc.                */
      offset = 0; zz = 1;
      for (jj = 0 ; jj < *nblks; jj++) {
	for (j = jj; j < nprocs; j += (*nblks) ) {
	  if (j == myid) zz = 0;
	  if (zz == 1   ) offset += total_rows[j];
	}
      }

      if ( myid < *nblks ) {
	for (jj = *nblks ; jj < nprocs; jj++) 
	  if (  (jj%(*nblks)) == myid ) {
	    ii    += total_nz[jj]; 
	    Nrows += total_rows[jj];
	  }
      }
    }
#endif

    adjncy = (idxtype *) ML_allocate( (ii+1) * sizeof(idxtype) );
    if (adjncy == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");
    xadj = (idxtype *) ML_allocate( (Nrows+1) * sizeof(idxtype) );
    if (xadj == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");

    /* fill graph */

    if (which_partitioner == ML_USEMETIS) {
      ii = 0;
      for( jj = 0 ; jj < n ; jj++ ) {
	xadj[jj] = ii;
	ML_get_matrix_row(matrix,1,&jj,&allocated,&bindx,&val,&row_length,0);
	for (j = 0; j < row_length; j++) {
	  /* do not record diagonal or off-processor entries */
	  if (( bindx[j] != jj ) && (bindx[j] < n))
	    adjncy[ii++] = bindx[j];
	}
      }
      xadj[n] = ii;
    }
    else {
      ML_create_unique_id(n, &map, matrix->getrow->pre_comm, matrix->comm,
			  offset);
      ii = 0;
      for( jj = 0 ; jj < n ; jj++ )   {
	xadj[jj] = ii;
	ML_get_matrix_row(matrix,1,&jj,&allocated,&bindx,&val,&row_length,0);
	for (j = 0; j < row_length; j++)     {
	  if ( bindx[j] != jj ) {
	    adjncy[ii++] = map[bindx[j]];
	  }
	}
      }
      xadj[n] = ii;
    }
  }

#ifdef HAVE_ML_JOSTLE

  /* Ship rows from "idle" processors to "active" processors. */

  if ( (which_partitioner == ML_USEJOSTLE) && (*nblks != nprocs)) {
    if ( myid < *nblks ) {
      ii = n;
      zz = 0;
      
      request = (USR_REQ *) ML_allocate((2+(nprocs - *nblks)/(*nblks))*
					sizeof(USR_REQ)); /* estimated number*/
                                                          /* of rcv neighbors*/
      for (jj = *nblks ; jj < nprocs; jj++) {
	if (  (jj%(*nblks)) == myid ) {
	  /* receive xadj and adjncy */
	  if (total_rows[jj] != 0) {
	    msgtype = 95539;
	    matrix->comm->USR_irecvbytes((void *) &(xadj[ii+1]), 
					 total_rows[jj]*sizeof(int), &jj,
					 &msgtype, matrix->comm->USR_comm, 
					 request+zz);
	    matrix->comm->USR_cheapwaitbytes((void *) &(xadj[ii+1]), 
					 total_rows[jj]*sizeof(int), &jj,
					 &msgtype, matrix->comm->USR_comm, 
					 request+zz);
	    
	    /* fix xadj so that it is correct */
	    for (j = ii+1; j <= ii+total_rows[jj]; j++) {
	      xadj[j] += xadj[ii];
	    }
	    msgtype = 94539;
	    matrix->comm->USR_irecvbytes((void *) &(adjncy[xadj[ii]]), 
					 total_nz[jj]*sizeof(int), &jj,
					 &msgtype, matrix->comm->USR_comm, 
					 request+zz);
	    matrix->comm->USR_cheapwaitbytes((void *) &(adjncy[xadj[ii]]), 
					 total_nz[jj]*sizeof(int), &jj,
					 &msgtype, matrix->comm->USR_comm, 
					 request+zz);
	    zz++;
	    ii += total_rows[jj];
	  } /* total_rows[jj] != 0) */
	} /* (jj%nblks)) == myid */
      } /* for (jj = *nblks; jj < nprocs; jj++) */
    } /* if ( myid < *nblks ) */
    else {
      request = (USR_REQ *) ML_allocate(sizeof(USR_REQ));
                                                         
      /* send xadj and adjncy */
      if (n != 0) {
	msgtype = 95539;
	matrix->comm->USR_sendbytes((void *) &(xadj[1]), n*sizeof(int), 
				    myid%(*nblks), msgtype, 
				    matrix->comm->USR_comm);
	msgtype = 94539;
	matrix->comm->USR_sendbytes((void *) adjncy,
				    total_nz[myid]*sizeof(int), 
				    myid%(*nblks), msgtype, 
				    matrix->comm->USR_comm);
      } /* if (n != 0) */
      Nrows = 0;
    } /* if ( myid < *nblks ) */
    ML_free(total_nz);
    temp_part = (int *) ML_allocate( (Nrows+1)*sizeof(int));
  } /* (which_partitioner == ML_USEPARMETIS) && (*nblks != nprocs) */
  else temp_part = pnode_part;
#endif


  switch(which_partitioner) {
  case ML_USEZOLTAN:
#ifdef HAVE_ML_ZOLTAN
    if (Zoltan_Initialize(0, NULL, &ZoltanVersion) == ZOLTAN_OK)
      if (matrix->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 4) {
        printf("Repartitioning using Zoltan %3.2f\n",ZoltanVersion);
      }
    if (ML_DecomposeGraph_with_Zoltan(matrix, *nblks, pnode_part, NULL,
				      x_coord, y_coord, z_coord, matrix->to->levelnum,
                                      grid_info->zoltan_type,
                                      grid_info->zoltan_estimated_its,
				      grid_info->zoltan_timers,
                                      grid_info->smoothing_steps,
                                      num_PDE_eqns
                                      ) < 0)
      for (ii = 0; ii < n; ii++) pnode_part[ii] = myid;
#endif
    break;
  case ML_USEMETIS:
#ifdef HAVE_ML_METIS
    if (matrix->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 4)
      printf("Repartitioning using METIS\n");
    if (*nblks < 8)
      METIS_PartGraphRecursive( &n, xadj, adjncy, NULL, NULL, &weightflag,
				&Cstyle, nblks, options, &dummy, pnode_part );
    else
      METIS_PartGraphKway( &n, xadj, adjncy, NULL, NULL, &weightflag, &Cstyle,
			   nblks, options, &dummy, pnode_part );
#endif

    /* It looks like this handles an empty block */

    blks = (int *) ML_allocate((*nblks+1)*sizeof(int));
    if (blks == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");
    for (j = 0; j < *nblks; j++) blks[j] = -1;
    for (j = 0; j < n; j++) blks[pnode_part[j]] = -2;
    ii = 0;
    for (j = 0; j < *nblks; j++) {
      if ( blks[j] == -2) {
	blks[j] = ii;
	ii++;
      }
    }
    for (j = 0; j < n; j++) pnode_part[j] = blks[pnode_part[j]];
    *nblks = ii;
    break;
  case ML_USEPARMETIS:
    tpwts   = (idxtype *) ML_allocate(sizeof(idxtype)*ML_max(*nblks,nprocs));
    vtxdist = (idxtype *) ML_allocate(sizeof(idxtype)*(nprocs+1));
    for (ii = 0; ii <= nprocs; ii++) vtxdist[ii] = 0;
    vtxdist[myid] = n;
    ML_gsum_vec_int(&vtxdist,&tpwts,nprocs,matrix->comm);

    itemp1 = 0;
    for (ii = 0; ii < nprocs ; ii++) {
      itemp2 = vtxdist[ii];
      vtxdist[ii] = itemp1;
      itemp1 += itemp2;
    }
    vtxdist[nprocs] = itemp1;
    
    for (ii = 0; ii < n; ii++) pnode_part[ii] = matrix->comm->ML_mypid;

#ifdef HAVE_ML_PARMETIS
    node_wt = (idxtype *) ML_allocate( (Nrows+1) * sizeof(idxtype) );
    for (j = 0; j < Nrows; j++) {
      node_wt[j] = xadj[j+1] - xadj[j] + 1;
    }
    weightflag = 2;  /*node weights*/

    if (matrix->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 4)
      printf("Repartitioning using ParMETIS3x\n");
    ParMETIS_V3_AdaptiveRepart( vtxdist,xadj,adjncy, node_wt,
				NULL, NULL, &weightflag,
				&Cstyle, &ncon, nblks, NULL, &ubvec, 
				&itr, options, &dummy, pnode_part,
				&(matrix->comm->USR_comm));
#endif
    break;

#ifdef HAVE_ML_JOSTLE
  case ML_USEJOSTLE:
    if (matrix->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 4)
      printf("Repartitioning using jostle\n");
    if (jostle_called == 0) {
      pjostle_init(&nprocs, &myid);
      jostle_called = 1;
    }
    output_level = 2;
    nnodes = ML_Comm_GsumInt(matrix->comm,n);
    jostle_env("format = contiguous");
    proc_ids = (int *) ML_allocate( nprocs*sizeof(int));
    if (itmp == NULL) 
      itmp  = (int *) ML_allocate( nprocs*sizeof(int));
    for (j = 0; j < nprocs; j++) proc_ids[j] = 0;
    proc_ids[myid] = Nrows;
    ML_gsum_vec_int(&proc_ids, &itmp, nprocs, matrix->comm);
    node_wt = (idxtype *) ML_allocate( (Nrows+1) * sizeof(idxtype) );

    for (j = 0; j < Nrows; j++) {
      xadj[j] = xadj[j+1] - xadj[j];
      node_wt[j] = xadj[j] + 1;
    }
    ML_free(itmp);
    jostle_env("matching = local");
    /*
    jostle_env("threshold = 2");
    jostle_env("imbalance = 0");
    */
    /*    jostle_env("connect on"); */
    /*    jostle_env("check = on"); */

    dummy = 0;

    /* check if some processors are idle */

    if ( *nblks < nprocs) {
      sprintf(str,"idle = %d",nprocs - *nblks);
      jostle_env(str); 
    }
    pjostle(&nnodes,&Cstyle,&Nrows,&dummy,proc_ids,xadj,node_wt,temp_part,
	    /*    pjostle(&nnodes,&Cstyle,&Nrows,&dummy,proc_ids,xadj,(int *) NULL,temp_part, */
	    &(xadj[Nrows]),adjncy,(int *) NULL, &nprocs, (int *) NULL,
	    &output_level, &dummy, (double *) NULL);

    /* ship over temp pnode part */
    if ( (which_partitioner == ML_USEPARMETIS) && (*nblks != nprocs)) {
      msgtype = 93539;
      if ( myid < *nblks) {
	for (jj = 0; jj < n; jj++) pnode_part[jj] = temp_part[jj];
	ii = n;
	for (jj = *nblks ; jj < nprocs; jj++) {
	  if (  (jj%(*nblks)) == myid ) {
	    /* send pnode part */
	    if (total_rows[jj] != 0) {
	      matrix->comm->USR_sendbytes((void *) &(temp_part[ii]),
				    total_rows[jj]*sizeof(int), 
				    jj, msgtype, 
				    matrix->comm->USR_comm);
	      ii += total_rows[jj];
	    }
	  }
	}
      }
      else {
	/*receive pnode part stuff */
	zz = myid%(*nblks);
	if (n != 0) {
	  matrix->comm->USR_irecvbytes((void *) pnode_part,
				     n*sizeof(int), &zz,
				     &msgtype, matrix->comm->USR_comm, 
				     request);
	  matrix->comm->USR_cheapwaitbytes((void *) pnode_part,
					 n*sizeof(int), &zz,
					 &msgtype, matrix->comm->USR_comm, 
					 request);
	}
      } /* end else */
      ML_free(temp_part);
      ML_free(total_rows);
    }
    if (request != NULL) ML_free(request);

    break;
#endif
/* -------------------------------------- */
  default:
    printf("ML_partitionBlocksNodes: Unknown partitioner requested %d\n",
	   which_partitioner);
  }
  if (map     != NULL) ML_free(map);
  if (val     != NULL) ML_free(val);
  if (bindx   != NULL) ML_free(bindx);
  if (blks    != NULL) ML_free(blks);
  if (vtxdist != NULL) ML_free(vtxdist);
  if (tpwts   != NULL) ML_free(tpwts);
  if (adjncy  != NULL) ML_free(adjncy);
  if (xadj    != NULL) ML_free(xadj);
  if (node_wt != NULL) ML_free(node_wt);

  return 0;
}
/* ******************************************************************** */
/* Take the Transpose of an ML_Operator and stick it in a new matrix    */
/* -------------------------------------------------------------------- */

int ML_Operator_Transpose(ML_Operator *Amat, ML_Operator *Amat_trans )
{

   int  *row_ptr, *colbuf, *cols;
   int isize, osize, i, j, N_nzs, flag, length, sum, new_sum;
   int Nneighbors, *neigh_list, *send_list, *rcv_list, Nsend, Nrcv;
   void *data = NULL;
   double *valbuf, *vals;
   int (*getrow)(ML_Operator* , int , int *, int , int *, double *, int *) = NULL;
   struct ML_CSR_MSRdata *temp;
   int Nghost = 0, Nghost2 = 0;
   int *remap, remap_leng;
   ML_CommInfoOP *c_info, **c2_info;

   temp = (struct ML_CSR_MSRdata *) Amat->data;
   isize = Amat->outvec_leng;
   osize = Amat->invec_leng;
   data   = (void *) Amat;
   getrow = Amat->getrow->func_ptr;

   /* transpose Amat's communication list. This means that PRE communication */
   /* is replaced by POST, ML_OVERWRITE is replaced by ML_ADD, and the send  */
   /* send and receive lists are swapped.                                    */

   c_info     = Amat->getrow->pre_comm;
   Nneighbors = ML_CommInfoOP_Get_Nneighbors(c_info);
   neigh_list = ML_CommInfoOP_Get_neighbors(c_info);
   remap_leng = osize;
   Nrcv = 0;
   Nsend = 0;
   for (i = 0; i < Nneighbors; i++) {
      Nrcv  += ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nsend += ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
   }
   remap_leng = osize + Nrcv + Nsend;
   remap = (int *) ML_allocate( remap_leng*sizeof(int));
   for (i = 0; i < osize; i++) remap[i] = i;
   for (i = osize; i < osize+Nrcv+Nsend; i++) 
      remap[i] = -1;
 
   c2_info     = &(Amat_trans->getrow->post_comm);
   ML_CommInfoOP_Set_neighbors(c2_info, Nneighbors,
 			      neigh_list,ML_ADD,remap,remap_leng);
   ML_free(remap);
   for (i = 0; i < Nneighbors; i++) {
      Nsend      = ML_CommInfoOP_Get_Nsendlist(c_info, neigh_list[i]);
      send_list  = ML_CommInfoOP_Get_sendlist (c_info, neigh_list[i]);
      Nrcv       = ML_CommInfoOP_Get_Nrcvlist (c_info, neigh_list[i]);
      Nghost    += Nrcv;
      rcv_list   = ML_CommInfoOP_Get_rcvlist(c_info, neigh_list[i]);
      /* handle empty rows ... i.e. ghost variables not used */
      if (rcv_list != NULL) {
	for (j = 0; j < Nrcv; j++) {
            if (rcv_list[j] > Nghost2 + osize - 1)
               Nghost2 = rcv_list[j] - osize + 1;
         }
      }
 
      ML_CommInfoOP_Set_exch_info(*c2_info, neigh_list[i], Nsend, send_list,
 				 Nrcv,rcv_list);
      if (send_list != NULL) ML_free(send_list);
      if ( rcv_list != NULL) ML_free( rcv_list);
   }
   if (Nghost2 > Nghost) Nghost = Nghost2;
   if (neigh_list != NULL) ML_free(neigh_list);

   row_ptr = (int    *) ML_allocate(sizeof(int)*(Nghost+osize+1));
   colbuf  = (int    *) ML_allocate(sizeof(int)*(Nghost+osize+1));
   valbuf  = (double *) ML_allocate(sizeof(double)*(Nghost+osize+1));

   /* count the total number of nonzeros and compute */
   /* the length of each row in the transpose.       */
 
   for (i = 0; i < Nghost+osize; i++) row_ptr[i] = 0;

   N_nzs = 0;
   for (i = 0; i < isize; i++) {
      flag = getrow((ML_Operator *) data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      if (flag == 0) perror("ML_Transpose_Prolongator: sizes don't work\n");
      N_nzs += length;
      for (j = 0; j < length; j++)
         row_ptr[  colbuf[j] ]++;
   }

   cols    = (int    *) ML_allocate(sizeof(int   )*(N_nzs+1));
   vals    = (double *) ML_allocate(sizeof(double)*(N_nzs+1));
   if (vals == NULL) 
      pr_error("ML_Gen_Restrictor_TransP: Out of space\n");

   /* set 'row_ptr' so it points to the beginning of each row */

   sum = 0;
   for (i = 0; i < Nghost+osize; i++) {
      new_sum = sum + row_ptr[i];
      row_ptr[i] = sum;
      sum = new_sum;
   }
   row_ptr[osize+Nghost] = sum;

   /* read in the prolongator matrix and store transpose in Amat_trans */

   for (i = 0; i < isize; i++) {
      getrow((ML_Operator *) data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
      for (j = 0; j < length; j++) {
         cols[ row_ptr[ colbuf[j] ]   ] = i;
         vals[ row_ptr[ colbuf[j] ]++ ] = valbuf[j];
      }
   }

   /* row_ptr[i] now points to the i+1th row.    */
   /* Reset it so that it points to the ith row. */

   for (i = Nghost+osize; i > 0; i--)
      row_ptr[i] = row_ptr[i-1];
   row_ptr[0] = 0;

   ML_free(valbuf);
   ML_free(colbuf);

   /* store the matrix into ML */

   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
   temp->columns = cols;
   temp->values  = vals;
   temp->rowptr  = row_ptr;
   Amat_trans->N_nonzeros = N_nzs;
   Amat_trans->data_destroy = ML_CSR_MSRdata_Destroy;
   
   ML_Operator_Set_ApplyFuncData(Amat_trans, isize, osize, 
                                  temp, osize, NULL, 0);
   ML_Operator_Set_ApplyFunc(Amat_trans, CSR_matvec);
   ML_Operator_Set_Getrow(Amat_trans, Nghost+osize, CSR_getrow);

  return(1);
}

/************************************************************************/
/* Take a matrix that is effectively partitioned by columns and         */
/* transform it into one that is partitioned by rows. The original      */
/* matrix was most likely created by transposing a matrix partitioned   */
/* by row.                                                              */
/*----------------------------------------------------------------------*/

int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans)
{
 
  ML_Operator *eye1;
 
  eye1 = ML_Operator_Create(A->comm);
 
  ML_Operator_Set_ApplyFuncData(eye1, A->invec_leng, A->invec_leng,
            NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye1, A->invec_leng, eye_getrows);
 
  ML_2matmult(A, eye1, Atrans, ML_CSR_MATRIX);

  ML_Operator_Destroy(&eye1);

 
  return 1;
}  


#ifdef new2row
int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans)
{
  int         max_per_proc;
  ML_Operator *Acomm, *tptr;
 
  if (A->getrow->use_loc_glob_map == ML_YES)
     pr_error("ML_Operator_ColPartition2RowPartition: Matrix already has local"
              "column indices mapped to global indices\n");
  if (A->getrow->pre_comm != NULL)
     pr_error("ML_Operator_ColPartition2RowPartiion: Matrix has a"
              "pre-communication structure?\n");
 
  ML_create_unique_col_id(A->invec_leng, &(A->getrow->loc_glob_map),
                           NULL, &max_per_proc, A->comm);
 
  A->getrow->use_loc_glob_map = ML_YES;
 
  if (A->getrow->post_comm != NULL)
      ML_exchange_rows( A, &Acomm, A->getrow->post_comm);
  else Acomm = A;
 
  ML_back_to_csrlocal(Acomm, Atrans, max_per_proc);
 
  ML_free(A->getrow->loc_glob_map); A->getrow->loc_glob_map = NULL;
  A->getrow->use_loc_glob_map = ML_NO;
 
  if (A->getrow->post_comm != NULL) {
      tptr = Acomm;
      while ( (tptr!= NULL) && (tptr->sub_matrix != A))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(Acomm);
      ML_Operator_Destroy(&Acomm);
   }
 
  return 1;
} 
#endif


/************************************************************************/
/* Getrow function for the identity matrix.                             */
/*----------------------------------------------------------------------*/

int eye_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
                int allocated_space, int columns[], double values[],
				int row_lengths[])
{
   int    i;

   
   if (allocated_space < N_requested_rows) {
     ML_avoid_unused_param( data);
     return(0);
   }

   for (i = 0; i < N_requested_rows; i++) {
      row_lengths[i] = 1;
      columns[i]     = requested_rows[i];
      values[i]      = 1.;
   }
   return(1);
}

/************************************************************************/
/* Matvec function for the identity matrix.                             */
/*----------------------------------------------------------------------*/

int eye_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{
  int i;

  if (ilen == -57) ML_avoid_unused_param( Amat_in);


  for (i = 0; i < olen; i++) ap[i] = p[i];

  return(1);
}
/************************************************************************/
/* Take the transpose of an ML_Operator and realign resulting matrix    */
/* so that it is partitioned by rows.                                   */
/*----------------------------------------------------------------------*/

int ML_Operator_Transpose_byrow(ML_Operator *A, ML_Operator *Atrans)
{
  ML_Operator *temp;

  temp = ML_Operator_Create(A->comm);
  ML_Operator_Transpose(A, temp);
  ML_Operator_ColPartition2RowPartition(temp, Atrans);
  ML_Operator_Destroy(&temp);
  return 1;
}
#include "ml_utils.h"
#include "ml_xyt.h"
int ML_Operator_Dump(ML_Operator *Ke, double *x, double *rhs,
		     char *istr, int print_matrix)	
{
  double *global_nodes, *global_rows, *colVal = NULL;
  int    N_nodes, node_offset, row_offset;
  int *colInd = NULL, i, j, ncnt, allocated = 0;
  char str[80];
  FILE *fid;
  ML_Comm *comm;
  int Nnodes_global, Nrows_global;
  int Nghost_nodes;
  int Nrows;
  

  comm = Ke->comm;
  if (Ke->getrow->pre_comm == NULL) Nghost_nodes = 0;
  else {
    if (Ke->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Ke->getrow->pre_comm);
    Nghost_nodes = Ke->getrow->pre_comm->total_rcv_length;
  }


  N_nodes = Ke->invec_leng;
  node_offset = ML_gpartialsum_int(N_nodes, comm);
  Nnodes_global = N_nodes;
  ML_gsum_scalar_int(&Nnodes_global, &i, comm);

  Nrows = Ke->outvec_leng;
  row_offset = ML_gpartialsum_int(Nrows, comm);
  Nrows_global = Nrows;
  ML_gsum_scalar_int(&Nrows_global, &i, comm);

  global_nodes  =(double *) ML_allocate(sizeof(double)*(N_nodes+Nghost_nodes));
  global_rows   =(double *) ML_allocate(sizeof(double)*(Nrows));

  for (i = 0 ; i < N_nodes; i++) global_nodes[i] = (double) (node_offset + i);
  for (i = 0 ; i < Nrows; i++) global_rows[i] = (double) (row_offset + i);

  for (i = 0 ; i < Nghost_nodes; i++) global_nodes[i+N_nodes] = -1;

  ML_exchange_bdry(global_nodes,Ke->getrow->pre_comm, 
 		 Ke->invec_leng,comm,ML_OVERWRITE,NULL);

  /* spit out Ke  */

  if (print_matrix) {
    sprintf(str,"%s_mat.%d",istr,comm->ML_mypid);
    fid = fopen(str,"w");
    for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Ke, 1, &i, &allocated, &colInd, &colVal,
                        &ncnt, 0);

      for (j = 0; j < ncnt; j++) {
	if (colVal[j] != 0.0) {
	  fprintf(fid,"%5d %5d %20.13e\n",(int) global_rows[i]+1,
		  (int) global_nodes[colInd[j]]+1, colVal[j]);
	}
      }
    }
    fclose(fid);
    ML_free(colVal); ML_free(colInd);
  }


  /* spit out x */

  if (x != NULL) {
    sprintf(str,"%s_xxx.%d",istr,comm->ML_mypid);
    fid = fopen(str,"w");
    for (i = 0; i < Ke->invec_leng; i++) {
      fprintf(fid,"%5d %20.13e\n",(int) global_nodes[i]+1,x[i]);
    }
    fclose(fid);
  }

  /* spit out rhs */

  if (rhs != NULL) {
    sprintf(str,"%s_rhs.%d",istr,comm->ML_mypid);
    fid = fopen(str,"w");
    for (i = 0; i < Ke->outvec_leng; i++) {
      fprintf(fid,"%5d %20.13e\n",(int) global_rows[i]+1,rhs[i]);
    }
    fclose(fid);
  }


  ML_free(global_nodes);
  ML_free(global_rows);
  return 0;
}

int ML_Operator_Getrow_Diag(ML_Operator *Amat, double **diagonal)
{
  int allocated_space, *cols, i, j, n, Nghost;
  double *vals, *tdiag;
   if (Amat->diagonal == NULL) 
   {
      if (Amat->getrow->func_ptr == NULL) 
         pr_error("Error(ML_Operator_Getrow_Diag): No getrow function\n");
      else 
      {
         Nghost = ML_CommInfoOP_Compute_TotalRcvLength(Amat->getrow->pre_comm);
         allocated_space = 30;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         tdiag = (double *) ML_allocate((Amat->outvec_leng+Nghost+1)*sizeof(double));
         for (i = 0; i < Amat->outvec_leng; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&n) == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols); 
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL)
               {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < n; j++) 
               if (cols[j] == i) tdiag[i] = vals[j];
         }
         if ( Amat->getrow->pre_comm != NULL )
           ML_exchange_bdry(tdiag,Amat->getrow->pre_comm,Amat->getrow->Nrows, 
                            Amat->comm, ML_OVERWRITE,NULL);

         ML_free(cols); ML_free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         ML_free(tdiag);
      } 
   }
   ML_DVector_GetDataPtr( Amat->diagonal, diagonal);
   return 0;
}

/*******************************************************************
 *  Take an ML_Operator and make a new copy of its data corresponding
 *  to a scaled ML_Operator. 
 *******************************************************************/
ML_Operator *ML_Operator_ExplicitlyScale(ML_Operator *matrix,
					 double scalar)
{
  int i, k, Nrows, Nnz, allocated = 0, *columns = NULL, row_length;
  int *row_ptr, *col_ptr;
  double *values = NULL;
  double  *val_ptr;
  struct ML_CSR_MSRdata *temp;
  ML_Operator *new_matrix;

  /* first count how many nonzeros are in the old matrix */

  Nnz = 0;
  if (matrix->getrow == NULL) return(NULL);
  Nrows = matrix->getrow->Nrows;

  for (i = 0 ; i < Nrows; i++) {
    ML_get_matrix_row(matrix, 1, &i, &allocated, &columns, &values,
                        &row_length, 0);
    Nnz += row_length;
  }

  /* allocate space */

  row_ptr = (int   *) ML_allocate(sizeof(int  )*(Nrows + 1));
  col_ptr = (int   *) ML_allocate(sizeof(int  )*(Nnz + 1));
  val_ptr = (double *) ML_allocate(sizeof(double)*(Nnz + 1));
  temp    = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  
  /* getrow everything and scale it */

   row_ptr[0] = 0;
   Nnz = 0;
   for (i = 0 ; i < Nrows; i++) {
     ML_get_matrix_row(matrix, 1, &i, &allocated, &columns, &values,
		       &row_length, 0);
     for (k = 0; k < row_length; k++) {
         val_ptr[Nnz  ] = scalar*values[k];
         col_ptr[Nnz++] = columns[k];
     }
     row_ptr[i+1] = Nnz;
   }
   temp->rowptr  = row_ptr;
   temp->columns = col_ptr;
   temp->values  = val_ptr;


   /* Get rid of the old data pointer */

   new_matrix = ML_Operator_Create(matrix->comm);
   /*   ML_Operator_Init(new_matrix, matrix->comm); */

   ML_Operator_Set_ApplyFuncData(new_matrix,matrix->invec_leng, 
				 matrix->outvec_leng,temp,
				 matrix->matvec->Nrows, CSR_matvec,
				 matrix->from_an_ml_operator);
   ML_Operator_Set_Getrow(new_matrix,matrix->getrow->Nrows,CSR_getrow);
   ML_CommInfoOP_Clone(&(new_matrix->getrow->pre_comm),matrix->getrow->pre_comm);
   new_matrix->data_destroy   = ML_CSR_MSRdata_Destroy;
   if (values  != NULL) ML_free(values);
   if (columns != NULL) ML_free(columns);

   return new_matrix;
}

/*******************************************************************
 *  Take an ML_Operator and using getrow() make a new copy of the
 *  matrix in CSR format using single precision numbers. Then,
 *  get rid of the data in the old matrix (by calling the 
 *  data destroy function) and replace the old data with the
 *  new data.
 *******************************************************************/
int ML_Operator_ChangeToSinglePrecision(ML_Operator *matrix)
{
  int i, k, Nrows, Nnz, allocated = 0, *columns = NULL, row_length;
  int *row_ptr, *col_ptr;
  double *values = NULL;
  float  *val_ptr;
  struct ML_CSR_MSRdata *temp;

  if (ML_Use_LowMemory() != ML_TRUE)
  return 1;

  /* Do not do anything if there is no destroy function for    */
  /* this matrix as we have no way to destroy the old data.    */
  if ((matrix->data_destroy == NULL) || (matrix->data == NULL))
    return 1;


  /* first count how many nonzeros are in the old matrix */

  Nnz = 0;
  if (matrix->getrow == NULL) return(1);
  Nrows = matrix->getrow->Nrows;

  for (i = 0 ; i < Nrows; i++) {
    ML_get_matrix_row(matrix, 1, &i, &allocated, &columns, &values,
                        &row_length, 0);
    Nnz += row_length;
  }

  /* allocate space */

  row_ptr = (int   *) ML_allocate(sizeof(int  )*(Nrows + 1));
  col_ptr = (int   *) ML_allocate(sizeof(int  )*(Nnz + 1));
  val_ptr = (float *) ML_allocate(sizeof(float)*(Nnz + 1));
  temp    = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  
  /* getrow everything and copy it to single precision */

   row_ptr[0] = 0;
   Nnz = 0;
   for (i = 0 ; i < Nrows; i++) {
     ML_get_matrix_row(matrix, 1, &i, &allocated, &columns, &values,
		       &row_length, 0);
     for (k = 0; k < row_length; k++) {
         val_ptr[Nnz  ] = (float) values[k];
         col_ptr[Nnz++] = columns[k];
     }
     row_ptr[i+1] = Nnz;
   }
   temp->rowptr = row_ptr;
   temp->columns = col_ptr;
   temp->values = (double *) val_ptr;

   /* Get rid of the old data pointer */

   if ((matrix->data_destroy != NULL) && (matrix->data != NULL)) {
      matrix->data_destroy(matrix->data);
      matrix->data = NULL;
   }

   ML_Operator_Set_ApplyFuncData(matrix,matrix->invec_leng, 
				 matrix->outvec_leng,temp,
				 matrix->matvec->Nrows, sCSR_matvec,
				 matrix->from_an_ml_operator);

   ML_Operator_Set_Getrow(matrix,matrix->getrow->Nrows,sCSR_getrows);
   matrix->data_destroy   = ML_CSR_MSRdata_Destroy;
   if (values  != NULL) ML_free(values);
   if (columns != NULL) ML_free(columns);

   return 0;
}
int ML_Operator_ChangeToChar(ML_Operator *matrix)
{
  int i, k, Nrows, Nnz, allocated = 0, *columns = NULL, row_length;
  int *row_ptr, *col_ptr;
  double *values = NULL;
  char  *val_ptr;
  struct ML_CSR_MSRdata *temp;

  if (ML_Use_LowMemory() != ML_TRUE)
    return(1);

  /* Do not do anything if there is no destroy function for    */
  /* this matrix as we have no way to destroy the old data.    */
  if ((matrix->data_destroy == NULL) || (matrix->data == NULL))
    return 1;

  /* first count how many nonzeros are in the old matrix */

  Nnz = 0;
  if (matrix->getrow == NULL) return(1);
  Nrows = matrix->getrow->Nrows;

  for (i = 0 ; i < Nrows; i++) {
    ML_get_matrix_row(matrix, 1, &i, &allocated, &columns, &values,
                        &row_length, 0);
    Nnz += row_length;
  }

  /* allocate space */

  row_ptr = (int   *) ML_allocate(sizeof(int  )*(Nrows + 1));
  col_ptr = (int   *) ML_allocate(sizeof(int  )*(Nnz + 1));
  val_ptr = (char  *) ML_allocate(sizeof(char)*(Nnz + 1));
  temp    = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  
  /* getrow everything and copy it to single precision */

   row_ptr[0] = 0;
   Nnz = 0;
   for (i = 0 ; i < Nrows; i++) {
     ML_get_matrix_row(matrix, 1, &i, &allocated, &columns, &values,
		       &row_length, 0);
     for (k = 0; k < row_length; k++) {
       if      (values[k] == -1.)  val_ptr[Nnz] = (int) 2;
       else if (values[k] == 1.)   val_ptr[Nnz] = (int) 1;
       else if (values[k] == 0.)   val_ptr[Nnz] = (int) 0;
       else pr_error("ML_Operator_ChangeToChar(%d): T(%d,%d) = %e! It must be 1,-1,or 0!!!",matrix->comm->ML_mypid,i,columns[k],values[k]);
       col_ptr[Nnz++] = columns[k];
     }
     row_ptr[i+1] = Nnz;
   }
   temp->rowptr = row_ptr;
   temp->columns = col_ptr;
   temp->values = (double *) val_ptr;

   /* Get rid of the old data pointer */

   if ((matrix->data_destroy != NULL) && (matrix->data != NULL)) {
      matrix->data_destroy(matrix->data);
      matrix->data = NULL;
   }

   ML_Operator_Set_ApplyFuncData(matrix,matrix->invec_leng, 
				 matrix->outvec_leng,temp,
				 matrix->matvec->Nrows, cCSR_matvec,
				 matrix->from_an_ml_operator);

   ML_Operator_Set_Getrow(matrix,matrix->getrow->Nrows,cCSR_getrows);
   matrix->data_destroy   = ML_CSR_MSRdata_Destroy;
   if (values  != NULL) ML_free(values);
   if (columns != NULL) ML_free(columns);

   return 0;
}

/*******************************************************************
 *  Take an ML_Operator single precision or char CSR matrix in Pmat and use it 
 *  to implicitly define a transpose matvec that is stored in Rmat. This 
 *  implies that Rmat will not have a getrow() function ... so it is important
 *  that this function only be called on matrices where getrow() is
 *  no longer needed. If there is already data in Rmat, get rid of
 *  it. 
 *
 *  Note: The main tricky thing about this function is that we must
 *  set up a post commuication data structure. This post communication
 *  structure is fairly fragile. We also allow for the possibility
 *  that the post communciation structure is already set up. In
 *  this case, we leave the existing communication structure in 
 *  Rmat.
 *******************************************************************/
int ML_Operator_ImplicitTranspose(ML_Operator *Rmat, 
				  ML_Operator *Pmat,
				  int PostCommAlreadySet)
{

  if (ML_Use_LowMemory() != ML_TRUE)
    return 1;

  if ( (Pmat == NULL) || (Rmat == NULL)) return 1;

  if ( Pmat->getrow == NULL) return 1;

  if ( (Pmat->getrow->func_ptr != sCSR_getrows) &&
       (Pmat->getrow->func_ptr != cCSR_getrows)) return 1;

  if (PostCommAlreadySet == ML_FALSE) {
    if (Rmat->getrow->post_comm != NULL)
      ML_CommInfoOP_Destroy(&(Rmat->getrow->post_comm));
    ML_CommInfoOP_TransComm(Pmat->getrow->pre_comm,&(Rmat->getrow->post_comm),
			    Pmat->invec_leng);
  }

  if (Pmat->getrow->func_ptr == sCSR_getrows)
    ML_Operator_Set_ApplyFuncData(Rmat, Pmat->outvec_leng,
				Pmat->invec_leng, 
				Pmat->data, -1, sCSR_trans_matvec, 0);
  else
    ML_Operator_Set_ApplyFuncData(Rmat, Pmat->outvec_leng,
				Pmat->invec_leng, 
				Pmat->data, -1, cCSR_trans_matvec, 0);

  Rmat->getrow->func_ptr = NULL;
  Rmat->data_destroy = NULL;
  return 0;
}
/*******************************************************************
 * Build a communication structure for newly received matrix rows.
 * This entails looking at the new rows and deciding which columns
 * (global ids) are not already within this overlapped matrix. These
 * ghost columns are then appended to extern[] and a communication object
 * is built to obtain them.
 *******************************************************************/
int ML_build_overlapped_pre_comm(ML_Operator *tempA, ML_CommInfoOP
				 *old_comm, int max_per_proc,
				 int *hash_list, int hash_length, 
				 int *hash_used, int old_Nrows, 
				 int *Nexternal, int *external[],
				 int *Nexternal_allocated)
{
  int    i, j, k, NGhost, index;
  int allocated = 0, *column = NULL, row_length, current, proc_id;
  double *val = NULL;
  int    oldNexternal, *temp;

  oldNexternal = *Nexternal;
  proc_id      = tempA->comm->ML_mypid;
  NGhost       = ML_CommInfoOP_Compute_TotalRcvLength(old_comm);

  if ( *Nexternal_allocated - oldNexternal < 2*NGhost) {
    /* estimate space for new externals.      */
    /* if not enough, we allocate more later. */
    k = *Nexternal_allocated;
    *Nexternal_allocated = oldNexternal + (NGhost+5)*5;
    temp = (int *) ML_allocate(*Nexternal_allocated*sizeof(int));
    if (temp==NULL) perror("ML_build_overlapped_pre_comm: Not enough space\n");
    for (i = 0; i < *Nexternal; i++) temp[i] = (*external)[i];
    if (k != 0) ML_free(*external);
    *external = temp;
  }

  /* look at the newly received matrix rows associated with the ghost */
  /* unknowns and see who is needed from other processors             */

  for (i = old_Nrows; i < NGhost+old_Nrows ; i++) {
    ML_get_matrix_row(tempA, 1, &i, &allocated, &column, &val,
                        &row_length, 0);

    for (j = 0; j < row_length; j++) {
      current = column[j];

      /* check if a nonlocal guy */
      if ( (current < max_per_proc*proc_id) ||
	   (current >= max_per_proc*(proc_id+1))) {
	/* now check if it is already in the hash table       */
	/* if not, add to hash table and mark gid as external */
	ML_hash_it( current, hash_list, hash_length, hash_used, &index);
	if ( hash_list[index] == -1 ) {
	  hash_list[index] = current;

	  /* we must add current to the list of guys that we need */

	  if (*Nexternal == *Nexternal_allocated) {
	    *Nexternal_allocated += (NGhost+10+ *Nexternal_allocated-oldNexternal);
	    temp = (int *) ML_allocate(*Nexternal_allocated*sizeof(int));
	    if (temp == NULL) perror("ML_build_overlapped_pre_comm: Not enough space\n");	  
	    for (k = 0; k < *Nexternal; k++) temp[k] = (*external)[k];
	    ML_free(*external);
	    *external = temp;
	  }
	  (*external)[(*Nexternal)++] = current;
	}
      }
    }
  }
  /* external[] contains ghost nodes from previous overlapping steps as   */
  /* well as the ghost nodes from the current overlapping step. It is     */
  /* only these ones that we wish to sort and consider when building the  */
  /* communication object.                                                */

  ML_az_sort(&((*external)[oldNexternal]),*Nexternal-oldNexternal,NULL,NULL);
  tempA->invec_leng = tempA->outvec_leng;
  ML_CommInfoOP_GenUsingGIDExternals(*Nexternal - oldNexternal,
			     &((*external)[oldNexternal]),max_per_proc,tempA);

  if (val != NULL) ML_free(val);
  if (column != NULL) ML_free(column);

  return 0;
}
/*******************************************************************
 *  Go through the receive list and put the associated global ids
 *  in hash_list[] and extern[]. At the end 'extern' is sorted.
 *  This routine is used for ML_Operator_ProcessorSubdomainOverlap().
 *  We use hash_list[] as a quick and easy way to see if a global id
 *  has already been processed and extern is used to assign local ids
 *  for the ghost nodes.
 *******************************************************************/
int ML_Operator_HashGlobalRcvList(ML_CommInfoOP *pre_comm, int Nrows, 
				  int hash_list[], int hash_length, 
				  int *hash_used, int Gid_assigned_to_proc[], 
				  ML_Comm *comm, 
				  int *Nexternal, int **external,
				  int *Nexternal_allocated) 

{
  double *global_ids;
  int    i, j, k, Nneighbors, Nrcv, *neighbors, NGhost, *rlist, index;
  int    oldNexternal, *temp;

  oldNexternal = *Nexternal;
  Nneighbors   = ML_CommInfoOP_Get_Nneighbors(pre_comm);
  neighbors    = ML_CommInfoOP_Get_neighbors(pre_comm);
  NGhost       = ML_CommInfoOP_Compute_TotalRcvLength(pre_comm);

  /* Estimate space needed for new externals.   */
  /* is allocated later if there is not enough. */
  /* Note: we allow for the possibility that    */
  /* extern[] already has some data that we     */
  /* wish to preserve.                          */

  if ( *Nexternal_allocated - oldNexternal < 2*NGhost) {
    j    = oldNexternal + (NGhost+5)*5;
    temp = (int *) ML_allocate(j*sizeof(int));
    if (temp==NULL) perror("ML_Operator_HashGlobalRcvlist: Out of space\n");
    for (i = 0; i < *Nexternal; i++) temp[i] = (*external)[i];
    if (*Nexternal_allocated != 0) ML_free(*external);
    *external = temp;
    *Nexternal_allocated = j;
  }

  /* Create a vector with global ids for local and ghost nodes */
  /* This is used to associate global ids with the rcv list    */

  global_ids = (double *) ML_allocate(sizeof(double)*(NGhost+Nrows));
  if (global_ids == NULL) perror("ML_Operator_HashGlobalRcvlist: No space\n");
  for (i = 0; i < Nrows; i++) global_ids[i] = Gid_assigned_to_proc[i];
  ML_exchange_bdry(global_ids,pre_comm, Nrows,comm,ML_OVERWRITE,NULL);

  /* Go through Rcv list, if GID associated with receive is not   */
  /* already in the hash table, this means that this is the first */
  /* time that this processor have encountered this GID. We add   */
  /* it to the hash table and put it in the external list.        */

  for (i = 0; i < Nneighbors; i++) {
    Nrcv = ML_CommInfoOP_Get_Nrcvlist(pre_comm, neighbors[i]);
    rlist = ML_CommInfoOP_Get_rcvlist(pre_comm, neighbors[i]);
    for (j = 0; j < Nrcv; j++) {
      ML_hash_it((int)(global_ids[rlist[j]]),hash_list,hash_length,
			 hash_used, &index);
      if (hash_list[index] == -1) {
	if (*Nexternal == *Nexternal_allocated) {/* Need to increase extern[]*/
	  *Nexternal_allocated +=(*Nexternal_allocated+NGhost+10-oldNexternal);
	  temp = (int *) ML_allocate(*Nexternal_allocated*sizeof(int));
	  if (temp == NULL) perror("ML_build_overlapped_pre_comm: No space\n");
	  for (k = 0; k < *Nexternal; k++) temp[k] = (*external)[k];
	  ML_free(*external);
	  *external = temp;
	}
	(*external)[(*Nexternal)++] = (int) global_ids[rlist[j]];
	hash_list[index]            = (int) global_ids[rlist[j]];
      }
    }
    ML_free(rlist);
  }
  ML_free(neighbors);
  ML_free(global_ids);
  ML_az_sort(&((*external)[oldNexternal]),*Nexternal-oldNexternal,NULL,NULL);

  return 0;
}

/*******************************************************************
 *  Take an ML_Operator corresponding to a square matrix distributed
 *  over several processors and make a new ML_Operator corresponding
 *  to an overlapped version of the same matrix. Specifically, each
 *  processor's piece of the original matrix is considered as a subdomain.
 *  For overlap = 1, we do the following on each subdomain
 *      a) make up a set of global indices for the matrix columns
 *      b) append to the matrix any rows corresponding to the 
 *         recv list of the subdomain.
 *      c) Create new send and receive lists to indicate column
 *         indices within the new appended rows that reside on other 
 *         processors.
 *  For overlap > 1, we repeat b) and c) 'overlap' times. When finished,
 *  we post-process the overlapped matrix in the following way:
 *      1) Convert the global indices to local indices. This conversion
 *         is different from the ones in 'ML_back2****' as global column
 *         indices corresponding to overlap regions should use the correct
 *         local index within the processor as opposed to an index on another
 *         processor (that corresponds to this unknown in the 
 *         nonoverlapped matrix). 
 *      2) Throw away any column indices that do not reside on processor.
 *      3) Build a communication object so that we can take unoverlapped 
 *         vectors and map them to overlapped vectors. 
 *
 *
 *******************************************************************/

int ML_overlap(ML_Operator *oldA, ML_Operator *newA, int overlap,
	       ML_CommInfoOP **nonOverlapped_2_Overlapped)
{

  ML_Operator   *tempA, *bogus, *orig_A, *tptr;
  int           N_input_vector, max_per_proc, i, j, k;
  ML_CommInfoOP *getrow_comm = NULL;
  ML_Comm       *comm;
  int           *hash_list = NULL, hash_length, hash_used, Nrows_orig;
  int           *map = NULL, old_Nrows;
  int           nz_ptr, allocated = 0, *column = NULL, row_length;
  double        *val = NULL, *newval = NULL;
  int           current, *newcol = NULL, proc_id, index;
  int           *external = NULL, Nexternal = 0, Nexternal_allocated = 0;
  int           *permute_array = NULL, *rowptr = NULL;
  struct ML_CSR_MSRdata *temp;

  orig_A         = oldA;
  N_input_vector = oldA->invec_leng;
  getrow_comm    = oldA->getrow->pre_comm;
  proc_id        = oldA->comm->ML_mypid;
  comm           = oldA->comm;
  Nrows_orig = oldA->outvec_leng;

  /* Compute global ids for column indices. This makes        */
  /* it easier to keep track of things. Once we are filling   */
  /* the final overlapped matrix, we can replace these global */ 
  /* columns with the correct local column index and build a  */
  /* proper communication object.                             */

  ML_create_unique_col_id(N_input_vector, &(oldA->getrow->loc_glob_map),
			  getrow_comm, &max_per_proc, comm);
  oldA->getrow->use_loc_glob_map = ML_YES;
  map = oldA->getrow->loc_glob_map;

  if ((overlap >= 1) && (getrow_comm != NULL) && (comm->ML_nprocs!=1)) {
    ML_CommInfoOP_Compute_TotalRcvLength(getrow_comm);
    /* estimated space for new ghost variables (will */
    /* resize later if needed).                     */
    Nexternal_allocated = getrow_comm->total_rcv_length;
    Nexternal_allocated += 10;
    Nexternal_allocated *= 2;
    external = (int *) ML_allocate(sizeof(int)*Nexternal_allocated);

    if (max_per_proc == 0 && comm->ML_mypid == 0) {
      printf("WARNING: In ML_overlap, maximum number of local unknowns\n       on any processor (max_per_proc) is zero !\n");
      if (external != NULL) ML_free(external);
      orig_A->getrow->loc_glob_map = NULL;
      orig_A->getrow->use_loc_glob_map = ML_NO;
      ML_free(map);
      return 1;
    }
    hash_length = ML_CommInfoOP_Compute_TotalRcvLength(getrow_comm);
    hash_length = (hash_length + 10)*overlap*overlap;
    hash_list   = (int *) ML_allocate(hash_length*sizeof(int));
    ML_hash_init(hash_list, hash_length, &hash_used);

    /* Put the RcvList in the Hash table and record in the extern list  */
    /* the Gids of the associated Ghost nodes. The hash table will be   */
    /* used later to determine which nonzeros in newly received rows are*/
    /* already contained within the subdomain. The external list is used*/
    /* to build communication structures for the overlapped matrices.   */

    ML_Operator_HashGlobalRcvList(getrow_comm, oldA->outvec_leng,
				  hash_list, hash_length, 
				  &hash_used, map, oldA->comm,
				  &Nexternal, &external,
				  &Nexternal_allocated);

    /* Increase the overlap in oldA one at a time. We do this by appending */
    /* the new rows to oldA, building a new communication object which is  */
    /* needed by ML_exchange_rows(). This object marks all columns that we */
    /* do not yet own and sets up send and receive lists for them.         */

    for (i = 0; i < overlap; i++) {
      old_Nrows = oldA->outvec_leng;
      ML_exchange_rows( oldA, &tempA, getrow_comm);
      oldA = tempA;
      if (i != overlap - 1) {
	ML_build_overlapped_pre_comm(tempA, getrow_comm, max_per_proc,
				     hash_list, hash_length, &hash_used,
				     old_Nrows, &Nexternal, 
				     &external, &Nexternal_allocated);

	getrow_comm = tempA->getrow->pre_comm;
      }
    }
  }
  else {
    tempA = oldA;
    hash_length = 10;
    hash_list = (int *) ML_allocate(hash_length*sizeof(int));
    ML_hash_init(hash_list, hash_length, &hash_used);
  }

  /* clean up the resulting matrix. This means the following:            */
  /*    1) make new copy of the matrix (so sub_matrix is not used)       */
  /*       In the new copy reorder the matrix rows so that all rows      */
  /*       that are received from the same processor are consecutive.    */
  /*       I believe that we need to do this to be consistent with       */
  /*       receive list numbering (in ML and Aztec).                     */
  /*    2) convert global indices to local indices                       */
  /*    3) make a pre_comm structure that can be used to convert         */
  /*       nonoverlapped vectors to overlapped ones.                     */
  /*    4) throw away nonlocal matrix column enties and make an empty    */
  /*       pre_comm structure for the matrix itself.                     */

 if (tempA->N_nonzeros == -1) 
   tempA->N_nonzeros = ML_Operator_ComputeNumNzs(tempA);

 permute_array = (int *) ML_allocate((Nexternal+1)*sizeof(int));
 for (i = 0; i < Nexternal; i++) permute_array[i] = i;
 ML_az_sort(external,Nexternal,permute_array,NULL);
 nz_ptr = 0;
 newcol = (int    *) ML_allocate(sizeof(int   )*(tempA->N_nonzeros+2));
 newval = (double *) ML_allocate(sizeof(double)*(tempA->N_nonzeros+2));
 rowptr = (int    *) ML_allocate(sizeof(int   )*(tempA->outvec_leng+1));
 rowptr[0] = 0;


 for (i = 0; i < tempA->outvec_leng; i++ ) {
    if (i < Nrows_orig) 
      ML_get_matrix_row(tempA, 1, &i, &allocated, &column, &val,
                        &row_length, 0);
    else { /* ghost rows are reordered so that those from the same */
           /* processor are consecutive.                           */
      j = permute_array[i - Nrows_orig] + Nrows_orig;
      ML_get_matrix_row(tempA, 1, &j, &allocated, &column, &val,
                        &row_length, 0);
    }
    for (j = 0; j < row_length; j++) {
      current = column[j];

      /* check if a local guy */
      if ( (current >= max_per_proc*proc_id) &&
	   (current < max_per_proc*(proc_id+1))) {
	newval[nz_ptr  ] = val[j];
	newcol[nz_ptr++] = current - max_per_proc*proc_id;
      }
      else {

	/* If in the hash table, then it is an element that will be */
	/* kept in overlapped matrix. Othewise, this element gets   */
	/* thrown away as it extends too far.                       */
	ML_hash_it( current, hash_list, hash_length, &hash_used, &index);
	if ( hash_list[index] == -1 ) {
	  hash_used--;   /* ML_hash_it assumes that the element will be */
	                 /* added to the hash table and so it increments*/
                         /* hash_used. Since this is not the case, we   */
                         /* must decrement hash_used.                   */
	}
	else {
            k = ML_find_index(current,external,Nexternal);
	    newval[nz_ptr  ] = val[j];
	    newcol[nz_ptr++] = k + Nrows_orig;
	}
      }
    }
    rowptr[i+1] = nz_ptr;
  }

 /* set up new ML operator */

 temp= (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata) );
 temp->columns       = newcol;
 temp->values        = newval;
 temp->rowptr        = rowptr;

 newA->data_destroy = ML_CSR_MSRdata_Destroy;
 ML_Operator_Set_ApplyFuncData(newA, tempA->outvec_leng,
			       tempA->outvec_leng, (void*)temp,
			       tempA->outvec_leng, NULL, 0);
 ML_Operator_Set_Getrow(newA, tempA->outvec_leng, CSR_getrow);
 newA->max_nz_per_row = tempA->max_nz_per_row;
 newA->N_nonzeros     = nz_ptr;
 ML_Operator_Set_ApplyFunc (newA, CSR_matvec);

 /* Create the communication structure necessary to convert nonoverlapped    */
 /* vectors to overlapped vectors. Note:ML_CommInfoOP_GenUsingGIDExternals() */
 /* needs a matrix ... so we create  bogus matrix & set invec_leng (as this  */
 /* is used within the routine).                                             */
 bogus = ML_Operator_Create(tempA->comm);
 bogus->invec_leng = Nrows_orig;
 ML_CommInfoOP_GenUsingGIDExternals(Nexternal, external, max_per_proc, bogus);
 *nonOverlapped_2_Overlapped = bogus->getrow->pre_comm;
 bogus->getrow->pre_comm = NULL;
 ML_Operator_Destroy(&bogus);

 /* The final overlapped matrix is square and requires no communication. */

 newA->getrow->pre_comm = ML_CommInfoOP_Create();

 if (hash_list != NULL) ML_free(hash_list);
 if (external != NULL) ML_free(external);
 if (column   != NULL) ML_free(column);
 if (map      != NULL) ML_free(map);
 orig_A->getrow->loc_glob_map = NULL;
 orig_A->getrow->use_loc_glob_map = ML_NO;
 if (val      != NULL) ML_free(val);
 if (permute_array != NULL) ML_free(permute_array);

 /* Now that we have made a copy of the overlapped matrices, we can clear */
 /* out the memory associated with the intermediate overlapped matrices.  */
 /* The main tricky thing that we must avoid is getting rid of the        */
 /* original nonoverlapped matrix. This matrix appears as a submatrix     */
 /* within tempA .... so we must find this pointer and set it to NULL     */
 /* before we destroy tempA.                                              */

 if (tempA != orig_A) {
   tptr = tempA;
   while ( (tptr!= NULL) && (tptr->sub_matrix != orig_A))
     tptr = tptr->sub_matrix;
   if (tptr != NULL) tptr->sub_matrix = NULL;
   ML_RECUR_CSR_MSRdata_Destroy(tempA);
   ML_Operator_Destroy(&tempA);
  }
  return 0;

}    
 
/*******************************************************************************
  ML_Operator_ReportStatistics() dumps out static, communication, and
  performance statistics for a particular operator.  It's meant to be used in
  conjunction with ML_Operator_Profile().
*******************************************************************************/
#define ML_FUNCTION_NAME "ML_Operator_ReportStatistics"
void ML_Operator_ReportStatistics(ML_Operator *mat, char *appendlabel,
                                  int perfAndCommStats)
{
  double t1;
  double j, Nglobrows, Nglobcols;
  int i,k,NumActiveProc, proc_active;
  ML_Comm *comm = mat->comm;
  int mypid = mat->comm->ML_mypid;
  int maxrows,maxproc,minrows,minproc;
  int total_rcv_leng = 0;
  char *origlabel = NULL, *modlabel = NULL;
  ML_CommInfoOP *c_info;
  int minnzs,maxnzs;
  char eqLine[76];

  if (ML_Get_PrintLevel() == 0)
    return;
  modlabel = (char *) ML_allocate(80 * sizeof(char));
  origlabel = mat->label;
  if (mat->label == NULL) {
    if (appendlabel != NULL)
      sprintf(modlabel,"unlabeled_operator_%s",appendlabel);
    else
      sprintf(modlabel,"unlabeled_operator");
  }
  else {
    if (appendlabel != NULL)
      sprintf(modlabel,"%s_%s",mat->label,appendlabel);
    else
      sprintf(modlabel,"%s",mat->label);
  }

  mat->label = modlabel;

  if (mat->invec_leng > 0 || mat->outvec_leng > 0)
    proc_active = 1;
  else proc_active = 0;
  NumActiveProc = ML_gsum_int(proc_active, comm);

  i = mat->invec_leng; 
  Nglobcols = ML_gsum_double((double)i, comm);
  i = mat->outvec_leng; 
  Nglobrows = ML_gsum_double((double)i, comm);

  for (i=0; i<75; i++) eqLine[i] = '=';
  eqLine[75] = '\0';

  if  (NumActiveProc > 0)
  {

    /* static operator statistics */
    if (mat->getrow->pre_comm != NULL && proc_active)
      i = mat->getrow->pre_comm->N_neighbors;
    else i = 0;
    i = ML_Comm_GsumInt(comm, i);
    j = ML_Comm_GsumDouble(comm, (double) mat->N_nonzeros);
    maxnzs = ML_gmax_int(mat->N_nonzeros, comm);
    maxproc = ML_gmax_int( (maxnzs == mat->N_nonzeros ? mypid:0), comm);
    if (proc_active) k=mat->N_nonzeros;
    else k = maxnzs;
    minnzs = ML_gmin_int(k, comm);
    minproc = ML_gmax_int((minnzs == mat->N_nonzeros ? mypid:0), comm);
    t1 = ML_gsum_double( (proc_active ? (double) mat->N_nonzeros : 0.0), comm);
    t1 = t1/((double) NumActiveProc);
    if (mypid == 0) {
       printf("= %s =\n",eqLine);
       printf("%s: %1.0f rows, %1.0f cols, %1.0f global nonzeros\n",
              mat->label,Nglobrows,Nglobcols, j);
       printf("%s: num PDES = %d, lambda_max = %2.3e, lambda_min = %2.3e\n",
              mat->label, mat->num_PDEs, mat->lambda_max, mat->lambda_min);
       printf("%s: %2.3e avg nbrs, %d active proc\n\n",
              mat->label,((double) i)/((double) NumActiveProc), NumActiveProc);
       printf("%s: max nonzeros (pid %d) \t= %d\n",
              mat->label,maxproc,maxnzs);
       printf("%s: min nonzeros (pid %d) \t= %d\n",
              mat->label,minproc,minnzs);
       printf("%s: avg nonzeros \t\t= %2.3e\n",
              mat->label,t1);
    }

    maxrows = ML_gmax_int( mat->outvec_leng , comm);
    maxproc = ML_gmax_int( (maxrows == mat->outvec_leng ? mypid:0), comm);
    minrows = ML_gmin_int( (mat->outvec_leng > 0 ? mat->outvec_leng: maxrows),
                             comm);
    minproc = ML_gmax_int((minrows == mat->outvec_leng ? mypid:0), comm);
    if (mypid == 0) {
      printf("%s: max number of rows (pid %d) \t= %d\n",
             mat->label,maxproc,maxrows);
      printf("%s: min number of rows (pid %d) \t= %d\n",
             mat->label,minproc,minrows);
    }
    t1 = ML_gsum_double( (proc_active ? (double) mat->outvec_leng : 0.0), comm);
    t1 = t1/((double) NumActiveProc);
    if (mypid == 0)
       printf("%s: avg number of rows \t\t= %2.3e\n",mat->label,t1);

    if (perfAndCommStats == ML_TRUE) {
      /* communication statistics */
      if (mypid == 0) printf("\n");
      if (mat->getrow->pre_comm != NULL) {
        c_info = mat->getrow->pre_comm;
        total_rcv_leng = ML_CommInfoOP_Compute_TotalRcvLength(c_info);
        maxrows = ML_gmax_int( total_rcv_leng , comm);
        maxproc = ML_gmax_int( (maxrows == total_rcv_leng ? mypid:-1), comm);
        minrows = ML_gmin_int( (total_rcv_leng > 0 ? total_rcv_leng: maxrows),
                               comm);
        minproc = ML_gmax_int((minrows == total_rcv_leng ? mypid:-1), comm);
        if (mypid == 0) {
          printf("%s: max pre_comm recv length (pid %d) \t= %d\n",
                 mat->label,maxproc,maxrows);
          printf("%s: min pre_comm recv length (pid %d) \t= %d\n",
                 mat->label,minproc,minrows);
        }
        t1 = ML_gsum_double( (proc_active ? (double) total_rcv_leng :0.0),comm);
        t1 = t1/((double) NumActiveProc);
        if (mypid == 0)
          printf("%s: avg pre_comm recv length \t\t= %2.3e\n",mat->label,t1);
        if (mypid == 0) printf("\n");
  
        t1 = ML_gmax_double( (proc_active ? c_info->time: 0.0), comm);
        i =ML_gmax_int((t1 == c_info->time ? mypid:0),comm);
        if (mypid == 0)
           printf("%s: max pre_comm exch bdry time (pid %d) \t= %2.3e\n",
                  mat->label,i,t1);
        t1 = - c_info->time;
        t1 = ML_gmax_double( (proc_active ? t1: -1.0e20), comm);
        t1 = - t1;
        i =ML_gmax_int((t1 == c_info->time ? mypid:0),comm);
        if (mypid == 0)
           printf("%s: min pre_comm exch bdry time (pid %d) \t= %2.3e\n",
                  mat->label,i,t1);
        t1 = ML_gsum_double( (proc_active ? c_info->time : 0.0), comm);
        t1 = t1/((double) NumActiveProc);
        if (mypid == 0)
           printf("%s: avg pre_comm exch bdry time \t\t= %2.3e\n\n",
                  mat->label,t1);
      }
  
      if (mat->getrow->post_comm != NULL) {
        c_info = mat->getrow->post_comm;
        total_rcv_leng = ML_CommInfoOP_Compute_TotalRcvLength(c_info);
        maxrows = ML_gmax_int( total_rcv_leng , comm);
        maxproc = ML_gmax_int( (maxrows == total_rcv_leng ? mypid:-1), comm);
        minrows = ML_gmin_int( (total_rcv_leng > 0 ? total_rcv_leng: maxrows),
                               comm);
        minproc = ML_gmax_int((minrows == total_rcv_leng ? mypid:-1), comm);
        if (mypid == 0) {
          printf("%s: max post_comm recv length (pid %d) \t= %d\n",
                 mat->label,maxproc,maxrows);
          printf("%s: min post_comm recv length (pid %d) \t= %d\n",
                 mat->label,minproc,minrows);
        }
        t1 = ML_gsum_double( (proc_active ? (double) total_rcv_leng : 0.0), comm);
        t1 = t1/((double) NumActiveProc);
        if (mypid == 0)
          printf("%s: avg post_comm recv length \t= %2.3e\n",mat->label,t1);
        if (mypid == 0) printf("\n");
        t1 = ML_gsum_double( (proc_active ? c_info->time : 0.0), comm);
        t1 = t1/((double) NumActiveProc);
        if (mypid == 0)
           printf("%s: avg post_comm exch boundary time \t= %2.3e\n",mat->label,t1);
        t1 = ML_gmax_double( (proc_active ? c_info->time: 0.0), comm);
        i =ML_gmax_int((t1 == c_info->time ? mypid:0),comm);
        if (mypid == 0)
           printf("%s: max post_comm exch bdry time (pid %d) \t= %2.3e\n",mat->label,i,t1);
        t1 = - c_info->time;
        t1 = ML_gmax_double( (proc_active ? t1: -1.0e20), comm);
        t1 = - t1;
        i =ML_gmax_int((t1 == c_info->time ? mypid:0),comm);
        if (mypid == 0)
           printf("%s: min post_comm exch bdry time (pid %d) \t= %2.3e\n",mat->label,i,t1);
        if (mypid == 0) printf("\n");
      }

      /* performance statistics */
      if (mypid == 0)
         printf("%s: number of applies \t= %d\n",mat->label,mat->ntimes);
      t1 = ML_gmax_double( (proc_active ? mat->apply_time : 0.0 ), mat->comm);
      i =ML_gmax_int((t1 == mat->apply_time ? mypid:0),mat->comm);
      if (mypid == 0)
         printf("%s: max apply+comm time (pid %d) \t= %2.3e\n",mat->label,i,t1);
      t1 = ML_gsum_double( (mypid==i ? mat->apply_without_comm_time:0.0),
                           mat->comm);
      if (mypid == 0)
         printf("%s:     apply only time (pid %d) \t= %2.3e\n",mat->label,i,t1);
      t1 = - mat->apply_time;
      t1 = ML_gmax_double( (proc_active ? t1: -1.0e20), mat->comm);
      t1 = - t1;
      i =ML_gmax_int((t1 == mat->apply_time ? mypid:0),mat->comm);
      if (mypid == 0)
         printf("%s: min apply+comm time (pid %d) \t= %2.3e\n",mat->label,i,t1);
      t1 = ML_gsum_double( (mypid==i ? mat->apply_without_comm_time:0.0),
                           mat->comm);
      if (mypid == 0)
         printf("%s:     apply only time (pid %d) \t= %2.3e\n",mat->label,i,t1);
      t1 = ML_gsum_double( (proc_active ? mat->apply_time : 0.0), mat->comm);
      t1 = t1/((double) NumActiveProc);
      if (mypid == 0)
         printf("%s: avg apply+comm time \t\t= %2.3e\n\n",mat->label,t1);
  
      t1 = ML_gmax_double( (proc_active ?mat->apply_without_comm_time:0.0),mat->comm);
      i =ML_gmax_int((t1 == mat->apply_without_comm_time ? mypid:0),mat->comm);
      if (mypid == 0)
         printf("%s: max apply only time (pid %d) \t= %2.3e\n",mat->label,i,t1);
      t1 = - mat->apply_without_comm_time;
      t1 = ML_gmax_double( (proc_active ? t1: -1.0e20), mat->comm);
      t1 = - t1;
      i =ML_gmax_int((t1 == mat->apply_without_comm_time ? mypid:0),mat->comm);
      if (mypid == 0)
         printf("%s: min apply only time (pid %d) \t= %2.3e\n",mat->label,i,t1);

      t1 = ML_gsum_double( (proc_active ? mat->apply_without_comm_time:0.0),
                            mat->comm);
      t1 = t1/((double) NumActiveProc);
      if (mypid == 0)
         printf("%s: avg apply only time \t\t= %2.3e\n",mat->label,t1);

    } /* if perfAndCommStats == ML_TRUE*/
    if (mypid == 0) {
       printf("= %s =\n",eqLine);
       fflush(stdout);
    }
  }

  mat->label = origlabel;
  ML_free(modlabel);
}
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

/*******************************************************************************
  ML_Operator_Profile() profiles a particular ML_Operator.  The apply time,
  communication time, and number of applies are saved.  All processor are
  sync'd up, then the operator is applied to a vector N times.  The results of
  this are reported by ML_Operator_ReportStatistics().  Finally, the old apply
  time, communication time, and number of applies are set for the operator.

  ML_Operator *A            operator to profile
  char *appendlabel         string to append to operator's label for id purposes
  int numits                number of applies to perform
*******************************************************************************/

void ML_Operator_Profile(ML_Operator *A, char *appendlabel)
{
  int numits = ML_Operator_Profile_GetIterations();
#if defined(ML_TIMING)
  int j, ntimes;
  double *xvec,*bvec;
  double apply_time, apply_without_comm_time, pre_time=0.0, post_time=0.0;
  double t0;
  double totalPreTime=0.0;
  double totalPostTime=0.0;

  if (ML_Get_PrintLevel() == 0) return;

  if (numits == 0) {
    ML_Operator_ReportStatistics(A,appendlabel,ML_FALSE);
    return;
  }

  xvec = (double *) ML_allocate((A->invec_leng) * sizeof(double));
  ML_random_vec(xvec, A->invec_leng, A->comm);
  bvec = (double *) ML_allocate((A->outvec_leng) * sizeof(double));

  apply_time = A->apply_time; A->apply_time = 0.0;
  apply_without_comm_time = A->apply_without_comm_time; A->apply_without_comm_time = 0.0;
  ntimes = A->ntimes; A->ntimes = 0;
  if (A->getrow->pre_comm != NULL) {
    pre_time = A->getrow->pre_comm->time;
    A->getrow->pre_comm->time = 0.0;
  }
  if (A->getrow->post_comm != NULL) {
    post_time = A->getrow->post_comm->time;
    A->getrow->post_comm->time = 0.0;
  }

  if (numits < 0) {
    /* Figure out how many iterations we need to do get apply+comm time
       equal to abs(numits) seconds, but don't let the user shoot herself
       in the foot. */
    int proc_active, NumActiveProc;
    numits = -numits;
    if (numits > 30) numits = 30;
    ML_Comm_Barrier(A->comm);
    ML_Operator_Apply(A,A->invec_leng,xvec,A->outvec_leng,bvec);
    if (A->invec_leng > 0 || A->outvec_leng > 0) proc_active = 1;
    else                                         proc_active = 0;
    NumActiveProc = ML_gsum_int(proc_active, A->comm);
    if (NumActiveProc > 0) {
      /* this could be skewed if max >> avg */
      t0 = ML_gsum_double( (proc_active ? A->apply_time : 0.0), A->comm);
      t0 = t0/NumActiveProc;
      numits = (int) (1.1 * floor(numits / t0));
    }
    else numits = 0;
    A->apply_time = 0.0;
    A->apply_without_comm_time = 0.0;
    A->ntimes = 0;
    if (A->getrow->pre_comm != NULL) A->getrow->pre_comm->time = 0.0;
    if (A->getrow->post_comm != NULL) A->getrow->post_comm->time = 0.0;
  }

  ML_Comm_Barrier(A->comm);
  t0 = GetClock();
  for (j=0; j<numits; j++) {
    ML_Operator_Apply(A,A->invec_leng,xvec,A->outvec_leng,bvec);
  ML_Comm_Barrier(A->comm);
  }
  /* If communication time is zero, the matrix was probably created by Epetra */
  /* (i.e. outside ML). Thus, ML cannot time only communication. However, ML  */
  /* always creates its own communication widget via ML_CommInfoOP_Generate() */
  /* which can be used to estimate communication time.                        */

  if (A->getrow->pre_comm != NULL)
    totalPreTime = ML_gsum_double(A->getrow->pre_comm->time,A->comm);
  if (A->getrow->post_comm != NULL)
    totalPostTime = ML_gsum_double(A->getrow->post_comm->time,A->comm);

  if ( (A->comm->ML_nprocs > 1) && (A->getrow->pre_comm != NULL) &&
       (totalPreTime == 0.0) &&
       ((A->getrow->post_comm == NULL) || (totalPostTime == 0.0))){

    if ( (ML_Get_PrintLevel() != 0) && (A->comm->ML_mypid == 0))
       printf("==> %s communication time that follows is only an estimate!!\n",
              A->label);

    ML_free(bvec);
    bvec = (double *) ML_allocate((A->getrow->pre_comm->minimum_vec_size+
                                   A->invec_leng+1)*sizeof(double));
    if (bvec == NULL)
       pr_error("ML_Operator_Profile(%d): out of space\n",A->comm->ML_mypid);

    for (j = 0; j < A->invec_leng; j++) bvec[j] = xvec[j];

    for (j=0; j<numits; j++) {
       ML_exchange_bdry(bvec,A->getrow->pre_comm, A->invec_leng, 
                        A->comm, ML_OVERWRITE,NULL);
    ML_Comm_Barrier(A->comm);
    }

    A->apply_without_comm_time = A->apply_time - A->getrow->pre_comm->time;
  }

  ML_Operator_ReportStatistics(A,appendlabel,ML_TRUE);

  A->apply_time = apply_time;
  A->apply_without_comm_time = apply_without_comm_time;
  A->ntimes = ntimes;
  if (A->getrow->pre_comm != NULL)
    A->getrow->pre_comm->time = pre_time;
  if (A->getrow->post_comm != NULL)
    A->getrow->post_comm->time = post_time;

  ML_free(xvec);
  ML_free(bvec);
#else
  if (ML_Get_PrintLevel() == 0) return;
  if (A->comm->ML_mypid == 0 && numits > 0)
    printf("ML_Operator_Profile: not compiled with -DML_TIMING\n");
  ML_Operator_ReportStatistics(A,appendlabel,ML_FALSE);
#endif
}

int ML_profile_num_its=0;
void ML_Operator_Profile_SetIterations(int numits)
{
  ML_profile_num_its = numits;
}

int ML_Operator_Profile_GetIterations()
{
  return ML_profile_num_its;
}

int ML_Operator_Get_Nnz(ML_Operator *A)
{
  int i;
  int space=0, *columns=NULL, row_lengths;
  double *values = NULL;


  if (A == NULL) return 0;
  if (A->getrow == NULL) return 0;
  if (A->getrow->func_ptr == NULL) return 0;

  if (A->N_nonzeros == -1) {
    A->N_nonzeros = 0;
    for (i=0; i<A->outvec_leng; i++) {
      ML_get_matrix_row(A, 1, &i, &space, &columns, &values, &row_lengths, 0);
      A->N_nonzeros += row_lengths;
    }
    if (columns != NULL) ML_free(columns);
    if (values != NULL) ML_free(values);
  }
  return A->N_nonzeros; 
}

/* ******************************************************************** */
/* Functions to color matrices to cheapen ml_agg_min_energy.c           */
/* ******************************************************************** */
/* Author        : Ray Tuminaro (SNL)                                   */
/* Date          : Sept, 2005                                           */
/* ******************************************************************** */

int ML_Operator_MisRootPts( ML_Operator *Amatrix,  int num_PDE_eqns,
			    int **aggr_index)
{
  /*
   * Take Amatrix and perform a maximal independent set on the associated
   * block matrix (where the block size is num_PDE_eqns). The resulting
   * independent set is returned in aggr_index[] where 
   *            aggr_index[k] = -1  ==> not part of independent set
   *            aggr_index[k] > -1  ==> aggr_index[k]th local point in 
   *                                    independent set
   * NOTE: This routine has not yet been checked in parallel.
   *
   */

  int Nrows, index, i, j, kk, aggr_count;
  int *rcvleng= NULL, *sndleng= NULL;
  char                  *vtype, *state, *bdry;
   int                   nvertices, *vlist, ntotal, Nneigh;
   int                   *neigh = NULL, max_element, Nghost;
   int                   **sndbuf = NULL, **rcvbuf = NULL;
   int                   **sendlist = NULL, **recvlist = NULL;
   ML_CommInfoOP         *mat_comm;
   int                   allocated = 0, *rowi_col = NULL, rowi_N;
   double                *rowi_val = NULL, *dtemp;
   int                   *templist, **proclist, *rowptr, *columns;
   ML_Comm               *comm;

   Nrows        = Amatrix->outvec_leng;
   comm         = Amatrix->comm;

   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(1);
   }

   /*= `amalgamate' refers to the conversion from block matrices
     to point matrices */
   
   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, 0.0);

   Nrows     /= num_PDE_eqns;
   nvertices  = Amatrix->outvec_leng;

   mat_comm  = Amatrix->getrow->pre_comm;
   Nneigh    = ML_CommInfoOP_Get_Nneighbors(mat_comm);
   neigh     = ML_CommInfoOP_Get_neighbors(mat_comm);
   sendlist  = (int **) ML_allocate(sizeof(int *)*Nneigh);
   recvlist  = (int **) ML_allocate(sizeof(int *)*Nneigh);
   rcvbuf    = (int **) ML_allocate(sizeof(int *)*Nneigh);
   sndbuf    = (int **) ML_allocate(sizeof(int *)*Nneigh);
   rcvleng   = (int  *) ML_allocate(sizeof(int  )*Nneigh);
   sndleng   = (int  *) ML_allocate(sizeof(int  )*Nneigh);

   max_element = nvertices - 1;
   for (i = 0; i < Nneigh; i++) {
      recvlist[i]  = ML_CommInfoOP_Get_rcvlist(mat_comm, neigh[i]);
      rcvleng[i]   = ML_CommInfoOP_Get_Nrcvlist (mat_comm, neigh[i]);
      sendlist[i]  = ML_CommInfoOP_Get_sendlist (mat_comm, neigh[i]);
      sndleng[i]   = ML_CommInfoOP_Get_Nsendlist(mat_comm, neigh[i]);
      rcvbuf[i]    = (int *) ML_allocate(sizeof(int)*(rcvleng[i]+1));
      sndbuf[i]    = (int *) ML_allocate(sizeof(int)*(sndleng[i]+1));
                           /* +1 needed inside ML_Aggregate_LabelVertices */
      for (j = 0; j < rcvleng[i]; j++) 
         if (recvlist[i][j] > max_element ) 
            max_element = recvlist[i][j];
   }
   Nghost = max_element - nvertices + 1;
   ntotal = nvertices + Nghost;

   templist = (int *) ML_allocate(sizeof(int)*nvertices);
   for ( i = 0; i < nvertices; i++ ) templist[i] = 0;
   for ( i = 0; i < Nneigh; i++ ) {
      for ( j = 0; j < sndleng[i]; j++ ) {
         index = sendlist[i][j];
         if ( index >= nvertices || index < 0 ) {
	   printf("%d : Error : in sendlist.\n", comm->ML_mypid);
            exit(0);
         }
         templist[index]++;
	 /*= templist[j] is the number of processors who need `j' */	 
      }
   }

   /* Allocate proclist to record the processors and indices each of */
   /* my local vertices are to send.  The first element of the array */
   /* is a counter of how many processors, followed by a number of   */
   /* processor and index pairs.                                     */

   proclist = (int **) ML_allocate(ntotal * sizeof( int *));
   for ( i = 0; i < nvertices; i++ ) {
      proclist[i]    =(int *) ML_allocate( (2*templist[i]+1) * sizeof( int ) );
      proclist[i][0] = 0;
      templist[i]    = 0;
   }
   for ( i = 0; i < Nneigh; i++ ) {
      for ( j = 0; j < sndleng[i]; j++ ) {
         index = sendlist[i][j];
         proclist[index][templist[index]+1] = i;
         proclist[index][templist[index]+2] = j;
         templist[index] += 2;
         proclist[index][0]++;
      }
   }
   for ( i = nvertices; i < ntotal; i++ ) {
      proclist[i] = (int *) ML_allocate( sizeof( int ) );
   }
   for ( i = 0; i < Nneigh; i++ ) {
      for ( j = 0; j < rcvleng[i]; j++ ) {
         index = recvlist[i][j];
         proclist[index][0] = neigh[i];
      }
   }
   ML_free(templist);

   /* record the Dirichlet boundary and count number of nonzeros */

   bdry = (char *) ML_allocate(sizeof(char)*(ntotal + 1));
   for (i = Nrows ; i < ntotal; i++) bdry[i] = 'F';
   j = 0;
   for (i = 0; i < Nrows; i++) {
      bdry[i] = 'T';
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      if (rowi_N > 1) bdry[i] = 'F';
      j += rowi_N;
   }

   /* communicate the boundary information */

   dtemp = (double *) ML_allocate(sizeof(double)*(ntotal+1));
   for (i = nvertices; i < ntotal; i++) dtemp[i] = 0;
   for (i = 0; i < nvertices; i++) {
      if (bdry[i] == 'T') dtemp[i] = 1.;
      else  dtemp[i] = 0.;
   }
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,nvertices,comm,
                    ML_OVERWRITE,NULL);
   for (i = nvertices; i < ntotal; i++) {
      if (dtemp[i] == 1.) bdry[i] = 'T';
      else bdry[i] = 'F';
   }
   ML_free(dtemp);

   /* make a csr version of nonzero pattern */
   rowptr  = (int *) ML_allocate(sizeof(int)*(nvertices+1));
   columns = (int *) ML_allocate(sizeof(int)*(j+1));
   j = 0;
   rowptr[0] = 0;
   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      for (kk = 0; kk < rowi_N; kk++) columns[j++] = rowi_col[kk];
      rowptr[i+1] = j;
   }
   if (rowi_col != NULL) ML_free(rowi_col);
   if (rowi_val != NULL) ML_free(rowi_val);

   (*aggr_index) = (int *) ML_allocate(sizeof(int)* ntotal*num_PDE_eqns);
   for (i = 0; i < ntotal; i++) (*aggr_index)[i] = -1;

   vlist  = (int *) ML_allocate(sizeof(int)* nvertices);
   state  = (char *) ML_allocate(sizeof(char)* ntotal);
   vtype  = (char *) ML_allocate(sizeof(char)* ntotal);
   for (i = 0; i < nvertices; i++)  vlist[i] =  i;
   for (i = 0; i < ntotal;    i++)  state[i] = 'F';
   for (i = 0; i < ntotal;    i++)  vtype[i] = 'x';

   /* delete nodes that are just isolated Dirichlet points */

   for (i = 0; i < nvertices ; i++) {
     if (bdry[i] == 'T') state[i] = 'D'; 
   }

   aggr_count = ML_Aggregate_LabelVertices(nvertices, vlist,'x',state, vtype, 
					   nvertices, rowptr, 
					   columns, comm->ML_mypid,
					   proclist, 
                      Nneigh,sndbuf,neigh, sndleng,
                      Nneigh,rcvbuf, neigh, rcvleng,
                      recvlist, 1532, comm, *aggr_index);

   /* free memory used for doing the MIS stuff */

   for ( i = 0; i < ntotal; i++ ) ML_free(proclist[i]);
   ML_free(proclist);
   ML_free(vlist); ML_free(state); ML_free(vtype);
   for (i = 0; i < Nneigh; i++) {
      ML_free(recvlist[i]);
      ML_free(sendlist[i]);
      ML_free(rcvbuf[i]);
      ML_free(sndbuf[i]);
   }
   ML_free(sndleng); ML_free(rcvleng);  ML_free(sndbuf);
   ML_free(rcvbuf);  ML_free(recvlist); ML_free(sendlist);
   ML_free(neigh);


   ML_free(bdry); 
   ML_free(columns); ML_free(rowptr);

   /* communicate aggregate information so that ghost nodes which */
   /* are root points are identified. Note: local ids ?associated */

   dtemp = (double *) ML_allocate(sizeof(double)*(ntotal+1));
   for (i = 0; i < aggr_count; i++) dtemp[i] = (double) (*aggr_index)[i];
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,nvertices,comm,
                    ML_OVERWRITE,NULL);
   for (i = nvertices; i < ntotal; i++) {
      if (dtemp[i] != -1.) (*aggr_index)[i] = (int) dtemp[i];
   }
   ML_free(dtemp);
   ML_Operator_UnAmalgamateAndDropWeak(Amatrix, num_PDE_eqns, 0.0);

   for (i = ntotal - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
	if ( (*aggr_index)[i] == -1) 
	  (*aggr_index)[i*num_PDE_eqns+j] = -1;
	else 
	  (*aggr_index)[i*num_PDE_eqns+j] = num_PDE_eqns*(*aggr_index)[i] + j;
      }
   }

   /* Renumber the aggregate indices so that they are consecutive.  */
   /* This will make it easier to determine ghost/communication     */
   /* information when these MIS points are used to compress columns*/
   /* out of a matrix.                                              */
   j = 0;
   for (i = 0; i < ntotal*num_PDE_eqns; i++) {
     if ( (*aggr_index)[i] != -1.) {
       (*aggr_index)[i] = j++;
     }
   }

   return(aggr_count*num_PDE_eqns);
}

/*************************************************************************
Modify matrix by dropping small entries.
 *************************************************************************/
int ML_CSR_DropSmall(ML_Operator *Pe, double AbsoluteDrop,
		     double RelativeRowDrop, double RelativeColDrop)

{
  /* Drop entry Pe(I,J) when all of the following conditions are true:
   *
   *     1.  | Pe(I,J) | < AbsoluteDrop
   *     2.  | Pe(I,J) | < RelativeRowDrop * max | Pe(I,:) |
   *     3.  | Pe(I,J) | < RelativeColDrop * max | Pe(:,J) |
   *                Note: if RelativeColDrop != 0, we do not
   *                      drop columns associated with ghost nodes.
   */

  struct ML_CSR_MSRdata *csr_data;
  int            lower, nz_ptr, i, j, nn;
  double         row_max = 0., dtemp;
  double         *col_maxs = NULL;

  if ( (Pe->getrow == NULL) ||  (Pe->getrow->func_ptr != CSR_getrow)) {
     printf("ML_CSR_DropSmall can only be used with CSR matrices\n");
     return -1;
  }
  csr_data = (struct ML_CSR_MSRdata *) Pe->data;

  AbsoluteDrop    = ML_abs(AbsoluteDrop);
  RelativeRowDrop = ML_abs(RelativeRowDrop);
  RelativeColDrop = ML_abs(RelativeColDrop);

  if (RelativeColDrop != 0.0) {
    nn = Pe->invec_leng + ML_CommInfoOP_Compute_TotalRcvLength(
				       Pe->getrow->pre_comm);
    col_maxs = (double *) ML_allocate(sizeof(double)*(nn+1));
    for (i = 0; i < nn; i++) col_maxs[i] = 0.;
    for (i = 0; i < Pe->outvec_leng; i++) {
      for (j = csr_data->rowptr[i]; j < csr_data->rowptr[i+1]; j++) {
	nn = csr_data->columns[j];
	if (ML_abs(csr_data->values[j]) > col_maxs[nn] )
	  col_maxs[nn] = ML_abs(csr_data->values[j]);
      }
    }
    for (i = 0; i < nn; i++) col_maxs[nn] *= RelativeColDrop; 
    for (i = Pe->invec_leng ; i < nn; i++) col_maxs[nn] = 0.;
  }

  lower = csr_data->rowptr[0];
  nz_ptr = 0;
  for (i = 0; i < Pe->outvec_leng; i++) {
    if (RelativeRowDrop != 0.0) {
      row_max = 0.;
      for (j = lower; j < csr_data->rowptr[i+1]; j++) {
	if (ML_abs(csr_data->values[j]) > row_max) 
	  row_max = ML_abs(csr_data->values[j]);
      }
      if (row_max > 1.) row_max = 1.;
                  /* Ugh: For now ml_agg_min_energy seems   */
                  /* to work better with the above line in. */
      row_max *= RelativeRowDrop;
      if (AbsoluteDrop < row_max) row_max = AbsoluteDrop;
    }
    else row_max = AbsoluteDrop;

    if (RelativeColDrop == 0.) {
      for (j = lower; j < csr_data->rowptr[i+1]; j++) {
	dtemp = ML_abs(csr_data->values[j]);
	if  (dtemp > row_max) {
	  csr_data->values[nz_ptr] = csr_data->values[j];
	  csr_data->columns[nz_ptr] = csr_data->columns[j];
	  nz_ptr++;
	}
      }
    }
    else {
      for (j = lower; j < csr_data->rowptr[i+1]; j++) {
	dtemp = ML_abs(csr_data->values[j]);
	if  ( (dtemp > row_max) || (dtemp > col_maxs[csr_data->columns[j]])) {
	  csr_data->values[nz_ptr] = csr_data->values[j];
	  csr_data->columns[nz_ptr] = csr_data->columns[j];
	  nz_ptr++;
	}
      }
    }
    lower = csr_data->rowptr[i+1];
    csr_data->rowptr[i+1] = nz_ptr;
  }
  Pe->N_nonzeros = nz_ptr;

  if (col_maxs != NULL) ML_free(col_maxs);
  return 0;

}
/* Take a CSR matrix and remove all columns which are labelled as */
/* -1 in the array 'subset' */

ML_Operator *ML_CSRmatrix_ColumnSubset(ML_Operator *Amat, int Nsubset,
					       int subset[])
{
  struct ML_CSR_MSRdata *csr_data;
  int    *row_ptr, *columns, *newrowptr, *newcolumns;
  double *values, *newvalues;
  int    i, j, n, count, max_nz_per_row, nnz_in_row;
  ML_Operator *Amat_subset;

  /* This only works for CSR matrices */
  if (Amat->getrow->func_ptr != CSR_getrow) return NULL;

  n = Amat->outvec_leng;

  csr_data = (struct ML_CSR_MSRdata *) Amat->data;
  row_ptr = csr_data->rowptr;
  columns = csr_data->columns;
  values  = csr_data->values;

  /* Count the number of nonzeros in the new compressed matrix */

  count = 0;
  for (i = 0; i < n ; i++) {
    for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
      if ( subset[columns[j]] != -1) count++;
    }
  }
  newrowptr  = (int    *)  ML_allocate(sizeof(int)*(n+1));
  newcolumns = (int    *)  ML_allocate(sizeof(int)*(count+1));
  newvalues  = (double *)  ML_allocate(sizeof(double)*(count+1));
  count      = 0;
  newrowptr[0]   = 0;
  max_nz_per_row = 0;
  for (i = 0; i < n ; i++) {
    nnz_in_row = 0;
    for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
      if ( subset[columns[j]] != -1) {
	nnz_in_row++;
	newvalues[count   ] = values[j];
	newcolumns[count++] = subset[columns[j]];
      }
    }
    if (nnz_in_row > max_nz_per_row) max_nz_per_row = nnz_in_row;
    newrowptr[i+1] = count;
  }
  csr_data =(struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  if (csr_data == NULL) pr_error("no space for csr_data\n");
  csr_data->columns = newcolumns;
  csr_data->values  = newvalues;
  csr_data->rowptr  = newrowptr;
  Amat_subset = ML_Operator_Create(Amat->comm);

  ML_Operator_Set_ApplyFuncData(Amat_subset, Nsubset, n, csr_data, n, NULL, 0);
  ML_Operator_Set_Getrow(Amat_subset, n, CSR_getrow);
  ML_Operator_Set_ApplyFunc (Amat_subset, CSR_matvec);
  Amat_subset->getrow->pre_comm = ML_CommInfoOP_SqueezeColumns(
               Amat->getrow->pre_comm, Nsubset, subset);
  Amat_subset->data_destroy   = ML_CSR_MSRdata_Destroy;
  Amat_subset->max_nz_per_row = max_nz_per_row;
  Amat_subset->N_nonzeros     = count;
return Amat_subset;
}

/*
  Function: ML_Operator_IdentifyDirichletRows()
  Returns a character array that identifies local Dirichlet rows.  For local
  row i,
     A->DirichletRows[i] = 'T' if row i has no nonzero off-diagonals
     A->DirichletRows[i] = 'F' if row i has at least one nonzero off-diagonal.

  input:   pointer to ML_Operator

  comments:  The memory referenced by the returned pointer is owned
             by the ML_Operator and should not be deleted.

             The first time through, this function will populate the field
             A->DirichletRows.  Any subsequent call will return the address
             contained in this field.
*/

char* ML_Operator_IdentifyDirichletRows(ML_Operator *A)
{
  if (A==NULL) return NULL;
  if (A->ML_id != ML_ID_OP) {
    pr_error("ML_Operator_IdentifyDirichletRows: not an ML_Operator.\n");
  }

  if (A->DirichletRows == NULL) {
    int    i,j,allocated=50,row_lengths;
    int    *cols;
    double *vals;

    A->DirichletRows = (char *) ML_allocate(A->outvec_leng * sizeof(char));
    cols = (int *) ML_allocate(allocated*sizeof(int));
    vals = (double *) ML_allocate(allocated*sizeof(double));
    for (i=0; i<A->outvec_leng; i++) {
      ML_get_matrix_row(A,1,&i,&allocated,&cols,&vals,&row_lengths,0);
      A->DirichletRows[i] = 'T';
      for (j=0; j<row_lengths; j++)
        if (cols[j] != i  &&  vals[j] != 0.0) {  /*only consider off-diags*/
          A->DirichletRows[i] = 'F';
          break;
        }
    }

    ML_free(cols);
    ML_free(vals);
  } /*if (A->DirichletRows == NULL)*/

  return A->DirichletRows;

} /*ML_Operator_IdentifyDirichletRows()*/

/******************************************************/
/* Create an empty block matrix. This matirx will     */
/* be populated with ML_Operator_BlkMatInsert() and   */
/* then processed with ML_Operator_BlkMatFinalize()   */
/*                                                    */
/* Input:                                             */
/*      NBlockRows     number of blocks rows in matrix*/
/*                                                    */
/*      NBlockCols     number of blocks columns       */
/*                                                    */
/*      destroy_level  indicates whether subblock     */
/*                     destructors are invoked when   */
/*                     from BlkMat's destructor.      */
/*                                                    */
/*                     ML_DESTROY_SHALLOW             */
/*                     ML_CLEAN_EVERYTHING            */
/*                     ML_DESTROY_EVERYTHING          */
/*                                                    */
/* Output:                                            */
/*      BlkMat         Empty block matrix of dimension*/
/*                         NBlockRows x NBlockCols    */
/*                                                    */
/******************************************************/
int  ML_Operator_BlkMatInit(ML_Operator *BlkMat, ML_Comm *comm,
       int NBlockRows, int NBlockCols, int destroy_level)
{
  struct  MLBlkMat  *widget;
  int     i;


  /* fill data structure */ 

  ML_Operator_Clean(BlkMat);
  ML_Operator_Init(BlkMat,comm);

  widget = (struct MLBlkMat *) ML_allocate(sizeof(struct MLBlkMat));
  widget->final_called = 0;
  widget->NBlockRows   = NBlockRows;
  widget->NBlockCols   = NBlockCols;
  widget->RowStart     = (int *) ML_allocate(sizeof(int)*(NBlockRows+1));
  widget->ColStart     = (int *) ML_allocate(sizeof(int)*(NBlockCols+1));
  for (i=0; i <= NBlockRows; i++) widget->RowStart[i] = -1;
  for (i=0; i <= NBlockCols; i++) widget->ColStart[i] = -1;
  widget->matrix     = (ML_Operator **) ML_allocate((NBlockRows*NBlockCols+1)*
                                                  sizeof(ML_Operator *));
  widget->matdata    = (ML_BlkMatData **) ML_allocate((NBlockRows*NBlockCols+1)*
                                           sizeof(ML_BlkMatData *));

  for (i = 0; i < NBlockRows*NBlockCols; i++) { 
     widget->matrix[i] = NULL;
     widget->matdata[i] = NULL;
  }
  widget->destroy_level = destroy_level;


  BlkMat->data_destroy = ML_Operator_BlkMatDestroy;
 
  ML_Operator_Set_ApplyFuncData(BlkMat, 0, 0,
                                widget, 0,
                                ML_Operator_BlkMatMatvec,0);
  return(0);
}


/***************************************************************/
/* Insert Entry into the (Row,Col)th position of BlkMat. It is */
/* assummed that ML_Operator_BlkMatInit(BlkMat) has been       */
/* invoked but not ML_Operator_BlkMatFinalize(BlkMat).         */
/***************************************************************/
int ML_Operator_BlkMatInsert(ML_Operator *BlkMat, ML_Operator *Entry,
                        int Row, int Col)
{
  struct  MLBlkMat  *widget;
  int i;

  widget = (struct MLBlkMat *) BlkMat->data;

      
  if (widget == NULL) pr_error("ML_Operator_BlkMatInsert: empty data field\n");
  if (widget->final_called) pr_error("ML_Operator_BlkMatInsert: finalize already invoked\n");
  if (   Row < 0    ) pr_error("ML_Operator_BlkMatInsert: Row index < 0\n"); 
  if (   Col < 0    ) pr_error("ML_Operator_BlkMatInsert: Col index < 0\n");
  if (Row >= widget->NBlockRows) pr_error("ML_Operator_BlkMatInsert: Row index too large\n");
  if (Col >= widget->NBlockCols) pr_error("ML_Operator_BlkMatInsert: Col index too large\n");
  i = Row+Col*widget->NBlockRows;
  if (widget->matrix[i] != NULL) pr_error("ML_Operator_BlkMatInsert: matrix entry already assigned\n");
  widget->matrix[i] = Entry;
  widget->matdata[i] = (ML_BlkMatData *) ML_allocate(
                                           sizeof(ML_BlkMatData));
  return(0);
}

/***************************************************************/
/* Extract an entry from the (Row,Col)th position of BlkMat. It*/
/* is assummed that ML_Operator_BlkMatInit(BlkMat) and         */
/* ML_Operator_BlkMatFinalize(BlkMat) have been invoked.       */
/***************************************************************/
ML_Operator *ML_Operator_BlkMatExtract(ML_Operator *BlkMat, int Row, int Col)
{
  struct  MLBlkMat  *widget;
  int i;
  ML_Operator *Entry;

  widget = (struct MLBlkMat *) BlkMat->data;

      
  if (widget == NULL) pr_error("ML_Operator_BlkMatExtract: empty data field\n");
  if (widget->final_called==0) pr_error("ML_Operator_BlkMatExtract: finalize not invoked\n");
  if (   Row < 0    ) pr_error("ML_Operator_BlkMatExtract: Row index < 0\n"); 
  if (   Col < 0    ) pr_error("ML_Operator_BlkMatExtract: Col index < 0\n");
  if (Row >= widget->NBlockRows) pr_error("ML_Operator_BlkMatExtract: Row index too large\n");
  if (Col >= widget->NBlockCols) pr_error("ML_Operator_BlkMatExtract: Col index too large\n");
  i = Row+Col*widget->NBlockRows;
  return(widget->matrix[i]);
}


/*************************************************************/
/* Finish the block matrix definition. At this point it is   */
/* assumed that ML_Operator_BlkMatInit() has been invoked    */
/* and that the matrix has been populated via the function   */
/* ML_Operator_BlkMatInsert(). To complete the definition of */
/* the matrix we need to figure out the row/column           */
/* dimensions of the matrix. The main thing, however, that   */
/* must be computed is a ML_ComminfoOp structure and a       */
/* renumbering of indiviual rows within subblocks so that    */
/* they correspond to a proper numbering in the total system.*/
/* This is complicated by the presence of ghost nodes and the*/
/* fact that ML has no notion of global numbering. Thus, this*/
/* function first assigns a global numbering and then uses   */
/* this global numbering to populate the ML_ComminfoOp       */
/* structure and renumber ghost columns.                     */
int  ML_Operator_BlkMatFinalize(ML_Operator *BlkMat)
{
  struct MLBlkMat *widget;
  ML_BlkMatData   *SubData;
  ML_CommInfoOP        *getrow_comm;
  ML_Operator          *SubMat;
  USR_REQ              *request;

  int NBlkRows, NBlkCols, *RowStart, *ColStart, *LocalId, *MaxCols;
  int MasterGhostIndex=0, *MasterGhostGid;
  int NneighTotal, type, *NrcvFromProc, *NSndToProc, *MyNeighbors, **SendIds;
  int ColOffset = 0, *GhostCopy, *Indices, *Offsets;
  int NNonGhosts, NGhosts, *OwnProc;
  int i,j,k,kk, count, *ibuf;
  int    *GlobalId;
  double *dGlobalId; /* ML_exchange_bdry() only works with doubles */
                     /* so this is need temporarily.               */


  widget = (struct MLBlkMat *) BlkMat->data;
  if (widget == NULL) return(1);
  if (widget->matrix == NULL) return(1);


  NBlkRows = widget->NBlockRows;
  NBlkCols = widget->NBlockCols;
  RowStart = widget->RowStart;
  ColStart = widget->ColStart;

  /* Check that all sub-block dimensions are consistent. That is        */
  /*                                                                    */
  /*  1) forall i,k  SubBlk(i,j)->invec_leng =  SubBlk(k,j)->invec_leng */
  /*                                                                    */
  /*  2) forall j,k SubBlk(i,j)->outvec_leng = SubBlk(i,k)->outvec_leng */
  /*                                                                    */
  /* Additionally, determine the starting local indices for each        */
  /* SubBlk(i,j) and store them in (RowStart[i],ColStart[j]).           */
  /*                                                                    */
  /* This is done by first putting the block dimensions in (RowStart[i],*/
  /* ColStart[j]) and then converting them to offsets.                  */


  for (i =0; i < NBlkRows; i++) {
     for (j =0; j < NBlkCols; j++) {
        SubMat = widget->matrix[i+j*NBlkRows];
        if (SubMat != NULL) {
           if (RowStart[i] == -1) RowStart[i] = SubMat->outvec_leng;
           else {
              if (RowStart[i] != SubMat->outvec_leng) {
                 printf("ML_Operator_BlkMatFinalize(%d): The (%d,%d)th\n",BlkMat->comm->ML_mypid,i,j);
                 printf("block has %d rows as opposed to %d rows from\n",
                        SubMat->outvec_leng, RowStart[i]);
                 printf("an earlier block.\n");
                 exit(1);
              }
           }
           if (ColStart[j] == -1) ColStart[j] = SubMat->invec_leng;
           else {
              if (ColStart[j] != SubMat->invec_leng) {
                 printf("ML_Operator_BlkMatFinalize(%d): The (%d,%d)th\n",BlkMat->comm->ML_mypid,i,j);
                 printf("block has %d cols as opposed to %d cols from\n",
                        SubMat->invec_leng, ColStart[j]);
                 printf("an earlier block.\n");
                 exit(1);
              }
           }
        }
     }
  }
  for (i =0; i < NBlkRows; i++) 
     if (RowStart[i] == -1) RowStart[i] = 0;

  for (j =0; j < NBlkCols; j++) 
     if (ColStart[j] == -1) ColStart[j] = 0;

  k = RowStart[0];
  RowStart[0] = 0;
  for (i =0; i < NBlkRows; i++) {
     kk = RowStart[i+1];
     RowStart[i+1] = RowStart[i] + k;
     k = kk;
  }
    
  k = ColStart[0];
  ColStart[0] = 0;
  for (j =0; j < NBlkCols; j++) {
     kk = ColStart[j+1];
     ColStart[j+1] = ColStart[j] + k;
     k = kk;
  }

  widget->outvec = RowStart[NBlkRows];
  widget->invec = ColStart[NBlkCols];
  widget->final_called = 1;

  /*********************************************************/
  /* Count the total number of possible nonghosts, ghosts, */
  /* and neighboring processors in all submatrices.        */
  /*********************************************************/

  NneighTotal = 0;
  NGhosts     = 0;
  NNonGhosts  = 0;
  for (j = 0; j < NBlkCols; j++) {
    for (i = 0; i < NBlkRows; i++) {
      SubMat  = widget->matrix[i+j*NBlkRows];
      if (SubMat != NULL) {
          NNonGhosts += SubMat->invec_leng;
          if (SubMat->getrow->pre_comm != NULL) {
             NneighTotal += ML_CommInfoOP_Get_Nneighbors(SubMat->getrow->pre_comm);
             NGhosts += ML_CommInfoOP_Compute_TotalRcvLength(SubMat->getrow->pre_comm);
          }
      }
    }
  }
       
  /***************************************************************************/
  /* 1) Assign global ids to each subblock:                                  */
  /*                                                                         */
  /*   global id = k + ML_mypid*MaxCols[j] + ColOffset                       */
  /*                                                                         */
  /*   where                                                                 */
  /*     MaxCols[j] is the maximum # of point columns in the jth block column*/
  /*     ML_mypid is the processor id.                                       */
  /*     ColOffset is sum of MaxCols(0:j-1)*Nprocs                           */
  /*     k is the local index.                                               */
  /*                                                                         */
  /* 2) Communicate so that ghost dofs also have global ids available.       */ 
  /*                                                                         */
  /* 3) Build list (MasterGhostGid) holding all the ghost global ids for all */
  /*    the blocks as well as OwnProc which records the owning processor.    */
  /*                                                                         */
  /***************************************************************************/


  MasterGhostGid    = (int *) ML_allocate(sizeof(int)*(NGhosts+1));
  OwnProc           = (int *) ML_allocate(sizeof(int)*(NGhosts+1));
  Offsets = (int *) ML_allocate(sizeof(int)*(NBlkCols+1));
  MaxCols = (int *) ML_allocate(sizeof(int)*(NBlkCols+1));

  for (j = 0; j < NBlkCols; j++) {
    MaxCols[j] = ML_gmax_int(ColStart[j+1]-ColStart[j], BlkMat->comm);
    for (i = 0; i < NBlkRows; i++) {
      SubMat  = widget->matrix[i+j*NBlkRows];
      SubData = widget->matdata[i+j*NBlkRows];
      if (SubMat != NULL) {
         if (SubMat->getrow->pre_comm == NULL) 
           SubMat->getrow->pre_comm = ML_CommInfoOP_Create();


         dGlobalId= (double *) ML_allocate(sizeof(double)*(1+SubMat->invec_leng+
                                   SubMat->getrow->pre_comm->total_rcv_length));
         GlobalId= (int *) ML_allocate(sizeof(int)*(1+SubMat->invec_leng +
                                   SubMat->getrow->pre_comm->total_rcv_length));
         SubData->GlobalId = GlobalId;
         for (k = 0; k < SubMat->invec_leng; k++) {
           GlobalId[k] = k + BlkMat->comm->ML_mypid*MaxCols[j] + ColOffset;
           dGlobalId[k]= (double) GlobalId[k];
         }
         if (SubMat->getrow->pre_comm != NULL) 
            ML_exchange_bdry(dGlobalId, SubMat->getrow->pre_comm, SubMat->invec_leng, 
                             SubMat->comm, ML_OVERWRITE,NULL);
         for (k = 0; k < SubMat->getrow->pre_comm->total_rcv_length; k++) {
           GlobalId[k+SubMat->invec_leng] = (int) dGlobalId[k+SubMat->invec_leng];
           OwnProc[MasterGhostIndex] = (GlobalId[k+SubMat->invec_leng] - ColOffset)/MaxCols[j];
           MasterGhostGid[MasterGhostIndex++] = GlobalId[k+SubMat->invec_leng];
         }
         ML_free(dGlobalId);
      }
    }
    Offsets[j] = ColOffset;
    ColOffset += (MaxCols[j]*BlkMat->comm->ML_nprocs);
  }
  Offsets[NBlkCols] = ColOffset;

  /*********************************************************/
  /* Sort and remove duplicates from MasterGhostGid.       */
  /* Note: MasterGhostGid is a little messy because we     */
  /*       actually have 2 arrays (MasterGhostGid,OwnProc) */
  /*       We first sort on OwnProc to force contiguous    */
  /*       ordering of ghost nodes from the same processor */
  /*       (which is an ML requirement). We then sort      */
  /*       within a processor group so that these GIDs     */
  /*       are in order. Finally, we remove duplicates     */
  /*       but we can't use ML_rm_duplicates() because we  */
  /*       have to update both MasterGhostGid & OwnProc.   */
  /*********************************************************/

  ML_az_sort(OwnProc,MasterGhostIndex,MasterGhostGid,NULL);

  k = 0;  kk = 1;
  for (kk = 1; kk < MasterGhostIndex; kk++) {
    if (OwnProc[kk-1] != OwnProc[kk]) {
      ML_az_sort(&(MasterGhostGid[k]),kk-k,&(OwnProc[k]),NULL);
      k = kk;
    }
    else if (kk == MasterGhostIndex-1) {
      ML_az_sort(&(MasterGhostGid[k]),kk-k+1,&(OwnProc[k]),NULL);
      k = kk;
    }
  }
  
  /* remove duplicates */
  kk = 0;
  for (k = 1; k < MasterGhostIndex; k++) {
    if (MasterGhostGid[kk] != MasterGhostGid[k]) {
      kk++;
      MasterGhostGid[kk] = MasterGhostGid[k];
      OwnProc[kk] = OwnProc[k];
    }
  }
  if (MasterGhostIndex != 0) kk++;
  MasterGhostIndex = kk;


  /***************************************************/
  /* Make a copy of ghost gids and then sort         */
  /* This sorted list will be used when identifying  */
  /* local ids for ghost nodes.                      */
  /***************************************************/

  GhostCopy = (int *) ML_allocate(sizeof(int)*(MasterGhostIndex+1));
  Indices   = (int *) ML_allocate(sizeof(int)*(MasterGhostIndex+1));
  for (kk = 0; kk < MasterGhostIndex; kk++) {
    GhostCopy[kk] = MasterGhostGid[kk];
    Indices[kk] = kk;
  }
  ML_az_sort(GhostCopy,MasterGhostIndex,Indices,NULL);


  /************************************************************************/
  /* Compute the local id for each ghost global id by searching GhostCopy */
  /************************************************************************/

  for (j = 0; j < NBlkCols; j++) {
    for (i = 0; i < NBlkRows; i++) {
      SubMat   = widget->matrix[i+j*NBlkRows];
      SubData  = widget->matdata[i+j*NBlkRows];
      if (SubMat != NULL) {
         LocalId = (int *)  ML_allocate(sizeof(int)*(1+
                                   SubMat->getrow->pre_comm->total_rcv_length));
         SubData->LocalId = LocalId;
         GlobalId = SubData->GlobalId;
         for (k = 0; k < SubMat->getrow->pre_comm->total_rcv_length;k++) {
            kk = ML_find_index( GlobalId[k+SubMat->invec_leng], GhostCopy, MasterGhostIndex);
            if (kk == -1) {
               printf("ML_Operator_BlkMatFinalize: Global Id got lost?\n");
               exit(1);
            }
            else LocalId[k] = Indices[kk]+widget->invec;
         }
      }
    }
  }
  ML_free(Indices);
  ML_free(GhostCopy);

  /*********************************************************/
  /* Store neighboring procs, sort, and rm duplicates.     */
  /* Finally, set the neighbors in an ML_CommInfoOP object */
  /*********************************************************/

  MyNeighbors = (int *) ML_allocate(sizeof(int)*(NneighTotal+1));

  count = 0;
  for (j=0; j < NBlkCols; j++) {
     for (i=0; i < NBlkRows; i++) {
        SubMat = widget->matrix[i+j*NBlkRows];
        if ( (SubMat != NULL) && (SubMat->getrow->pre_comm != NULL)) {
           getrow_comm = SubMat->getrow->pre_comm;
           ibuf = ML_CommInfoOP_Get_neighbors(getrow_comm);
           for (kk = 0; kk < ML_CommInfoOP_Get_Nneighbors(getrow_comm); kk++) {
              MyNeighbors[count++] = ibuf[kk];
           }
           ML_free(ibuf);
        }
     }
  }
  ML_az_sort(MyNeighbors,NneighTotal,NULL,NULL);
  ML_rm_duplicates(MyNeighbors, &NneighTotal);

  ML_CommInfoOP_Set_neighbors(&(BlkMat->getrow->pre_comm), NneighTotal, 
                   MyNeighbors, ML_OVERWRITE, NULL, 0);


  /*******************************************************************/
  /* Count the number of items that we receive from each processor   */
  /* Send this number to the corresponding neighbors so that now     */
  /* each processor knows how much data it needs to send and receive */
  /* with each neigbhors.                                            */
  /*******************************************************************/

  NrcvFromProc = (int      *) ML_allocate(sizeof(int    )*(NneighTotal+1));
  NSndToProc   = (int      *) ML_allocate(sizeof(int    )*(NneighTotal+1));
  request      = (USR_REQ  *) ML_allocate(sizeof(USR_REQ)*(NneighTotal+1));

  for (kk = 0; kk < NneighTotal; kk++) NrcvFromProc[kk] = 0;
  kk = 0; 
  for (k = 0; k < NneighTotal; k++) {
    while ( (kk < MasterGhostIndex) && (OwnProc[kk] == MyNeighbors[k])) {
       NrcvFromProc[k]++;
       kk++;
    }
  }
  ML_free(OwnProc);

  type = 6234;

  for (i = 0; i < NneighTotal; i++)   /*** post receives ***/
    BlkMat->comm->USR_irecvbytes((void *) &(NSndToProc[i]), sizeof(int),
                &(MyNeighbors[i]), &type, BlkMat->comm->USR_comm, request+i);
  for (i = 0; i < NneighTotal; i++)   /*** send data *******/ 
    BlkMat->comm->USR_sendbytes((void *) &(NrcvFromProc[i]), sizeof(int), 
                MyNeighbors[i], type, BlkMat->comm->USR_comm);
  for (i = 0; i < NneighTotal; i++)  /*** recv data ********/ 
    BlkMat->comm->USR_cheapwaitbytes((void *) &(NSndToProc[i]), sizeof(int),
                &(MyNeighbors[i]), &type, BlkMat->comm->USR_comm, request+i);

  /*******************************************************************/
  /* Send the GIDs corresponding to ghost nodes so that now each     */
  /* processor knows the GIDs of nodes it must send to each neighbor.*/
  /*******************************************************************/

  SendIds = (int **) ML_allocate(sizeof(int *)*(NneighTotal+1));
  for (kk = 0; kk < NneighTotal; kk++)
     SendIds[kk] = (int *) ML_allocate(sizeof(int)*(NSndToProc[kk]+1));

  type = 5123;
  for (i = 0; i < NneighTotal; i++)    /*** post receives ***/
    BlkMat->comm->USR_irecvbytes((void *) SendIds[i], (unsigned int)
              (sizeof(int)*NSndToProc[i]),&(MyNeighbors[i]),
              &type,BlkMat->comm->USR_comm, request+i);
  kk = 0;
  for (i = 0; i < NneighTotal; i++) { /*** send data *******/ 
    BlkMat->comm->USR_sendbytes((void *) &(MasterGhostGid[kk]),
                sizeof(int)*NrcvFromProc[i], MyNeighbors[i],
                type, BlkMat->comm->USR_comm);
    kk += NrcvFromProc[i];
  }
  for (i = 0; i < NneighTotal; i++)  /*** recv data *******/ 
    BlkMat->comm->USR_cheapwaitbytes((void *) SendIds[i], (unsigned int)
              (sizeof(int)*NSndToProc[i]),&(MyNeighbors[i]),
              &type, BlkMat->comm->USR_comm, request+i);

  ML_free(request);
  ML_free(MasterGhostGid);

  /*************************************************************/
  /* Convert global ids back to local ids and set the exchange */
  /* information in BlkMat->getrow->pre_comm. Recall from the  */
  /* above comment that                                        */
  /*                                                           */
  /*   global id = k + ML_mypid*MaxCols[j] + ColOffset         */
  /*************************************************************/

  count = widget->invec;
  for (kk = 0; kk < NneighTotal; kk++) {
    for (k = 0; k < NSndToProc[kk]; k++) {
      i = 0; 
      while ( SendIds[kk][k]  >= Offsets[i])  i++; /* linear search */
      SendIds[kk][k] = SendIds[kk][k] - BlkMat->comm->ML_mypid*MaxCols[i-1] +
                                      ColStart[i-1] - Offsets[i-1];
    }
    ibuf = (int *) ML_allocate(sizeof(int)*(NrcvFromProc[kk]+1));
    for (k = 0; k < NrcvFromProc[kk]; k++) ibuf[k] = count++;
    ML_CommInfoOP_Set_exch_info(BlkMat->getrow->pre_comm, MyNeighbors[kk],
                       NrcvFromProc[kk], ibuf, NSndToProc[kk], SendIds[kk]);
    ML_free(ibuf);
    ML_free(SendIds[kk]);
  }
  ML_free(SendIds);
  ML_free(NSndToProc);
  ML_free(NrcvFromProc);
  ML_free(MyNeighbors);
  ML_free(MaxCols);
  ML_free(Offsets);


  BlkMat->data = NULL;  /* ML_Operator_Set_ApplyFuncData() destroys  */
                        /* BlkMat->data before setting it to widget  */
                        /* In this case, BlkMat->data already points */
                        /* to widget which should not be destroyed   */

  ML_Operator_Set_ApplyFuncData(BlkMat, widget->invec, widget->outvec,
                                widget, widget->outvec,
                                ML_Operator_BlkMatMatvec,0);

  ML_Operator_Set_Getrow(BlkMat, widget->outvec, ML_Operator_BlkMatGetrow);


}


/*******************************************************/
/* Data specific destructor for a block dense matrix   */
/* created via ML_Operator_BlkMatInit() and            */
/* ML_Operator_BlkMatFinalize()                        */
/*                                                     */
/* NOTE: widget->destroy_level indicates what part of  */
/* widget->matrix is freed.                            */
/*                                                     */
/*    ML_DESTROY_SHALLOW :   free(matrix)              */
/*                                                     */
/*    ML_DESTROY_EVERYTHING: destroy(matrix[i]);free(matrix)*/
/*                                                     */
/*    ML_CLEAN_EVERYTHING: clean(matrix);free(matrix)  */
/*                                                     */
/*******************************************************/
void  ML_Operator_BlkMatDestroy(void *data)
{
  struct  MLBlkMat  *widget;
  int i;

  widget = (struct MLBlkMat *) data;
  if (widget != NULL) {
    for (i = 0; i < widget->NBlockRows*widget->NBlockCols; i++) {
      if ( widget->matdata[i] != NULL) {
        ML_free(widget->matdata[i]->LocalId); 
        ML_free(widget->matdata[i]->GlobalId); 
        ML_free(widget->matdata[i]); 
      }
    }
    if (widget->destroy_level == ML_DESTROY_EVERYTHING) {
      for (i = 0; i < widget->NBlockRows*widget->NBlockCols; i++) {
        ML_Operator_Destroy( &(widget->matrix[i]) );
      }
    }
    else if (widget->destroy_level == ML_CLEAN_EVERYTHING) {
      for (i = 0; i < widget->NBlockRows*widget->NBlockCols; i++) {
        ML_Operator_Clean( widget->matrix[i] );
      }
    }
    else if (widget->destroy_level != ML_DESTROY_SHALLOW ) {
      printf("ML_Operator_BlkMatDestroy: Unknown destroy option?\n");
    }
    ML_free(widget->matdata);
    ML_free(widget->matrix);
    ML_free(widget->RowStart);
    ML_free(widget->ColStart);
    ML_free(widget);
  }
}

/**************************************************************/
/* Matvec corresponding to a block matrix created via         */
/* ML_Operator_BlkMatInit() and ML_Operator_BlkMatFinalize(). */
/**************************************************************/
int ML_Operator_BlkMatMatvec(ML_Operator *BlkMat, int ilen,
        double p[], int olen, double ap[])
{
  struct MLBlkMat *widget;
  ML_Operator **matrix;
  double *temp, *ptr;
  ML_Operator *SubMat;
  int NBlkRows, NBlkCols, i, j, k, *RowStart, *ColStart;


  widget = (struct MLBlkMat *) BlkMat->data;
  temp = (double *) ML_allocate(sizeof(double)*(BlkMat->outvec_leng+1));
  for (i = 0; i < BlkMat->outvec_leng; i++) ap[i] = 0.;

  NBlkRows = widget->NBlockRows;
  NBlkCols = widget->NBlockCols;
  RowStart = widget->RowStart;
  ColStart = widget->ColStart;
  matrix   = widget->matrix;

  for (i =0; i < NBlkRows; i++) {
     for (j =0; j < NBlkCols; j++) {
        SubMat = matrix[i+j*NBlkRows];
        if (SubMat != NULL) {
           ML_Operator_Apply(SubMat, SubMat->invec_leng,
            &(p[ColStart[j]]), SubMat->outvec_leng, temp);
           ptr = &(ap[RowStart[i]]);
           for (k = 0; k < SubMat->outvec_leng; k++) 
              ptr[k] += temp[k];
        }

     }
  }
  ML_free(temp);
  return(0);
}

/**************************************************************/
/* Getrow corresponding to a block matrix created via         */
/* ML_Operator_BlkMatInit() and ML_Operator_BlkMatFinalize(). */
/***************************************************************/
int ML_Operator_BlkMatGetrow(ML_Operator *BlkMat, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[],
   double values[], int row_lengths[])
{
  struct MLBlkMat *widget;
  ML_Operator *SubMat, **matrix;
  ML_BlkMatData   **matdata;
  int NBlkRows, NBlkCols, *RowStart, *ColStart;
  int BlockRow, row, *LocalId, i, j, k, SubInvecLeng;

  widget = (struct MLBlkMat *) BlkMat->data;
  NBlkRows = widget->NBlockRows;
  NBlkCols = widget->NBlockCols;
  RowStart = widget->RowStart;
  ColStart = widget->ColStart;
  matrix   = widget->matrix;
  matdata  = widget->matdata;

  /***************************************************/
  /* This is kind of a slow linear search. Shouldn't */
  /* matter if we have small block systems.          */
  /***************************************************/

  BlockRow = 0;
  while ( *requested_rows >= RowStart[BlockRow+1]) BlockRow++;

  row = *requested_rows - RowStart[BlockRow];
  k   = 0;
  for (j =0; j < NBlkCols; j++) {
     SubMat = matrix[BlockRow+j*NBlkRows];
     if (SubMat != NULL) {
        LocalId = matdata[BlockRow+j*NBlkRows]->LocalId;

        SubInvecLeng = SubMat->invec_leng;
        i = ML_Operator_Getrow(SubMat, 1, &row, allocated_space-k, 
                  &(columns[k]), &( values[k]), row_lengths);
        if (i != 1) return(0);
        /******************************************************/
        /* convert the local ids (both nonghost and ghost)    */
        /* within a block to local ids for the large system.  */
        /******************************************************/
        for (i = 0; i < *row_lengths ; i++) {
           if (columns[k+i] <  SubInvecLeng) columns[k+i] += ColStart[j];
           else columns[k+i] = LocalId[columns[k+i]-SubInvecLeng];
        }
        k += (*row_lengths);
     }
  }
  *row_lengths = k;
  return(1);
}

/*******************************************************************/
/* Like ML_rap but for block matrices. Form the explicit product   */
/* Rmat*Amat*Pmat and store the answer in Result.  Basically, we   */
/* invoke ML_rap() on the individual block matrices to get the     */
/* multiplication done. Currently, this function ASSUMES that both */
/* Rmat and Pmat are  block diagonal.                              */
/*                                                                 */
/* Note: there are cases where certain diagonal subblock have      */
/*       probably already been computed. It would be nice to add   */
/*       a capability so that these ones could be skipped. This    */
/*       would require some ability to keep track of what needs to */
/*       be destroyed when the destructor is eventually invoked.   */
/*******************************************************************/


int ML_Blkrap(ML_Operator *Rmat, ML_Operator *Amat, ML_Operator *Pmat, 
              ML_Operator *Result, int matrix_type)
{
  struct MLBlkMat    *Rwidget, *Pwidget, *Awidget;
  int    RNBlkRows,ANBlkRows,PNBlkRows,RNBlkCols,ANBlkCols,PNBlkCols;
  ML_Operator *RSubMat, *ASubMat, *PSubMat, *SubMat; 
  ML_Operator **Rmatrix, **Amatrix, **Pmatrix;
  int         i, j;


  /* Get all the block matrix information */

  Rwidget   = (struct MLBlkMat *) Rmat->data;
  RNBlkRows = Rwidget->NBlockRows;
  RNBlkCols = Rwidget->NBlockCols;
  Rmatrix   = Rwidget->matrix;

  Pwidget   = (struct MLBlkMat *) Pmat->data;
  PNBlkRows = Pwidget->NBlockRows;
  PNBlkCols = Pwidget->NBlockCols;
  Pmatrix   = Pwidget->matrix;

  Awidget   = (struct MLBlkMat *) Amat->data;
  ANBlkRows = Awidget->NBlockRows;
  ANBlkCols = Awidget->NBlockCols;
  Amatrix   = Awidget->matrix;

  if (RNBlkCols != ANBlkRows) {
    printf("Error ML_BlkRAP: Incompatible matrices. R must have the same\n");
    printf("                 number of block columns as A has block rows.\n");
    return 1;
  }
  if (PNBlkRows != ANBlkCols) {
    printf("Error ML_BlkRAP: Incompatible matrices. P must have the same\n");
    printf("                 number of block rows as A has block columns.\n");
    return 1;
  }

  /* Only matmat mult for block diagonal R and P are currently */
  /* implemented so we better check that this is what we have. */

  for (i =0; i < RNBlkRows; i++) {
     for (j =0; j < RNBlkCols; j++) {
        if ( (i != j) && (Rmatrix[i+j*RNBlkRows] != NULL)) {
          printf("Error ML_BlkRAP: R must be block diagonal\n");
          return 1;
        }
     }
  }

  for (i =0; i < PNBlkRows; i++) {
     for (j =0; j < PNBlkCols; j++) {
        if ( (i != j) && (Pmatrix[i+j*PNBlkRows] != NULL)) {
          printf("Error ML_BlkRAP: P must be block diagonal\n");
          return 1;
        }
     }
  }

  ML_Operator_BlkMatInit(Result, Rmat->comm, RNBlkRows, PNBlkCols, ML_DESTROY_EVERYTHING);

  for (i =0; i < RNBlkRows; i++) {
     RSubMat = Rmatrix[i+i*RNBlkRows];
     for (j =0; j < PNBlkCols; j++) {
       PSubMat = Pmatrix[j+j*PNBlkRows];
       ASubMat = Amatrix[i+j*ANBlkRows];
       if ( (RSubMat != NULL) && (ASubMat != NULL) && (PSubMat != NULL)) {
         SubMat = ML_Operator_Create(Rmat->comm);
         ML_rap(RSubMat, ASubMat, PSubMat, SubMat, matrix_type);
         ML_Operator_BlkMatInsert(Result, SubMat, i, j);
       }
     }
  }
  ML_Operator_BlkMatFinalize(Result);

  return 0;
}
/*****************************************************************/
/* Take Matrix and project it from level 'FromLevel' to the next */
/* coarser level using the multigrid transfers associated with   */
/* BlockLocation in mlptr.                                       */
/*                                                               */
/* Note: Right now ProjectMe() only does things corresponding to */
/* to the block diagonal. It would be pretty easy to change this */
/* for off-diagonals which would require a BlockRowLocation and  */
/* a BlockColLocation. I'm not sure if this is really needed?    */
/*****************************************************************/
ML_Operator *ProjectMe(ML *mlptr, int BlockLocation, int FromLevel, 
                       ML_Operator *Matrix,int matrix_type)
{
  ML_Operator *Answer, *CompR, *CompP, *Rmat, *Pmat;
  /* Get composite R and P */

  CompR = &(mlptr->Rmat[FromLevel]);
  CompP = CompR->to->Pmat;

  /* Now get the blks we want */

  Rmat = ML_Operator_BlkMatExtract(CompR, BlockLocation, BlockLocation);
  Pmat = ML_Operator_BlkMatExtract(CompP, BlockLocation, BlockLocation);
  Answer = ML_Operator_Create(Rmat->comm);
  ML_rap(Rmat, Matrix, Pmat, Answer, matrix_type);

  return(Answer);
}

