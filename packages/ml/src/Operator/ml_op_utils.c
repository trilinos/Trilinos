/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_struct.h"
#include "ml_op_utils.h"
#include "ml_agg_genP.h"
#include "ml_memory.h"

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
 * -------------------------------------------------------------------- */

int ML_Gen_Restrictor_TransP(ML *ml_handle, int level, int level2)
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


   Pmat  = &(ml_handle->Pmat[level2]);
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
      flag = getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
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
      getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
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
#include "metis.h"
#endif
int ML_Operator_BlockPartition(ML_Operator *matrix, int nLocalNd, int *nblk,
                         int *pnode_part,
                         int *ndwts /*=NULL*/, int *egwts/*=NULL*/,
                         int nedges /*= 0*/, int metis_or_parmetis )
{
#ifdef HAVE_ML_METIS
  int locid, ii, numadjac, *bindx = NULL;
  idxtype *xadj = NULL, *adjncy = NULL, *blks = NULL;
  int options[5]={0,3,1,1,0};
  int weightflag = ndwts ? 2 : 0;
  int nmbng = 0, edgecut = -1, n = nLocalNd, np = *nblk;
  double *val = NULL;
  int allocated = 0, row_length, j;
  /* FILE *fp; */

#if defined(HAVE_ML_PARMETIS_2x) || defined(HAVE_ML_PARMETIS_3x)
   int *map = NULL, itemp1, itemp2, nprocs, myid;
   idxtype ncon, *tpwts = NULL, *vtxdist = NULL;
   float ubvec, itr;
#endif

  if( egwts ) weightflag++;

  if( *nblk == 1 || ((nLocalNd < 1) && (metis_or_parmetis == ML_USEMETIS))){
    for( ii = 0 ; ii < nLocalNd ; ii++ ) pnode_part[ii] = 0;
    return 0;
  }

#if ! (defined(HAVE_ML_PARMETIS_2x) || defined(HAVE_ML_PARMETIS_3x))
  if (metis_or_parmetis == ML_USEPARMETIS) {
    for( ii = 0 ; ii < nLocalNd ; ii++ ) pnode_part[ii] = 0;
    return 1;
}
#endif

  /* Set 'xadj' & 'adjncy' adjacency data.  This is the graph
     that Metis requires.  It corresponds to the matrix graph
     without self connections (diagonal entries). */


  xadj = (idxtype *) ML_allocate( (nLocalNd+1) * sizeof(idxtype) );
  if (xadj == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");

  ii = 0;
  for( locid = 0 ; locid < nLocalNd ; locid++ )
  {
    ML_get_matrix_row(matrix, 1, &locid, &allocated, &bindx, &val,
                      &row_length, 0);
    ii += row_length - 1;
  }
  numadjac = ii;
  adjncy = (idxtype *) ML_allocate( (numadjac+1) * sizeof(idxtype) );
  if (adjncy == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");

  if (metis_or_parmetis == ML_USEMETIS) {
    ii = 0;
    for( locid = 0 ; locid < nLocalNd ; locid++ )   {
      xadj[locid] = ii;
      ML_get_matrix_row(matrix,1,&locid,&allocated,&bindx,&val,&row_length,0);
      for (j = 0; j < row_length; j++)     {
	if (( bindx[j] != locid ) & (bindx[j] < nLocalNd))
	  adjncy[ii++] = bindx[j];
      }
    }
    xadj[nLocalNd] = ii;
  }
#if  (defined(HAVE_ML_PARMETIS_2x) || defined(HAVE_ML_PARMETIS_3x))
  else {
    ML_create_unique_id(nLocalNd, &map, matrix->getrow->pre_comm, matrix->comm);
    ii = 0;
    for( locid = 0 ; locid < nLocalNd ; locid++ )   {
      xadj[locid] = ii;
      ML_get_matrix_row(matrix,1,&locid,&allocated,&bindx,&val,&row_length,0);
      for (j = 0; j < row_length; j++)     {
	if ( bindx[j] != locid ) {
	  adjncy[ii++] = map[bindx[j]];
	}
      }
    }
    xadj[nLocalNd] = ii;
    if (map != NULL) ML_free(map);
  }
#endif
  if (val   != NULL) ML_free(val);
  if (bindx != NULL) ML_free(bindx);

  /*
   n                number of vertices in graph
   xadj,adjncy      adjacency structure of graph
   vwgt,adjwgt      weights of vertices and edges
   weightflag          0 = no weights, 1 = weights on edges only
                    2 = weights on nodes only, 3 = weights on both
   numflag          0 = C style numbering, 1 = Fortran style numbering
   np               desired number of parts
   options          options for matching type, initial partitioning,
                    refinement, and debugging
                    If options[0] == 0, then default values are used.
   edgecut          returns number of edges cut by partition
   part             returns the partition
   
   The METIS authors suggest using METIS_PartGraphRecursive if the desired
   partition should have 8 or less parts.  They suggest using
   METIS_PartGraphKway otherwise.

   (See the METIS 4.0 manual for more details.)
  */

  /* Get local partition. */

  if (metis_or_parmetis == ML_USEMETIS) {
    if (np < 8)
      METIS_PartGraphRecursive( &n, xadj, adjncy, ndwts, egwts, &weightflag,
				&nmbng, &np, options, &edgecut, pnode_part );
    else
      METIS_PartGraphKway( &n, xadj, adjncy, ndwts, egwts, &weightflag, &nmbng,
			   &np, options, &edgecut, pnode_part );

    blks = (int *) ML_allocate((np+1)*sizeof(int));
    if (blks == NULL) pr_error("ML_Operator_BlockPartition: out of space\n");
    for (j = 0; j < *nblk; j++) blks[j] = -1;
    for (j = 0; j < n; j++) blks[pnode_part[j]] = -2;
    ii = 0;
    for (j = 0; j < *nblk; j++) {
      if ( blks[j] == -2) {
	blks[j] = ii;
	ii++;
      }
    }
    for (j = 0; j < n; j++) pnode_part[j] = blks[pnode_part[j]];
    *nblk = ii;
    if (blks    != NULL) ML_free(blks);
  }

#if  (defined(HAVE_ML_PARMETIS_2x) || defined(HAVE_ML_PARMETIS_3x))
  else {
    nprocs = matrix->comm->ML_nprocs;
    myid   = matrix->comm->ML_mypid;
    
    ubvec = 1.05;
    tpwts   = (idxtype *) ML_allocate(sizeof(idxtype)*ML_max(np,nprocs));
    vtxdist = (idxtype *) ML_allocate(sizeof(idxtype)*(nprocs+1));
    for (ii = 0; ii <= nprocs; ii++) vtxdist[ii] = 0;
    vtxdist[myid] = nLocalNd;
    ML_gsum_vec_int(&vtxdist,&tpwts,nprocs,matrix->comm);

    itemp1 = 0;
    for (ii = 0; ii < nprocs ; ii++) {
      itemp2 = vtxdist[ii];
      vtxdist[ii] = itemp1;
      itemp1 += itemp2;
    }
    vtxdist[nprocs] = itemp1;
    
    ncon = 1;
    itr = 1000.;
    options[0] = 0;

    for (ii = 0; ii < nLocalNd; ii++) pnode_part[ii] = matrix->comm->ML_mypid;

   weightflag = 0;
#if ( defined(HAVE_ML_PARMETIS_3x) )
   ParMETIS_V3_AdaptiveRepart( vtxdist,xadj,adjncy, NULL,
				NULL, NULL, &weightflag,
				&nmbng, &ncon, &np, NULL, &ubvec, 
				&itr, options, &edgecut, pnode_part,
				&(matrix->comm->USR_comm));
#else
   /* This routine does not seem to allow having a different number of */
   /* partitions and processors.         
   ParMETIS_RepartGDiffusion(vtxdist, xadj, adjncy, NULL, NULL,&weightflag,
			     &nmbng, options, &edgecut, pnode_part, 
			     &(matrix->comm->USR_comm));
   */
   ParMETIS_PartKway( vtxdist,xadj,adjncy, NULL, NULL, &weightflag,
			   &nmbng, &np, options, &edgecut, pnode_part,
			   &(matrix->comm->USR_comm));
#endif
    *nblk = np;
    if (vtxdist != NULL) ML_free(vtxdist);
    if (tpwts   != NULL) ML_free(tpwts);
  }
  if (adjncy  != NULL) ML_free(adjncy);
  if (xadj    != NULL) ML_free(xadj);

#endif

#else
  printf("ML_partitionBlocksNodes: Metis not linked\n");
  ML_avoid_unused_param( (void *) matrix);
  ML_avoid_unused_param( (void *) &nLocalNd);
  ML_avoid_unused_param( (void *) nblk);
  ML_avoid_unused_param( (void *) pnode_part);
  ML_avoid_unused_param( (void *) ndwts);
  ML_avoid_unused_param( (void *) egwts);
  ML_avoid_unused_param( (void *) &nedges);
#endif

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
      flag = getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
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
      getrow(data, 1, &i, Nghost+osize+1, colbuf, valbuf, &length);
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
 
  ML_Operator *eye1, *eye2;
 
  eye1 = ML_Operator_Create(A->comm);
  eye2 = ML_Operator_Create(A->comm);
 
  ML_Operator_Set_ApplyFuncData(eye1, A->invec_leng, A->invec_leng,
            NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye1, A->invec_leng, eye_getrows);
 
  ML_Operator_Set_ApplyFuncData(eye2, A->invec_leng, A->invec_leng,
            NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye2, A->invec_leng, eye_getrows);
  ML_2matmult(A, eye1, Atrans, ML_CSR_MATRIX);

  ML_Operator_Destroy(&eye1);
  ML_Operator_Destroy(&eye2);

 
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
  int allocated_space, *cols, i, j, n;
  double *vals, *tdiag;
   if (Amat->diagonal == NULL) 
   {
      if (Amat->getrow->func_ptr == NULL) 
         pr_error("Error(ML_Jacobi): Need diagonal\n");
      else 
      {
         allocated_space = 30;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
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
         ML_free(cols); ML_free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         ML_free(tdiag);
      } 
   }
   ML_DVector_GetDataPtr( Amat->diagonal, diagonal);
   return 0;
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

#include "limits.h"
#include "float.h"
#include "ml_cg.h"
/* ******************************************************************** */
/* "Cheap" analysis of the main properties of an ML_Operator.           */
/************************************************************************/
int ML_Operator_Analyze(ML_Operator * Op, char* name)
{

  int i,j;
  double MyFrobeniusNorm=0.0; 
  double FrobeniusNorm=0.0;
  double MyMinElement=DBL_MAX; 
  double MinElement=DBL_MAX;
  double MyMaxElement=DBL_MIN; 
  double MaxElement=DBL_MIN;
  double MyMinAbsElement=DBL_MAX; 
  double MinAbsElement=DBL_MAX;
  double MyMaxAbsElement=0.0; 
  double MaxAbsElement=0.0;
  int NumNonzeros; /* nonzero elements in a row */
  int MyMinNumNonzeros = INT_MAX;
  int MyMaxNumNonzeros = 0;
  int MinNumNonzeros = INT_MAX;
  int MaxNumNonzeros = 0;
  int MyTotalNumNonzeros = 0;
  int TotalNumNonzeros;
  int MyTotalNumGivenNonzeros = 0;
  int TotalNumGivenNonzeros;

  int NumMyRows = 0, NumGlobalRows = 0;
  
  int  allocated = 10;
  double* colVal = NULL;
  int*    colInd = NULL;
  int MyPID;
  int NumProc;

  double* Diagonal = NULL;
  double* SumOffDiagonal = NULL;
  int* IsDiagonallyDominant = NULL;
  double Element, AbsElement; /* generic nonzero element and its abs value */
  int ierr; /* return error code for ML_Operator_Getrow() */
  int NumDiagonallyDominant = 0;
  int NumMyDiagonallyDominant = 0;
  int NumWeaklyDiagonallyDominant = 0;
  int NumMyWeaklyDiagonallyDominant = 0;
  int NumMyDirichletRows = 0;
  int NumDirichletRows = 0;
  int NumMyLowerNonzeros = 0;
  int NumMyUpperNonzeros = 0;
  int NumLowerNonzeros;
  int NumUpperNonzeros;
  double Min;
  double MyMin;
  double Max;
  double MyMax;
  int Equation;
  double time;
  int NumGivenNonzeros;
  int offset;
  double* global_id = NULL;
  int MyMaxBandwidth = 0;
  int MaxBandwidth;
  int Nghost;
  int grow;
  int gcol;

  ML_Krylov data; /* for power-method */

  /* ---------------------- execution begins ------------------------------ */

  time = GetClock();

  /* returns if Op is not set, or Op is a rectangular matrix */
  if( Op == NULL ) return(-1);
  if( Op->invec_leng != Op->outvec_leng ) return(-2);

  /* set up data and allocate memory */
  MyPID = Op->comm->ML_mypid;
  NumProc = Op->comm->ML_nprocs;
  NumMyRows = Op->invec_leng;
  NumGlobalRows = ML_Comm_GsumInt(Op->comm,NumMyRows);

  colInd = (int*)   ML_allocate(sizeof(int)*allocated);
  colVal = (double*)ML_allocate(sizeof(double)*allocated);

  Diagonal = (double*)ML_allocate(sizeof(double)*NumMyRows);
  SumOffDiagonal = (double*)ML_allocate(sizeof(double)*NumMyRows);
  
  for (i=0 ; i<NumMyRows ; ++i) {
    SumOffDiagonal[i] = 0.0;
  }
 
  /* global column number */

  if (Op->getrow->pre_comm == NULL) Nghost = 0;
  else {
    if (Op->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Op->getrow->pre_comm);
    Nghost = Op->getrow->pre_comm->total_rcv_length;
  }

  global_id = (double*) ML_allocate(sizeof(double) * (NumMyRows + Nghost));
  if (global_id == NULL) {
    fprintf(stderr,
	    "Not enough memory\n");
    exit(EXIT_FAILURE);
  }
  offset = ML_gpartialsum_int(NumMyRows, Op->comm);

  for (i = 0 ; i < NumMyRows ; ++i) {
    global_id[i] = (double) (offset + i);
  }

  ML_exchange_bdry(global_id,Op->getrow->pre_comm, 
		   Op->invec_leng,Op->comm,ML_OVERWRITE,NULL);

  /* cycle over all matrix rows */

  for( i=0 ; i<NumMyRows ; ++i ) {

    grow = (int)global_id[i];

    ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NumGivenNonzeros);

    if( ierr == 0 ) {
      do {
	ML_free(colInd);
	ML_free(colVal);
	allocated *= 2;
	colInd = (int*)   ML_allocate(sizeof(int)*allocated);
	colVal = (double*)ML_allocate(sizeof(double)*allocated);

	ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NumGivenNonzeros);
      } while( ierr == 0 );
    }

    MyTotalNumGivenNonzeros += NumGivenNonzeros;

    /* compute the real number of nonzero elements */
    
    NumNonzeros = 0;
    for (j = 0 ; j < NumGivenNonzeros ; ++j) {
      /* compute the real number of nonzeros only */
      if (colVal[j] != 0.0) {
	++NumNonzeros;
	gcol = (int)global_id[colInd[j]];
        if (gcol < grow) 
	  NumMyLowerNonzeros++;
	else if (gcol > grow) 
	  NumMyUpperNonzeros++;
	/* compute bandwidth */
	if (abs(gcol - grow) > MyMaxBandwidth)
	  MyMaxBandwidth = abs(gcol - grow);
      }
    }

    MyTotalNumNonzeros += NumNonzeros;

    if (NumNonzeros > MyMaxNumNonzeros)
      MyMaxNumNonzeros = NumNonzeros;
    if (NumNonzeros < MyMinNumNonzeros)
      MyMinNumNonzeros = NumNonzeros;
    
    if (NumNonzeros == 1) NumMyDirichletRows++;
    
    /* start looking for min/max element (with and without abs()) */
    /* I consider the following:
     * - min element;
     * - min abs element;
     * - max element;
     * - max abs element;
     * - diagonal element;
     * - sum off abs off-diagonal elements;
     * Here I compute the local ones, with prefix `My'. The global
     * ones (without prefix will be computed later
     */
    for( j=0 ; j<NumNonzeros ; ++j ) {
      Element = colVal[j];
      AbsElement = ML_dabs(Element);
      if( Element<MyMinElement ) MyMinElement = Element;
      if( Element>MyMaxElement ) MyMaxElement = Element;
      if( AbsElement<MyMinAbsElement ) MyMinAbsElement = AbsElement;
      if( AbsElement>MyMaxAbsElement ) MyMaxAbsElement = AbsElement;
      if( colInd[j] == i ) Diagonal[i] = AbsElement;
      else
	SumOffDiagonal[i] += ML_dabs(Element);
      MyFrobeniusNorm += Element*Element;
    }
  } /* for over all matrix rows */

  /* compute the min/max of important quantities over all processes */

  MinElement    = - ML_Comm_GmaxDouble(Op->comm, -MyMinElement);
  MaxElement    =   ML_Comm_GmaxDouble(Op->comm, MyMaxElement);
  MinAbsElement = - ML_Comm_GmaxDouble(Op->comm, -MyMinAbsElement);
  MaxAbsElement =   ML_Comm_GmaxDouble(Op->comm, MyMaxAbsElement);
  
  FrobeniusNorm =   ML_Comm_GmaxDouble(Op->comm, MyFrobeniusNorm);
  NumDirichletRows = ML_Comm_GsumInt(Op->comm, NumMyDirichletRows);
  
  MinNumNonzeros   = - ML_Comm_GmaxInt(Op->comm, -MyMinNumNonzeros);
  MaxNumNonzeros   =   ML_Comm_GmaxInt(Op->comm, MyMaxNumNonzeros);
  TotalNumNonzeros =   ML_Comm_GsumInt(Op->comm, MyTotalNumNonzeros);
  TotalNumGivenNonzeros =   ML_Comm_GsumInt(Op->comm, MyTotalNumGivenNonzeros);
  NumLowerNonzeros =   ML_Comm_GsumInt(Op->comm, NumMyLowerNonzeros);
  NumUpperNonzeros =   ML_Comm_GsumInt(Op->comm, NumMyUpperNonzeros);
  MaxBandwidth     =   ML_Comm_GmaxInt(Op->comm, MyMaxBandwidth);

  /* a test to see if matrix is diagonally-dominant */

  NumMyDiagonallyDominant = 0;
  NumMyWeaklyDiagonallyDominant = 0;

  for( i=0 ; i<NumMyRows ; ++i ) {
    if( Diagonal[i]>SumOffDiagonal[i] ) 
      ++NumMyDiagonallyDominant;
    else if( Diagonal[i]==SumOffDiagonal[i] ) 
      ++NumMyWeaklyDiagonallyDominant;
    /* else nothing to track */
  }

  NumDiagonallyDominant = ML_Comm_GsumInt(Op->comm, NumMyDiagonallyDominant);
  NumWeaklyDiagonallyDominant = ML_Comm_GsumInt(Op->comm, NumMyWeaklyDiagonallyDominant);

  /* simply no output for MyPID>0, only proc 0 write on os */
  if( MyPID == 0 ) {

    if( name != NULL ) 
      printf("\n\t*** Analysis of ML_Operator `%s' ***\n\n", name);
    else
      printf("\n\t*** Analysis of ML_Operator ***\n\n");
    printf("\t%-50s = %d\n", 
	   "Number of global rows", NumGlobalRows);
    printf("\t%-50s = %d\n",
	   "Number of equations", Op->num_PDEs);
    printf("\t%-50s = %d\n",
	   "Number of stored elements", TotalNumGivenNonzeros);
    printf("\t%-50s = %d\n",
	   "Number of nonzero elements", TotalNumNonzeros);
    printf("\t%-50s = %d\n",
	   "Mininum number of nonzero elements/row", MinNumNonzeros);
    printf("\t%-50s = %d\n",
	   "Maximum number of nonzero elements/row", MaxNumNonzeros);
    printf("\t%-50s = %f\n",
	   "Average number of nonzero elements/rows", 
	   1.0 * TotalNumNonzeros / NumGlobalRows);
    printf("\t%-50s = %d\n",
	   "Nonzero elements in strict lower part", NumLowerNonzeros);
    printf("\t%-50s = %d\n",
	   "Nonzero elements in strict upper part", NumUpperNonzeros);
    printf("\t%-50s = %d\n",
	   "Max |i-j|, a(i,j) != 0",MaxBandwidth);
    printf("\t%-50s = %d (= %5.2f%%)\n",
	   "Number of diagonally dominant rows",
	   NumDiagonallyDominant,
	   100.0*NumDiagonallyDominant/NumGlobalRows); 
    printf("\t%-50s = %d (= %5.2f%%)\n",
	   "Number of weakly diagonally dominant rows",
	   NumWeaklyDiagonallyDominant,
	   100.0*NumWeaklyDiagonallyDominant/NumGlobalRows);
    printf("\t%-50s = %d (= %5.2f%%)\n",
	   "Number of Dirichlet rows",
	   NumDirichletRows,
	   100.0*NumDirichletRows/NumGlobalRows);
    printf("\t%-50s = %f\n",
	   "||A||_F",sqrt(FrobeniusNorm));
    printf("\t%-50s = %f\n",
	   "Min_{i,j} ( a(i,j) )", MinElement);
    printf("\t%-50s = %f\n",
	   "Max_{i,j} ( a(i,j) )", MaxElement);
    printf("\t%-50s = %f\n",
	   "Min_{i,j} ( abs(a(i,j)) )", MinAbsElement);
    printf("\t%-50s = %f\n",
	   "Max_{i,j} ( abs(a(i,j)) )", MaxAbsElement);

  }

  /* Analyze elements on diagonal for the entire matrix */

  MyMin = DBL_MAX, MyMax = 0.0;
  for( i=0 ; i<NumMyRows ; ++i ) {
    if( Diagonal[i]<MyMin ) MyMin = Diagonal[i];
    if( Diagonal[i]>MyMax ) MyMax = Diagonal[i];
  }
  Min = - ML_Comm_GmaxDouble(Op->comm, -MyMin);
  Max =   ML_Comm_GmaxDouble(Op->comm, MyMax);

  if( MyPID == 0 ) {

    printf("\t%-50s = %f\n",
	   "Min_i ( abs(a(i,i)) )", Min);
    printf("\t%-50s = %f\n",
	   "Max_i ( abs(a(i,i)) )", Max);

  }

  /* Analyze elements off diagonal for the entire matrix */

  MyMin = DBL_MAX, MyMax = 0.0;
  for( i=0 ; i<NumMyRows ; ++i ) {
    if( SumOffDiagonal[i]<MyMin ) MyMin = SumOffDiagonal[i];
    if( SumOffDiagonal[i]>MyMax ) MyMax = SumOffDiagonal[i];
  }
  Min = - ML_Comm_GmaxDouble(Op->comm, -MyMin);
  Max =   ML_Comm_GmaxDouble(Op->comm, MyMax);

  if( MyPID == 0 ) {

    printf("\t%-50s = %f\n",
	   "Min_i ( \\sum_{j!=i} abs(a(i,j)) )", Min);
    printf("\t%-50s = %f\n",
	   "Max_i ( \\sum_{j!=i} abs(a(i,j)) )", Max);

  }

  /* cycle over all equations and analyze diagonal elements. 
   * This may show that the matrix is badly scaled */
  
  if (Op->num_PDEs > 1) {

    for( Equation=0 ; Equation<Op->num_PDEs ; ++Equation ) {

      /* Analyze elements on diagonal */

      MyMin = DBL_MAX, MyMax = 0.0;
      for( i=Equation ; i<NumMyRows ; i+=Op->num_PDEs ) {
	if( Diagonal[i]<MyMin ) MyMin = Diagonal[i];
	if( Diagonal[i]>MyMax ) MyMax = Diagonal[i];
      }
      Min = - ML_Comm_GmaxDouble(Op->comm, -MyMin);
      Max =   ML_Comm_GmaxDouble(Op->comm, MyMax);

      if( MyPID == 0 ) {

	printf("\t(Eq %d) %-50s = %f\n",
	       Equation,
	       "Min_i ( abs(a(i,i)) )", 
	       Min);
	printf("\t(Eq %d) %-50s = %f\n",
	       Equation,
	       "Max_i ( abs(a(i,i)) )", 
	       Max);

      }
    }
  } /* if (Op->num_PDEs > 1) */

  /* free memory */

  if (colInd != NULL) ML_free(colInd);
  if (colVal != NULL) ML_free(colVal);
  if (Diagonal != NULL) ML_free(Diagonal);
  if (SumOffDiagonal != NULL) ML_free(SumOffDiagonal);
  if (IsDiagonallyDominant != NULL) ML_free(IsDiagonallyDominant);
  if (global_id != NULL) ML_free(global_id);

  data.ML_id = ML_ID_KRYLOVDATA;
  data.ML_matrix = Op;
  data.ML_com = Op->comm;
  data.ML_print_freq = 1;
  ML_Power_ComputeEigenvalues(&data,NumMyRows,0);

  if( MyPID == 0 ) 
    printf("\t%-50s = %f\n",
	   "max eig(A) (using power method)",
	   data.ML_eigen_max);
	   
  ML_Power_ComputeEigenvalues(&data,NumMyRows,1);

  if( MyPID == 0 ) 
    printf("\t%-50s = %f\n",
	   "max eig(D^{-1}A) (using power method)",
	   data.ML_eigen_max);

  if( MyPID == 0 )
    printf("\n\tTotal time for analysis = %e (s)\n", GetClock() - time);
    
  return 0;

}

/****************************************************************
 * Largely inspired from Yousef Saad's SPARSKIT plot function.  *
 * Plots the sparsity pattern of an ML_Operator into a PS file. *
 * Op : ML_Operator pointer
 * title : char string containing the title. Can be NULL.
 * PrintDecomposition : if ML_TRUE, prints lines corresponding to
 *                      the row and column decomposition across
 *                      the processors. If ML_FALSE, do nothing.
 * NumPDEEqns : number of PDE equations. The function will plot
 *              the block structure of the matrix if NumPDEEqns > 1
 ****************************************************************/

int ML_Operator_PrintSparsity(ML_Operator* Op, char* title,
			      int PrintDecomposition,
			      int NumPDEEqns)
{

  int i;
  int NumMyRows;
  int m,nc,nr,maxdim,ltit;
  double lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz;
  double xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx,xx,yy;
  int square = ML_FALSE;
  /*change square to .true. if you prefer a square frame around
    a rectangular matrix */
  double haf = 0.5, zero = 0.0, conv = 2.54;
  char munt = 'E'; /* put 'E' for centimeters, 'U' for inches */
  int ptitle = 0; /* position of the title, 0 under the drawing,
		     else above */
  FILE* fp = NULL;
  int* colInd = NULL;
  double* colVal = NULL;
  int ierr;
  int j;
  int NumNonzeros;
  int allocated = 0;
  int NumMyCols;
  int NumGlobalRows;
  int NumGlobalCols;
  int Nghost;
  double* global_col_id = NULL;
  int row_offset;
  int col_offset;
  int pid;
  int MyPID;
  int NumProc;
  int grow;
  int gcol;
  int nlines;
  int* row_lines = NULL;
  int* col_lines = NULL;
  int* itmp = NULL;
  int row_isep;
  int col_isep;

  /* --------------------- execution begins ---------------------- */

  MyPID = Op->comm->ML_mypid;
  NumProc = Op->comm->ML_nprocs;

  NumMyRows = Op->outvec_leng;
  NumMyCols = Op->invec_leng;

  NumGlobalRows = ML_Comm_GsumInt(Op->comm,NumMyRows);
  NumGlobalCols = ML_Comm_GsumInt(Op->comm,NumMyCols);

  /* need to define the global column number */

  if (Op->getrow->pre_comm == NULL) 
    Nghost = 0;
  else {
    if (Op->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Op->getrow->pre_comm);
    Nghost = Op->getrow->pre_comm->total_rcv_length;
  }

  global_col_id = (double*) ML_allocate(sizeof(double) * (NumMyCols + Nghost));
  if (global_col_id == NULL) {
    fprintf(stderr,
	    "Not enough memory\n");
    exit(EXIT_FAILURE);
  }

  row_offset = ML_gpartialsum_int(NumMyRows, Op->comm);
  col_offset = ML_gpartialsum_int(NumMyCols, Op->comm);

  for (i = 0 ; i < NumMyCols ; ++i) {
    global_col_id[i] = (double) (col_offset + i);
  }

  ML_exchange_bdry(global_col_id,Op->getrow->pre_comm, 
		   Op->invec_leng,Op->comm,ML_OVERWRITE,NULL);

  /* to be changed for rect matrices */
  maxdim = (NumGlobalRows>NumGlobalCols)?NumGlobalRows:NumGlobalCols;
  maxdim /= NumPDEEqns;

  m = 1 + maxdim;
  nr = NumGlobalCols / NumPDEEqns + 1;
  nc = NumGlobalRows / NumPDEEqns + 1;

  if (munt == 'E') {
    u2dot = 72.0/conv;
    paperx = 21.0;
    siz = 10.0;
  }
  else {
    u2dot = 72.0;
    paperx = 8.5*conv;
    siz = siz*conv;
  }

  /* left and right margins (drawing is centered) */

  lrmrgn = (paperx-siz)/2.0;

  /* bottom margin : 2 cm */

  botmrgn = 2.0;
  /* c scaling factor */
  scfct = siz*u2dot/m;
  /* matrix frame line witdh */
  frlw = 0.25;
  /* font size for title (cm) */
  fnstit = 0.5;
  if (title) ltit = strlen(title);
  else       ltit = 0;
  /* position of title : centered horizontally */
  /*                     at 1.0 cm vertically over the drawing */
  ytitof = 1.0;
  xtit = paperx/2.0;
  ytit = botmrgn+siz*nr/m + ytitof;
  /* almost exact bounding box */
  xl = lrmrgn*u2dot - scfct*frlw/2;
  xr = (lrmrgn+siz)*u2dot + scfct*frlw/2;
  yb = botmrgn*u2dot - scfct*frlw/2;
  yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2;
  if (ltit == 0) {
    yt = yt + (ytitof+fnstit*0.70)*u2dot;
  } 
  /* add some room to bounding box */
  delt = 10.0;
  xl = xl-delt;
  xr = xr+delt;
  yb = yb-delt;
  yt = yt+delt;

  /* correction for title under the drawing */
  if ((ptitle == 0) && (ltit == 0)) {
    ytit = botmrgn + fnstit*0.3;
    botmrgn = botmrgn + ytitof + fnstit*0.7;
  }

  /* detect how to write the lines for process decomposition. */

  if (PrintDecomposition == ML_TRUE) {

    row_lines = (int*) ML_allocate(sizeof(int)*NumProc);
    col_lines = (int*) ML_allocate(sizeof(int)*NumProc);
    itmp  = (int*) ML_allocate(sizeof(int)*NumProc);
    for (i = 0 ; i < NumProc ; ++i ) {
      row_lines[i] = 0;
      col_lines[i] = 0;
    }
    row_lines[MyPID] = NumMyRows;
    col_lines[MyPID] = NumMyCols;
    nlines = NumProc;

    ML_Comm_GsumVecInt(Op->comm, row_lines, itmp, NumProc);
    ML_Comm_GsumVecInt(Op->comm, col_lines, itmp, NumProc);
    for (i = 1 ; i < NumProc ; ++i ) {
      row_lines[i] += row_lines[i-1];
      col_lines[i] += col_lines[i-1];
    }
  }
  
  /* begin of output */

  if( MyPID == 0 ) {

    fp = fopen("Amat.ps","w");

    fprintf(fp,"%%!\n");
    fprintf(fp,"%%Creator: ML_Operator_PrintSparsity\n");
    fprintf(fp,"%%BoundingBox: %f %f %f %f\n",
	    xl,yb,xr,yt);
    fprintf(fp,"%%EndComments\n");
    fprintf(fp,"/cm {72 mul 2.54 div} def\n");
    fprintf(fp,"/mc {72 div 2.54 mul} def\n");
    fprintf(fp,"/pnum { 72 div 2.54 mul 20 string ");
    fprintf(fp,"cvs print ( ) print} def\n");
    fprintf(fp,"/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def\n");

    /* we leave margins etc. in cm so it is easy to modify them if
       needed by editing the output file */
    fprintf(fp,"gsave\n");
    if (ltit != 0) {
      fprintf(fp,"/Helvetica findfont %e cm scalefont setfont\n",
	      fnstit);
      fprintf(fp,"%f cm %f cm moveto\n",
	      xtit,ytit);
      fprintf(fp,"(%s) Cshow\n", title);
      fprintf(fp,"%f cm %f cm translate\n",
	      lrmrgn,botmrgn);
    }
    fprintf(fp,"%f cm %d div dup scale \n",
	    siz,m);
    /* draw a frame around the matrix */

    fprintf(fp,"%f setlinewidth\n",
	    frlw);
    fprintf(fp,"newpath\n");
    fprintf(fp,"0 0 moveto ");
    if (square == ML_TRUE) {
      printf("------------------- %d\n",m);
      fprintf(fp,"%d %d lineto\n",
	      m,0);
      fprintf(fp,"%d %d lineto\n",
	      m, m);
      fprintf(fp,"%d %d lineto\n",
	      0, m);
    } 
    else {
      fprintf(fp,"%d %d lineto\n",
	      nc, 0);
      fprintf(fp,"%d %d lineto\n",
	      nc, nr);
      fprintf(fp,"%d %d lineto\n",
	      0, nr);
    }
    fprintf(fp,"closepath stroke\n");

    /*drawing the separation lines  */

    if (PrintDecomposition == ML_TRUE) {

      fprintf(fp,"0.2 setlinewidth\n");

      /* this can be used to draw the decomposition
       * into processors */
      for (i = 0 ; i < nlines - 1 ; ++i) {

	row_isep = row_lines[i];
	col_isep = col_lines[i];

	/*     horizontal lines  */

	yy =  1.0 * (NumGlobalCols / NumPDEEqns - col_isep) + haf;
	xx = 1.0 * (NumGlobalRows / NumPDEEqns + 1);
	printf("%e %e\n", xx, yy);
	fprintf(fp,"%f %f moveto \n",
		zero, yy);
	fprintf(fp,"%f %f lineto stroke\n",
		xx, yy);

	/*vertical lines  */

	xx = 1.0 * (row_isep) + haf;
	yy = 1.0 * (NumGlobalCols+1);
	fprintf(fp,"%f %f moveto \n",
		xx, zero);
	fprintf(fp,"%f %f lineto stroke\n",
		xx, yy);
      }
    }

    /* plotting loop */

    fprintf(fp,"1 1 translate\n");
    fprintf(fp,"0.8 setlinewidth\n");
    fprintf(fp,"/p {moveto 0 -.40 rmoveto \n");
    fprintf(fp,"           0  .80 rlineto stroke} def\n");

    fclose(fp);
  }
  
  if (PrintDecomposition == ML_TRUE) {
    ML_free(row_lines);
    ML_free(col_lines);
    ML_free(itmp);
  }

  allocated = 10;

  colInd = (int*)   ML_allocate(sizeof(int)*allocated);
  colVal = (double*)ML_allocate(sizeof(double)*allocated);

  for (pid = 0 ; pid < NumProc ; ++pid) {

    if (pid == MyPID) {

      fp = fopen("Amat.ps","a");
      if( fp == NULL ) {
	fprintf(stderr,"ERROR\n");
	exit(EXIT_FAILURE);
      }

      for (i = 0 ; i < NumMyRows ; ++i) {

	if (i % NumPDEEqns) continue;

	ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NumNonzeros);

	if( ierr == 0 ) {
	  do {
	    ML_free(colInd);
	    ML_free(colVal);
	    allocated *= 2;
	    colInd = (int*)   ML_allocate(sizeof(int)*allocated);
	    colVal = (double*)ML_allocate(sizeof(double)*allocated);

	    ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NumNonzeros);
	  } while( ierr == 0 );
	}

	grow = i + row_offset;
	for (j = 0 ; j < NumNonzeros ; ++j) {
	  if (colInd[j] % NumPDEEqns == 0) {
	    gcol = (int)(global_col_id[colInd[j]]);
	    grow /= NumPDEEqns;
	    gcol /= NumPDEEqns;
	    fprintf(fp,"%d %d p\n",
		    grow, NumGlobalCols - gcol - 1); 
	  }
	}
      }

      fprintf(fp,"%%end of data for this process\n");

      if( pid == NumProc - 1 )
	fprintf(fp,"showpage\n");

      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  ML_free(colInd);
  ML_free(colVal);
  ML_free(global_col_id);
  
  return 1;
}

/****************************************************************
 * Compute eigenvectors and eigenvalues of an ML_Operator.
 * This function converts the ML_Operator into a dense matrix,
 * then call LAPACK (DGEEV) to evaluate the eigenvalues (real and imaginary
 * part) and the right eigenvectors. (It is trivial to change
 * the function so that it computed the left eigenvectors.)
 *
 * Warning: the code may require a *lot* of memory (as the
 * ML_Operator must be stored as a dense FORTRAN matrix), and
 * may take a *lot* of time.
 *
 * parameters:
 * ===========
 * 
 * Amat :      ML_Operator to be analyzed.
 * Er :        double vector, or size Amat->outvec_leng.
 *             In output, it will contain the real part of the
 *             computed eigenvalues.
 * Ei :        double vector, of size Amat->outvec_leng.
 *             In output, it will contain the imaginary part of the
 *             computed eigenvalues.
 * vectors :   double vector, of size (Amat->outvec_leng) ^ 2.
 *             In output, it will contain the eigenvectors.
 *             See the LAPACK documentation for DGEEV for 
 *             details.
 *
 * return value:  
 * ============
 *
 * function returns 0 if ok, a negative value otherwise
 *
 * date:       19-Aug-03, MS
 ****************************************************************/
#include "ml_lapack.h"

int ML_Operator_Eigensolver_Dense(ML_Operator* Amat,
				  double* Er,
				  double* Ei,
				  double* vectors)
{

  int i;
  int j;
  int n;
  char jobvl, jobvr;
  int lda,  ldvl, ldvr, lwork, info;
  double *a, *work;
  int *ipiv;
  double time;
  int allocated = 1;
  int* colInd;
  double* colVal;
  int ierr;
  int ncnt;

  /* -------------- execution begins ---------------------------- */

  time = GetClock();

  if (Amat->comm->ML_nprocs != 1) {
    fprintf(stderr,
	    "*ML*ERR* Function ML_Operator_EigenSolver_Dense()'\n"
	    "*ML*ERR* can be called with 1 process only.\n");
    return(-1);
  }

  if (Amat->invec_leng != Amat->outvec_leng) {
    fprintf(stderr,
	    "*ML*ERR* Function `ML_Operator_EigenSolver_Dense()'\n"
	    "*ML*ERR* requires square matrices.\n");
    return(-2);
  }

  n = Amat->invec_leng;
  jobvl = 'N'; /* V/N to calculate/not calculate left eigenvectors
                  of matrix H.*/

  jobvr = 'V'; /* As above, but for right eigenvectors. */

  lda = n; /* The leading dimension of matrix a. */

#ifdef NOT_DEF
  /* crappy check, I got some problems when too much memory
   * is allocated. Better to return now than never.... */

  if (n * n * 8 > 536870912) {
    fprintf(stderr,
	    "ML*ERR* LAPACK analysis of finest matrix would\n"
	    "require %d Kbytes. This seems too\n"
	    "much to me. Now I return; maybe you can change the\n"
	    "source code (file %s, line %d)\n",
	    n * n * 8 /1024,
	    __FILE__,
	    __LINE__);
    return(-3);
  }
#endif

  a = (double*) ML_allocate(sizeof(double)*n*n);
  if( a == 0 ) {
    fprintf(stderr,
	    "*ML*ERR* not enough memory to allocate %d bytes\n"
	    "*ML*ERR* (file %s, line %d)\n",
	    (int)sizeof(double)*n*n,
	    __FILE__,
	    __LINE__);
    exit(EXIT_FAILURE);
  }
  
  /* set to zero the elements */
  for (i = 0; i < n; ++i)
   for (j = 0 ; j < n ; ++j ) 
     a[i+j*n] = 0.0;

  /* now insert the nonzero elements, row by row */
  
  colInd = (int*) ML_allocate(sizeof(int)*allocated);
  colVal = (double*) ML_allocate(sizeof(double)*allocated);
  if( (colInd == NULL) || (colVal == NULL) ) {
    fprintf(stderr,
	    "*ML*ERR* not enough memory for %d bytes\n"
	    "*ML*ERR* (file %s, line %d)\n",
	    allocated,
	    __FILE__,
	    __LINE__);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < n ; i++) {

    ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);

    if (ierr == 0) {
      do {
	ML_free(colInd);
	ML_free(colVal);
	allocated *= 2;
	colInd = (int*) ML_allocate(sizeof(int)*allocated);
	colVal = (double*) ML_allocate(sizeof(double)*allocated);
	if( (colInd == NULL) || (colVal == NULL) ) {
	  fprintf(stderr,
		  "*ML*ERR* not enough memory for %d bytes\n"
		  "*ML*ERR* (file %s, line %d)\n",
		  allocated,
		  __FILE__,
		  __LINE__);
	  exit(EXIT_FAILURE);
	}
	ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);
      } while (ierr == 0);
    }
    for (j = 0 ; j < ncnt ; ++j) {
      a[i+n*colInd[j]] = colVal[j];
    }
  }

  printf("Time to convert ML_Operator into LAPACK format = %e (s)\n",
	 GetClock() - time);

  time = GetClock();

  info = 0;

  ipiv = (int*) ML_allocate(sizeof(int) * n);
  ldvl = n;
  ldvr = n;
  work = (double*) ML_allocate(sizeof(double) * 4 * n);
  lwork = 4*n;

  for (i = 0 ; i < n ; i++ ) {
    Er[i] = 0.0, Ei[i] = 0.0;
  }

  /* largest and smallest eigenvalue (in module) of A */

  printf("Computing eigenvalues/eigenvectors of given ML_Operator\n"
	 "using LAPACK (DGEEV). This may take some time...\n");

  DGEEV_F77(&jobvl, &jobvr, &n, a, &n, Er, Ei, NULL,
            &ldvl, vectors, &ldvr, work, &lwork, &info);
 
  /* free memory */

  ML_free(colInd);
  ML_free(colVal);
  ML_free(a);
  ML_free(ipiv);
  ML_free(work);

  printf("Time for eigensolver = %e (s)\n", 
	 GetClock() - time);

  return(0);

} /* ML_Operator_EigenSolver_Dense */
