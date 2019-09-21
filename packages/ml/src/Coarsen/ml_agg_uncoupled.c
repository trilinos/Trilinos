/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators                                */
/*  (uncoupled aggregation)                                                  */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : August, 2000                                              */
/* ************************************************************************* */
/* Local Functions :                                                         */
/*    ML_Aggregate_CoarsenUncoupled                                          */
/*    ML_Aggregate_CoarsenUncoupledCore                                      */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_viz_stats.h"
#include "ml_op_utils.h"
#ifdef MB_MODIF_QR
#include "ml_qr_fix.h"
#endif

/* ************************************************************************* */
/* local defines                                                             */
/* ------------------------------------------------------------------------- */

#define ML_AGGR_READY      -11
#define ML_AGGR_NOTSEL     -12
#define ML_AGGR_SELECTED   -13
#define ML_AGGR_SELECTED2  -14
#define ML_AGGR_BDRY       -15
#define ML_AGGR_MINRANK      1
#define ML_AGGR_MAXLINK      2

/* ************************************************************************* */
/* ************************************************************************* */
/*          Uncoupled subroutines                                            */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* construct the tentative prolongator (local)                               */
/* This function assumes that the block information are stored in the        */
/* ML_Aggregate structure.                                                   */
/*  phase 1 : relax on the new seed point as Vanek                           */
/*  phase 2 : assign the rest of the nodes to one of the existing            */
/*            aggregate (attach_scheme), if possible.                        */
/*  phase 3 : see if the un-aggregated nodes have enough neighbors           */
/*            (min_nodes_per_aggregate) to form its own aggregate            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_CoarsenUncoupled(ML_Aggregate *ml_ag,
           ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     mypid, Nrows, nvblocks, *vblock_info = NULL, *vblock_info2 = NULL;
   int     i, j, k, m, m2, nullspace_dim, *col_ind, aggr_count, nvblockflag;
   int     nbytes, Ncoarse, *mat_indx=NULL,*aggr_index,nz_cnt, diff_level;
   int     *new_ia = NULL, *new_ja = NULL, maxnnz_per_row=500, minnnz_per_row=1e6;
   double  printflag;
   int     *amal_mat_indx=NULL, amal_count, **rows_in_aggs = NULL,ibeg;
   int     lwork, *agg_sizes = NULL, row, level, new_cnt, max_agg_size;
#ifndef MB_REMOVE
   int     offset, jnode, *agg_sizes_cum=NULL, index, info;
#else
   int     offset, jnode, index, info;
#endif
#ifdef MB_MODIF_QR
   int     numDeadNod; /* number of nodes with dead dofs on current coarse lev */
   int     dead;
   int     nCDofTrunc;
#endif
   int     iend, bdry_blk, *bdry_array;
   double  epsilon, rowsum_threshold, *col_val, *tmp_vect = NULL;
   double  dcompare1, dcompare2, rowsum, *new_val=NULL, *diagonal=NULL;
   double  *nullspace_vect=NULL, *new_null=NULL, *work=NULL, *qr_tmp=NULL;
   double  /*largest,*/ thesign, dtemp;
   char    *col_entered;
   ML_Operator *Cmatrix;
   struct  ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm *aggr_comm;
   ML_GetrowFunc     *getrow_obj;
   int               (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void              *getrowdata;
   /*MS*/
   ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
   int * graph_decomposition;
   int NumBlockRows;
   /*ms*/
   char *true_bdry;
   int diag_idx;
#ifdef ML_ROWSUM_DEBUG
   int num_points_reset_for_rowsum=0;   
#endif
#ifdef ML_AGGR_INAGGR
   char fname[80];
   FILE *fp;
   static int level_count = 0;
#endif
   double t0=0,delta1=0,delta2=0,delta3=0,delta4=0;

   StartTimer(&t0);
   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */
   mypid          = comm->ML_mypid;
   epsilon        = ml_ag->threshold;
   nullspace_dim  = ml_ag->nullspace_dim;
   nullspace_vect = ml_ag->nullspace_vect;
   Nrows          = Amatrix->outvec_leng;
   nvblocks       = ml_ag->nvblocks;
   vblock_info    = ml_ag->vblock_info;
   printflag      = ml_ag->print_flag;

   /*MS*/
   NumBlockRows = Nrows / ml_ag->num_PDE_eqns;
   /*ms*/
   /* ============================================================= */
   /* check that this function is called properly.                  */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   rowsum_threshold = ml_ag->rowsum_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && printflag  < ML_Get_PrintLevel())
     {
       printf("ML_Aggregate_CoarsenUncoupled : current level = %d\n",
	      ml_ag->cur_level);
       printf("ML_Aggregate_CoarsenUncoupled : current eps = %e\n",
	      epsilon);
       if(rowsum_threshold > 0.0) {
          printf("ML_Aggregate_CoarsenUncoupled : current rowsum eps = %e\n",
                 rowsum_threshold);
       }
     }
   epsilon = epsilon * epsilon;

   /* ============================================================= */
   /* generate MSR matrix from incoming A matrix                    */
   /* ============================================================= */

   /* ------------------------------------------------------------- */
   /* first find out whether the getrow function is available       */
   /* ------------------------------------------------------------- */

   getrow_obj = Amatrix->getrow;
   getrowfunc=getrow_obj->func_ptr;
   getrowdata = Amatrix;
   if ( getrowfunc == NULL )
     {
       printf("ML_Aggregate_CoarsenUncoupled ERROR : no getrow.\n");
       exit(-1);
     }

   /* ------------------------------------------------------------- */
   /* allocate initial temporary storage space for getrow           */
   /* also allocate space for storing the diagonal (if epsilon>0)   */
   /* ------------------------------------------------------------- */
   if (Amatrix->max_nz_per_row > maxnnz_per_row)
      maxnnz_per_row = Amatrix->max_nz_per_row;
   if (Amatrix->min_nz_per_row < minnnz_per_row)
      minnnz_per_row = Amatrix->min_nz_per_row;


   col_ind = (int *)    ML_allocate( maxnnz_per_row * sizeof(int) );
   col_val = (double *) ML_allocate( maxnnz_per_row * sizeof(double) );

   /* ------------------------------------------------------------- */
   /* find out about how much memory to allocate for the matrix     */
   /* ------------------------------------------------------------- */

   ML_Operator_Get_Diag(Amatrix,Amatrix->outvec_leng,&diagonal);
   nz_cnt = ML_Operator_ComputeNumNzs(Amatrix);

   /* ------------------------------------------------------------- */
   /* allocate memory for the entire matrix (only the column indices*/
   /* are needed since the matrix will be pruned here               */
   /* ------------------------------------------------------------- */

   nbytes = (nz_cnt + 1) * sizeof( int );
   ML_memory_alloc((void**) &mat_indx, (unsigned int) nbytes, "AVA");
   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, nz_cnt);

   if ( mypid == 0 && printflag  < ML_Get_PrintLevel())
     printf("Aggregation(UVB) : Total nonzeros = %d (Nrows=%d)\n",m,k);


   m2 = ML_Comm_GsumInt( comm, Amatrix->N_nonzeros);
   if ( ml_ag->operator_complexity == 0.0 ){
     ml_ag->fine_complexity = 1.0 * m2;
     ml_ag->operator_complexity = 1.0 * m2;
   }
   else ml_ag->operator_complexity += 1.0 * m2;

   /* ------------------------------------------------------------- */
   /* extract the matrix using the getrow function                  */
   /* (pruning is done at this stage)                               */
   /* ------------------------------------------------------------- */

   nz_cnt = Nrows + 1;
   mat_indx[0] = nz_cnt;

   for ( i = 0; i < Nrows; i++ )
   {
     int itmp = nz_cnt;
     int allocated_space = maxnnz_per_row;
     ML_get_matrix_row( (ML_Operator *)getrowdata, 1, &i, &allocated_space, &col_ind,
                        &col_val, &m, 0);
     if ( m > maxnnz_per_row ) Amatrix->max_nz_per_row = m;
     if ( m < minnnz_per_row && m>0 ) Amatrix->min_nz_per_row = m;

     rowsum = 0;
     diag_idx=0;
     for (j = 0; j < m; j++)
       {
         jnode = col_ind[j];
         rowsum += col_val[j];
         if(jnode ==i) diag_idx=j;

         if ( jnode != i && jnode < Nrows && epsilon > 0.0 )
           {             
             dcompare1 = col_val[j] * col_val[j];
             if ( dcompare1 > 0.0 )
               {
                 dcompare2 = diagonal[i] * diagonal[jnode];
                 dcompare1 = ML_dabs( dcompare1 );
                 dcompare2 = ML_dabs( dcompare2 );
                 if ( dcompare1 >= epsilon * dcompare2 )
                   mat_indx[nz_cnt++] = col_ind[j];
               }
           }
         else if ( jnode != i && jnode < Nrows && col_val[j] != 0.0)
           mat_indx[nz_cnt++] = col_ind[j];
       }

     /* Reaction-diffusion dropping.  If the rowsum is sufficiently than the diagonal, then
        we should be in a reaction-limited regime at this node and can afford to drpp *all* 
        connections, effectively turning this guy into a Dirichlet unknown.  We trust that
        the smoother is sufficient for these unknowns - CMS 7/23/19 */
     if(rowsum_threshold > 0.0 && fabs(rowsum) > fabs(diagonal[i]) * rowsum_threshold) {
#ifdef ML_ROWSUM_DEBUG
       num_points_reset_for_rowsum++;
#endif
       mat_indx[itmp]=col_ind[diag_idx];
       nz_cnt = itmp;
     }

     mat_indx[i+1] = nz_cnt;
     /* JJH 6/28/06
     printf("(%d) nz_cnt - itmp = %d, m = %d\n",i,nz_cnt-itmp,m);
     */
     if ( nz_cnt-itmp ==0 ) mat_indx[i] = -mat_indx[i];
   }
   ML_free(col_ind);
   ML_free(col_val);
   StopTimer(&t0,&delta1);


   /* ============================================================= */
   /* Construct the matrix that relates to the nodes by combining   */
   /* the rows of the matrix corresponding to different PDE's at    */
   /* the same node.                                                */
   /* ============================================================= */

   StartTimer(&t0);
   nvblockflag = 0;
   if ( nvblocks == 0 )
     {
       /* MB_MODIF */
       nvblocks = Nrows/ml_ag->num_PDE_eqns;
       nbytes   = nvblocks * sizeof(int);
       ML_memory_alloc((void**) &vblock_info, (unsigned int) nbytes,"AVE");
       /* MB_MODIF */
       for ( i = 0; i < nvblocks; i++ ) vblock_info[i] = ml_ag->num_PDE_eqns;
       nvblockflag = 1;
     }
   nbytes = (nvblocks + 1)* sizeof(int);
   if (nbytes > 0) ML_memory_alloc((void**) &vblock_info2, (unsigned int) nbytes,"AVC");
   vblock_info2[0] = vblock_info[0];
   for ( i = 1; i < nvblocks; i++ )
     vblock_info2[i] = vblock_info2[i-1] + vblock_info[i];

   bdry_array = (int *) ML_allocate(sizeof(int)*nvblocks);
   for ( i = 0; i < nvblocks; i++) bdry_array[i] = 0;


   if ((nvblockflag == 1) && (ml_ag->num_PDE_eqns == 1))
   {
     for ( i = 0; i < nvblocks; i++) {
       if ( mat_indx[i] < 0 ) bdry_array[i] = 1;
     }
   }
   else
   {
     nbytes = (nz_cnt + 1) * sizeof( int ); /* probably excessive */
     if (nbytes > 0) ML_memory_alloc((void**) &amal_mat_indx, (unsigned int) nbytes,"AVB");

     amal_count = nvblocks + 1;
     amal_mat_indx[0] = amal_count;
     row = 0;

     col_entered = (char *) ML_allocate(sizeof(char)*(1+ nvblocks) );
     if (col_entered == NULL) {
       printf("Not enough space in ML_aggregate\n");
       exit(1);
     }
     for ( i = 0; i < nvblocks; i++) col_entered[i] = 'F';
     for ( i = 0; i < nvblocks; i++)
       {
	 col_entered[i] = 'T';
	 bdry_blk = 1;
	 for ( j = 0; j < vblock_info[i]; j++)
	   {
	     if ( mat_indx[row+j] >= 0 ) bdry_blk = 0;
	     else                        mat_indx[row+j] = - mat_indx[row+j];
	   }
	 if ( bdry_blk == 1 ) {row += vblock_info[i]; bdry_array[i] = 1;}
	 else
     {
       for ( j = 0; j < vblock_info[i]; j++)
       {
		 ibeg = mat_indx[row];
		 iend = mat_indx[row+1];
		 if ( iend < 0 ) iend = - iend;
		 for ( k = ibeg; k < iend; k++)
         {
           if ( mat_indx[k] < vblock_info2[0] ) index = 0;
           else
           {
			 /* if we really have variable blocks we search
			    though ... this might be better:
			    index=ML_fastsorted_search(mat_indx[k],nvblocks,
			          vblock_info2,
			          mat_indx[k]/ml_ag->num_PDE_eqns-1);
			 */
			 if (nvblockflag == 1)
			   index = mat_indx[k]/ml_ag->num_PDE_eqns-1;
			 else index=ML_sorted_search(mat_indx[k],nvblocks,vblock_info2);

			 if ( index < 0 ) index = - index;
			 else             index++;
           }
		   if ( index < 0 || index >= nvblocks )
             printf("ERROR : in almalgamation %d => %d(%d).  Have you specified the correct number of DOFs?\n",mat_indx[k],
			      index,nvblocks);
           if (col_entered[index] == 'F')
           {
			 amal_mat_indx[ amal_count++] = index;
			 col_entered[index] = 'T';
           }
         }
         row++;
       }
     }
	 amal_mat_indx[i+1] = amal_count;
	 col_entered[i] = 'F';
	 for ( j = amal_mat_indx[i]; j < amal_mat_indx[i+1]; j++)
	   col_entered[ amal_mat_indx[j]] = 'F';
       }
     ML_free(col_entered);

     if ( mypid == 0 && printflag  < ML_Get_PrintLevel())
       printf("Aggregation(UVB) : Amalgamated matrix done \n");
   }
#ifdef ML_ROWSUM_DEBUG
   printf("ML_ROWSUM_DEBUG: # points reset via rowsum = %d\n",num_points_reset_for_rowsum);
#endif

   StopTimer(&t0,&delta2);

   /* ============================================================= */
   /* perform coarsening                                            */
   /* ============================================================= */

   StartTimer(&t0);
   csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(
						  struct ML_CSR_MSRdata));
   csr_data->values  = NULL;
   Cmatrix = ML_Operator_Create(Amatrix->comm);
   if ((nvblockflag == 1) && (ml_ag->num_PDE_eqns == 1)) {
     /*JJH Prevents aggregation error in Phase2_3 clean up.*/
     for ( i = 0; i < nz_cnt; i++ ) mat_indx[i] = abs(mat_indx[i]);

     csr_data->columns = mat_indx;
     ML_Operator_Set_ApplyFuncData(Cmatrix,Amatrix->invec_leng,
				   Amatrix->outvec_leng,
				   csr_data,
				   Amatrix->outvec_leng,NULL,0);
     ML_Operator_Set_Getrow(Cmatrix, Cmatrix->outvec_leng,
			    MSR_get_ones_rows);
     true_bdry = ML_Operator_IdentifyDirichletRows(Cmatrix);
     ML_Aggregate_CoarsenUncoupledCore(ml_ag,comm,Cmatrix,mat_indx,
				       bdry_array, &aggr_count, &aggr_index, true_bdry);
   }
   else {
     /*JJH same error as above could occur here!*/
     csr_data->columns = amal_mat_indx;
     ML_Operator_Set_ApplyFuncData(Cmatrix,nvblocks, nvblocks,
				   csr_data, nvblocks, NULL,0);
     ML_Operator_Set_Getrow(Cmatrix, nvblocks, MSR_get_ones_rows);
     true_bdry = ML_Operator_IdentifyDirichletRows(Cmatrix);
     ML_Aggregate_CoarsenUncoupledCore(ml_ag,comm,Cmatrix,amal_mat_indx,
				       bdry_array, &aggr_count, &aggr_index, true_bdry);
#ifdef ML_AGGR_INAGGR

   for ( i = 0; i < nvblocks; i++ ) aggr_index[i] = -1;
   sprintf(fname,"agg0_%d",level_count); level_count++;
   fp = fopen(fname,"r");
   if (fp == NULL) {
      printf("Cannot open aggregate file %s for reading.\n",fname);
      exit(1);
   }
   else {
      printf("Reading aggregate for level %d from file %s\n",
             level_count-1,fname);
      fflush(stdout);
   }
   aggr_count = 0;
   for (i = 0; i <nvblocks; i++) {
     if ( fscanf(fp,"%d%d",&j,&k) != 2) break;
      aggr_index[j/5] = k;
      if (k >= aggr_count) aggr_count = k+1;
      fscanf(fp,"%d%d",&j,&k);
      fscanf(fp,"%d%d",&j,&k);
      fscanf(fp,"%d%d",&j,&k);
      fscanf(fp,"%d%d",&j,&k);
   }
   fclose(fp);
   printf("Read in %d aggregates\n\n",aggr_count);

#endif
   }
   ML_Operator_Destroy(&Cmatrix);
   ML_free(csr_data);
   ML_free( bdry_array );

   /* ********************************************************************** */
   /* I allocate room to copy aggr_index and pass this value to the user,    */
   /* who will be able to analyze and visualize this after the construction  */
   /* of the levels. This way, the only price we have to pay for stats and   */
   /* viz is essentially a little bit of memory.                             */
   /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
   /* I set the pointers using the ML_Aggregate_Info structure. This is      */
   /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
   /* ********************************************************************** */

   if( Amatrix->to != NULL && Amatrix->to->Grid != NULL &&
       Amatrix->to->Grid->Grid != NULL ) {

     graph_decomposition = (int *)ML_allocate(sizeof(int)*(NumBlockRows+1));
     if( graph_decomposition == NULL ) {
       fprintf( stderr,
        "*ML*ERR* Not enough memory for %d bytes\n"
        "*ML*ERR* (file %s, line %d)\n",
        (int)sizeof(int)*NumBlockRows,
        __FILE__,
          __LINE__ );
       exit( EXIT_FAILURE );
     }

     for( i=0 ; i<NumBlockRows ; i++ ) graph_decomposition[i] = aggr_index[i];

     aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *)(Amatrix->to->Grid->Grid);
     aggr_viz_and_stats->graph_decomposition = graph_decomposition;
     aggr_viz_and_stats->Nlocal = NumBlockRows;
     aggr_viz_and_stats->Naggregates = aggr_count;
     aggr_viz_and_stats->local_or_global = ML_LOCAL_INDICES;
     aggr_viz_and_stats->is_filled = ML_YES;
     aggr_viz_and_stats->Amatrix = NULL;
   }

   StopTimer(&t0,&delta3);

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   StartTimer(&t0);

   Ncoarse = aggr_count * nullspace_dim;
   level   = ml_ag->cur_level;
   nbytes  = Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), (unsigned int) nbytes, "AVC");
   new_cnt = aggr_count;
   for ( i = 0; i < nvblocks; i++ )
   {
      if ( i == 0 ) offset = 0;
      else          offset = vblock_info2[i-1];
      if ( aggr_index[i] >= 0 )
      {
         for ( j = 0; j < vblock_info[i]; j++ )
            ml_ag->aggr_info[level][offset+j] = aggr_index[i];
      }
      else
      {
         for ( j = 0; j < vblock_info[i]; j++ )
            ml_ag->aggr_info[level][offset+j] = new_cnt;
         new_cnt++;
      }
   }
   /*ml_ag->aggr_count[level] = aggr_count; */
   ml_ag->aggr_count[level] = new_cnt; /* for relaxing boundary points */

   /* ------------------------------------------------------------- */
   /* Ray's debugging tools follows.                                */
   /* ------------------------------------------------------------- */

#ifdef ML_AGGR_TEST
   if ( ml_ag->cur_level == ml_ag->max_levels-1 && nprocs > 1 )
   {
      itmp_array  = (int *) ML_allocate( nprocs * sizeof(int));
      itmp_array2 = (int *) ML_allocate( nprocs * sizeof(int));
      for ( i = 0; i < nprocs; i++ ) itmp_array[i] = 0;
      itmp_array[mypid] = Nrows;
      ML_gsum_vec_int(&itmp_array, &itmp_array2, nprocs, comm);
      Nrows_offset = 0;
      for ( i = 0; i < mypid; i++ ) Nrows_offset += itmp_array[i];
      for ( i = 0; i < nprocs; i++ ) itmp_array[i] = 0;
      itmp_array[mypid] = count;
      ML_gsum_vec_int(&itmp_array, &itmp_array2, nprocs, comm);
      naggr_offset = 0;
      for ( i = 0; i < mypid; i++ ) naggr_offset += itmp_array[i];
      ML_free(itmp_array);
      ML_free(itmp_array2);
      sprintf(zfn, "mlaggr.out.%d",mypid);
      zfp = fopen(zfn,"w");
      fprintf(zfp, "aggr_count = %d \n", count);
      for ( i = 0; i < Nrows; i++ )
         fprintf(zfp, "%d \n", ml_ag->aggr_info[level][i]+naggr_offset);
      fclose(zfp);
      ML_gsum_scalar_int(&i, &j, comm);
      ML_gsum_scalar_int(&i, &j, comm);
      exit(1);
   }
   if ( ml_ag->cur_level == ml_ag->max_levels-1 && nprocs == 1 )
   {
      sprintf(zfn, "mlaggr.out");
      zfp = fopen(zfn,"r");
      for ( i = 0; i < Nrows; i++ )
         fscanf(zfp, "%d", &(ml_ag->aggr_info[level][i]));
      fclose(zfp);
      aggr_count = 0;
      for ( i = 0; i < Nrows; i++ )
         if ( ml_ag->aggr_info[level][i] > aggr_count )
            aggr_count = ml_ag->aggr_info[level][i];
      aggr_count++;
      Ncoarse = aggr_count * nullspace_dim;
      for ( i = 0; i < Nrows; i++ )
         aggr_index[i] = ml_ag->aggr_info[level][i];
      ml_ag->aggr_count[level] = aggr_count;
      printf("ML_Aggregate : total no. of aggregates input = %d\n",aggr_count);
   }
#endif

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new operator and null space  */
   /* ------------------------------------------------------------- */

   nbytes = ( Nrows + 1 ) * sizeof(int);
   ML_memory_alloc((void**)&(new_ia), (unsigned int) nbytes, "AVM");
   nbytes = Nrows * nullspace_dim * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), (unsigned int) nbytes, "AVN");
   nbytes = Nrows * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_val), (unsigned int) nbytes, "AVO");
   nbytes = Ncoarse * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null), (unsigned int) nbytes,"AVX");
   for (i = 0; i < Ncoarse*nullspace_dim; i++) new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each roll will have at most nullspace_dim nonzero entries)   */
   /* ------------------------------------------------------------- */

   for (i = 0; i <  Nrows*nullspace_dim; i++) new_ja[i] = 0;
   for (i = 0; i <  Nrows*nullspace_dim; i++) new_val[i]= 0.0;
   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* ------------------------------------------------------------- */
   /* temporary variables for use in subsequent processing          */
   /* ------------------------------------------------------------- */

   nbytes = aggr_count * sizeof(int);
   ML_memory_alloc((void**)&agg_sizes, (unsigned int) nbytes,"AVI");
#ifndef MB_REMOVE
/*-mb: this array is allocated and set, but never used! */
   ML_memory_alloc((void**)&agg_sizes_cum, (unsigned int) nbytes, "AVJ");
#endif

   /* ------------------------------------------------------------- */
   /* fill the temporary variables and also find the maximum        */
   /* aggregate size for allocation of qr_tmp                       */
   /* ------------------------------------------------------------- */

   for (i = 0; i < nvblocks; i++)
   {
      if (aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
         agg_sizes[aggr_index[i]] += vblock_info[i];
      else if (aggr_index[i] != -1)
      {
         printf("%d : CoarsenUncoupled - wrong index %d(%d)%d\n",mypid,
                      aggr_index[i], aggr_count,i);
         exit(1);
      }
   }
   max_agg_size = agg_sizes[0];
#ifndef MB_REMOVE
   if ( aggr_count > 0 ) agg_sizes_cum[0] = 0;
   for (i = 0; i < aggr_count-1; i++)
   {
      agg_sizes_cum[i+1] = agg_sizes_cum[i] + agg_sizes[i];
      if (agg_sizes[i+1] > max_agg_size) max_agg_size = agg_sizes[i+1];
   }
   ML_memory_free((void**)&agg_sizes_cum);
#else
   for (i = 0; i < aggr_count-1; i++)
   {
      if (agg_sizes[i+1] > max_agg_size) max_agg_size = agg_sizes[i+1];
   }
#endif

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLt");
   for (i = 0; i < aggr_count; i++)
      rows_in_aggs[i] = (int *) ML_allocate( (agg_sizes[i]+1)*sizeof(int) );
   if ( (aggr_count > 0) && (rows_in_aggs[aggr_count-1] == NULL) )
   {
      printf("Error: couldn't allocate memory in CoarsenUncoupledVB\n");
      exit(1);
   }
   for (i = 0; i < aggr_count; i++) agg_sizes[i] = 0;
   for (i = 0; i < nvblocks; i++)
   {
      if ( aggr_index[i] >= 0 )
      {
         for (j = 0; j < vblock_info[i]; j++)
         {
            index = agg_sizes[aggr_index[i]]++;
            if ( i == 0 ) offset = 0;
            else          offset = vblock_info2[i-1];
            rows_in_aggs[aggr_index[i]][index] = offset + j;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&qr_tmp, (unsigned int) nbytes, "AVU");
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&tmp_vect, (unsigned int) nbytes, "AVT");

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&work, (unsigned int) nbytes, "AVK");

   work[0] = lwork;
#ifdef MB_MODIF_QR
   ML_qr_fix_Create(aggr_count, nullspace_dim); /* alloc array in structure */
   numDeadNod = 0;  /* number of nodes with dead dofs on current coarse lev */
#endif
   for (i = 0; i < aggr_count; i++)
   {
      /* set up the matrix we want to decompose into Q and R: */

      if (nullspace_vect == NULL)
      {
         for (j = 0; j < agg_sizes[i]; j++)
            for (k = 0; k < nullspace_dim; k++)
            {
               row = rows_in_aggs[i][j];
               if ( row < vblock_info2[0] ) index = 0;
               else
               {
		 /* if we really have variable blocks we search. This might
		    be better ....
                  index=ML_fastsorted_search(row,nvblocks,vblock_info2,
				     row/ml_ag->num_PDE_eqns-1);
		 */
		 if (nvblockflag == 1)
		   index = row/ml_ag->num_PDE_eqns-1;
		 else index=ML_sorted_search(row,nvblocks,vblock_info2);
                  if ( index < 0 ) index = - index;
                  else             index++;
               }
               if ( index == 0 ) offset = row;
               else              offset = row-vblock_info2[index-1];
               if ( offset == k ) qr_tmp[k*agg_sizes[i] + j] = 1.0;
               else               qr_tmp[k*agg_sizes[i] + j] = 0.0;
            }
      }
      else
      {
         for (k = 0; k < nullspace_dim; k++)
            for (j = 0; j < agg_sizes[i]; j++)
               qr_tmp[k*agg_sizes[i] + j] =
                  nullspace_vect[ k*Nrows + rows_in_aggs[i][j] ];
      }

#ifdef MB_MODIF_QR
     /* Graciously treat the case where we have more kernel
      * components than freedom. For this, we need the xCDeadNodDof structure.
      */
      nCDofTrunc = nullspace_dim;
      dead       = 0;
      if (nullspace_dim > agg_sizes[i]) {
         /* mark the coarse degree(s) of freedom that will be dead */
          nCDofTrunc = agg_sizes[i];
          dead = 1;
          ML_qr_fix_markDOFsAsDead(i,nCDofTrunc, nullspace_dim);
          numDeadNod++;
      }

      /* now calculate QR using an LAPACK routine */

      if ( nullspace_dim == 1 )
      {
         dtemp = 0.0;
         for (j = 0; j < agg_sizes[i]; j++)
            dtemp += ( qr_tmp[j] * qr_tmp[j] );
         dtemp = sqrt( dtemp );
         tmp_vect[0] = qr_tmp[0];
         qr_tmp[0] = dtemp;
      }
      else
      {
         DGEQRF_F77(&(agg_sizes[i]), &nCDofTrunc, qr_tmp,
                           &(agg_sizes[i]), tmp_vect, work, &lwork, &info);
         if (info != 0) {
            pr_error("ERROR (CoarsenUncoupled) : dgeqrf returned a non-zero\n");
         }
      }

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGk");
      }
/*-mb: never truncate here
      else lwork=(int) work[0];
*/

     /* the upper triangle of qr_tmp is now R, so copy that into the
      * new nullspace.
      * treat separately the case where any coarse dofs are dead */

      if (dead) {
          for (k = 0; k < nullspace_dim; k++) {
            if (ML_qr_fix_isDOFDead(i,k)) {
               for (j = 0; j < k+1; j++)
                  new_null[i*nullspace_dim+j+k*Ncoarse] = 0.e0;
            } else {
               for (j = 0; j < k+1; j++)
                  new_null[i*nullspace_dim+j+k*Ncoarse] =
                      qr_tmp[j+agg_sizes[i]*k];
            }
          }
      } else {
          for (k = 0; k < nullspace_dim; k++)
            for (j = 0; j < k+1; j++)
              new_null[i*nullspace_dim+j+k*Ncoarse] = qr_tmp[j+agg_sizes[i]*k];
      }

      /* to get this block of P, need to run qr_tmp through another LAPACK
         function: */

      if ( nullspace_dim == 1 )
      {
         dtemp = qr_tmp[0];
         qr_tmp[0] = tmp_vect[0];
         dtemp = 1.0 / dtemp;
         for (j = 0; j < agg_sizes[i]; j++)
            qr_tmp[j] *= dtemp;
      }
      else
      {
         DORGQR_F77(&(agg_sizes[i]), &nCDofTrunc, &nCDofTrunc,
                 qr_tmp, &(agg_sizes[i]), tmp_vect, work, &lwork, &info);
         if (info != 0)
            pr_error("ERROR (CoarsenUncoupled): dorgqr returned a non-zero\n");
         if (dead) {
           /* modify dead columns of Q if any */
           for (k = 0; k < nullspace_dim; k++) {
             if (ML_qr_fix_isDOFDead(i,k)) {
               for (j = 0; j < agg_sizes[i]; j++)
                 qr_tmp[k*agg_sizes[i] + j] = 0.e0;
             }
           }
         }
      }
     /* -mb: we could check cols of Q locally to see if more dofs are dead */

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AVM");
      }
/*-mb: never truncate here
      else lwork=(int) work[0];
*/
#else /*MB_MODIF_QR*/

      /* now calculate QR using an LAPACK routine */

      if ( nullspace_dim == 1 )
      {
         dtemp = 0.0;
         for (j = 0; j < agg_sizes[i]; j++)
            dtemp += ( qr_tmp[j] * qr_tmp[j] );
         dtemp = sqrt( dtemp );
         tmp_vect[0] = qr_tmp[0];
         qr_tmp[0] = dtemp;
      }
      else
      {
         DGEQRF_F77(&(agg_sizes[i]), &nullspace_dim, qr_tmp,
                           &(agg_sizes[i]), tmp_vect, work, &lwork, &info);
         if (info != 0)
            pr_error("ERROR (CoarsenUncoupled) : dgeqrf returned a non-zero\n");
      }

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGk");
      }
      else lwork=(int) work[0];

      /* the upper triangle of qr_tmp is now R, so copy that into the
         new nullspace */

      for (j = 0; j < nullspace_dim; j++)
         for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse] = qr_tmp[j+agg_sizes[i]*k];

      /* to get this block of P, need to run qr_tmp through another LAPACK
         function: */

      if ( nullspace_dim == 1 )
      {
         dtemp = qr_tmp[0];
         qr_tmp[0] = tmp_vect[0];
         dtemp = 1.0 / dtemp;
         for (j = 0; j < agg_sizes[i]; j++)
            qr_tmp[j] *= dtemp;
      }
      else
      {
         DORGQR_F77(&(agg_sizes[i]), &nullspace_dim, &nullspace_dim,
                 qr_tmp, &(agg_sizes[i]), tmp_vect, work, &lwork, &info);
         if (info != 0)
            pr_error("ERROR (CoarsenUncoupled): dorgqr returned a non-zero\n");
      }

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AVM");
      }
      else lwork=(int) work[0];
#endif/*MO_MODIF*/

      /* now copy Q over into the appropriate part of P: */
      /* The rows of P get calculated out of order, so I assume the Q is
         totally dense and use what I know of how big each Q will be to
         determine where in ia, ja, etc each nonzero in Q belongs.  If I
         did not assume this, I would have to keep all of P in memory in
         order to determine where each entry should go */

      for (j = 0; j < agg_sizes[i]; j++)
      {
         /*largest = 0.0;*/ thesign = 1.;
/* this is a bad bug. -mb.
         for (k = 0; k < nullspace_dim; k++)
         {
            if ( ML_dabs(qr_tmp[k*agg_sizes[i]+j]) > largest )
            {
               largest = ML_dabs(qr_tmp[ k*agg_sizes[i] + j ]);
               if ( qr_tmp[ k*agg_sizes[i] + j ] < 0.0) thesign = -1.;
               else thesign = 1.;
            }
         }
*/
         for (k = 0; k < nullspace_dim; k++)
         {
            index = new_ia[rows_in_aggs[i][j]] + k;
            new_ja [index] = i * nullspace_dim + k;
            new_val[index] = thesign * qr_tmp[ k*agg_sizes[i] + j ];
         }
      }
   }
#ifdef MB_MODIF_QR
  /* set the number of nodes with dead dofs on current coarse grid */
   ML_qr_fix_setNumDeadNod(numDeadNod);
   k = numDeadNod;
   ML_gsum_scalar_int(&k, &i, comm);
   if (mypid == 0 && printflag  < ML_Get_PrintLevel())
     printf("Aggregation(UC) : QR factorization - too small aggregates = %d\n",k);
#endif

   ML_Aggregate_Set_NullSpace(ml_ag, nullspace_dim, nullspace_dim,
                              new_null, Ncoarse);
   ML_memory_free( (void **) &new_null);

   /* ------------------------------------------------------------- */
   /* compress the prolongation operator                            */
   /* ------------------------------------------------------------- */

   k     = new_ia[0];
   index = k;
   nz_cnt = 0;
   for (i = 0; i < Nrows; i++)
   {
      for (j = k; j < new_ia[i+1]; j++ )
      {
         if ( new_val[j] != 0.0 )
         {
            new_val[index]  = new_val[j];
            new_ja[index++] = new_ja[j];
            nz_cnt++;
         }
      }
	  /* JJH This code fragment forces at least one entry in each row,
			 even if that entry is zero.  This can cause failures in
             parallel.
      if ( index == new_ia[i] )
      {
         new_val[index] = new_val[k]; new_ja[index++] = new_ja[k];
      }
	  -- JJH */
      k = new_ia[i+1];
      new_ia[i+1] = index;
   }
   ML_memory_alloc((void**) &csr_data,sizeof(struct ML_CSR_MSRdata),"AVP");

   (*Pmatrix)->N_nonzeros = ML_Comm_GsumInt( comm, nz_cnt);

   /* this must be set so that the hierarchy generation does not abort early
      in adaptive SA */
   (*Pmatrix)->num_PDEs = nullspace_dim;

   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   ML_Operator_Set_ApplyFuncData( *Pmatrix, Ncoarse, Nrows,
                                  csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm), "AVQ");
   aggr_comm->comm = comm;
   aggr_comm->N_send_neighbors = 0;
   aggr_comm->N_recv_neighbors = 0;
   aggr_comm->send_neighbors = NULL;
   aggr_comm->recv_neighbors = NULL;
   aggr_comm->send_leng = NULL;
   aggr_comm->recv_leng = NULL;
   aggr_comm->send_list = NULL;
   aggr_comm->local_nrows = Ncoarse;

   m = 0;
   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm),
                           ML_Aggregate_ExchangeBdry, aggr_comm,
                           comm, Ncoarse, m);
   ML_Operator_Set_Getrow((*Pmatrix), Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc((*Pmatrix), CSR_matvec);

   /* ============================================================= */
   /* clean up                                                      */
   /* ============================================================= */

   ML_memory_free((void**) &vblock_info2);
   ML_memory_free((void**) &mat_indx);
   if (amal_mat_indx != NULL)   ML_memory_free((void**) &amal_mat_indx);
   ML_memory_free((void**) &aggr_index);
   ML_memory_free((void**)&agg_sizes);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp); ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);
   ML_memory_free((void**)&aggr_comm);

   /* tuminaro change */
   /* I think this is what Charles wanted */
   if ( nvblockflag == 1 ) ML_memory_free((void**)&vblock_info);

   StopTimer(&t0,&delta4);
   if (ML_Get_PrintLevel() > 9) {
#    ifdef ML_TIMING
     if (mypid == 0)
       printf("Detailed timing for forming tentative prolongator (level %d)\n", ml_ag->cur_level);
#    endif
     ReportTimer(delta1,"Drop weak entries      ",comm);
     ReportTimer(delta2,"Amalgamate             ",comm);
     ReportTimer(delta3,"Coarsen                ",comm);
     ReportTimer(delta4,"QR and matrix formation",comm);
#    ifdef ML_TIMING
     if (mypid == 0)
       printf("\n");
#    endif
   }

   return Ncoarse;
}

/* ************************************************************************* */
/* construct the tentative prolongator (local)                               */
/*  phase 1 : relax on the new seed point as Vanek                           */
/*  phase 2 : assign the rest of the nodes to one of the existing            */
/*            aggregate (attach_scheme), if possible.                        */
/*  phase 3 : see if the un-aggregated nodes have enough neighbors           */
/*            (min_nodes_per_aggregate) to form its own aggregate            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_CoarsenUncoupledCore(ML_Aggregate *ml_ag, ML_Comm *comm,
                      ML_Operator *Amat, int *mat_indx, int *bdry_array,
                      int *aggr_count_in, int **aggr_index_in, char *true_bdry)
{
   int     i, j, k, m, kk, inode = 0, jnode, nbytes, length, Nrows;
   int     select_flag, aggr_count, index, mypid, inode2;
   int     *aggr_index, *itmp_array = NULL, count;
   int     *aggr_stat, ordering;
   double  printflag;
   int     *randomVector, aggr_cnt_leng, *aggr_cnt_array;
   int     min_nodes_per_aggregate, max_neigh_selected;
#define newstuff
#ifndef newstuff
   int     attach_scheme;
#endif
   ML_Node       *node_head=NULL, *node_tail=NULL, *new_node=NULL;
   ML_SuperNode  *aggr_head=NULL, *aggr_curr=NULL, *supernode=NULL;
#ifndef newstuff
   int *int_buf = NULL, maxcount, mincount, search_flag;
#endif

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid                   = comm->ML_mypid;
   min_nodes_per_aggregate = ml_ag->min_nodes_per_aggregate;
   max_neigh_selected      = ml_ag->max_neigh_already_selected;
   ordering                = ml_ag->ordering;
#ifndef newstuff
   attach_scheme           = ml_ag->attach_scheme;
#endif
   printflag               = ml_ag->print_flag;
   Nrows                   = mat_indx[0] - 1;

   /* ============================================================= */
   /* aggr_stat indicates whether this node has been aggreated, and */
   /* aggr_index stores the aggregate number where this node has    */
   /* been aggregated into.                                         */
   /* ============================================================= */

   nbytes = Nrows * sizeof( int );
   if ( nbytes > 0 )
   {
      ML_memory_alloc((void**) &aggr_index, (unsigned int) nbytes, "AMA");
      ML_memory_alloc((void**) &aggr_stat, (unsigned int) nbytes, "AMB");
   } else aggr_index = aggr_stat = NULL;

   for ( i = 0; i < Nrows; i++ ) aggr_stat[i] = ML_AGGR_READY;
   for ( i = 0; i < Nrows; i++ ) aggr_index[i] = -1;

   /* ============================================================= */
   /* count number of boundary points                               */
   /* ============================================================= */

   m = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      if ( bdry_array[i] == 1 )
      {
         aggr_stat[i] = ML_AGGR_BDRY;
         m++;
      }
   }
   k = ML_Comm_GsumInt( comm, m);
   if ( mypid == 0 && printflag  < ML_Get_PrintLevel())
   {
      printf("Aggregation(UC) : Phase 0 - no. of bdry pts  = %d \n",k);
   }

   /* ============================================================= */
   /* Set up the data structures for aggregation                    */
   /* ============================================================= */

   aggr_count = 0;
   aggr_head = NULL;
   aggr_cnt_leng = Nrows / 5 + 2;
   nbytes = aggr_cnt_leng * sizeof( int );
   if ( nbytes > 0 )
   {
      ML_memory_alloc((void**) &aggr_cnt_array, (unsigned int) nbytes, "AME");
      for ( i = 0; i < aggr_cnt_leng; i++ ) aggr_cnt_array[i] = 0;
   } else
      aggr_cnt_array = NULL;

   /* ============================================================= */
   /* Phase 1  :                                                    */
   /*    for all nodes, form a new aggregate with its neighbors     */
   /*    if the number of its neighbors having been aggregated does */
   /*    not exceed a given threshold                               */
   /*    (max_neigh_selected = 0 ===> Vanek's scheme)               */
   /* ============================================================= */

   if ( ordering == 1 )       /* random ordering */
   {
      nbytes = Nrows * sizeof(int);
      ML_memory_alloc((void**) &randomVector, (unsigned int) nbytes, "AMF");
      for (i = 0; i < Nrows; i++) randomVector[i] = i;
      ML_randomize(Nrows, randomVector);
   }
   else if ( ordering == 2 )  /* graph ordering */
   {
      new_node = (ML_Node *) ML_allocate(sizeof(ML_Node));
      new_node->node_id = 0;
      node_head = new_node;
      node_tail = new_node;
      new_node->next = NULL;
   }

   inode2 = 0;
   while ( inode2 < Nrows)
   {
      /*------------------------------------------------------ */
      /* pick the next node to aggregate                       */
      /*------------------------------------------------------ */

      if      ( ordering == 0 ) inode = inode2++;
      else if ( ordering == 1 ) inode = randomVector[inode2++];
      else if ( ordering == 2 )
      {
         if ( node_head == NULL )
         {
            for ( jnode = 0; jnode < Nrows; jnode++ )
            {
               if ( aggr_stat[jnode] == ML_AGGR_READY )
               {
                  new_node = (ML_Node *) ML_allocate(sizeof(ML_Node));
                  new_node->node_id = jnode;
                  node_head = new_node;
                  node_tail = new_node;
                  new_node->next = NULL;
                  break;
               }
            }
         }
         if ( node_head == NULL ) break;
         new_node = node_head;
         inode = new_node->node_id;
         node_head = new_node->next;
         ML_free(new_node);
      }

      /*------------------------------------------------------ */
      /* consider further only if the node is in READY mode    */
      /*------------------------------------------------------ */

      if ( aggr_stat[inode] == ML_AGGR_READY )
      {
         length = mat_indx[inode+1] - mat_indx[inode] + 1;
         supernode = (ML_SuperNode *) ML_allocate(sizeof(ML_SuperNode));
         supernode->list = (int*) ML_allocate(length*sizeof(int));

         if ((supernode->list) == NULL)
         {
            printf("Error:couldn't allocate memory for supernode! %d\n",
                            length);
            exit(1);
         }

         supernode->maxlength = length;
         supernode->length = 1;
         supernode->list[0] = inode;
         select_flag = 1;

         /*--------------------------------------------------- */
         /* count the no. of neighbors having been aggregated  */
         /*--------------------------------------------------- */

         count = 0;
         for (jnode=mat_indx[inode];jnode<mat_indx[inode+1];jnode++)
         {
            index = mat_indx[jnode];
            if ( index < Nrows )
            {
               if ( aggr_stat[index] == ML_AGGR_READY ||
                    aggr_stat[index] == ML_AGGR_NOTSEL )
                  supernode->list[supernode->length++] = index;
               else if ( aggr_stat[index] != ML_AGGR_BDRY ) count++;
            }
         }

         /*--------------------------------------------------- */
         /* if there are too many neighbors aggregated or the  */
         /* number of nodes in the new aggregate is too few,   */
         /* don't do this one                                  */
         /*--------------------------------------------------- */

         if ( count > max_neigh_selected ) select_flag = 0;

         if (select_flag != 1 ||
             supernode->length < min_nodes_per_aggregate)
         {
            aggr_stat[inode] = ML_AGGR_NOTSEL;
            ML_free( supernode->list );
            ML_free( supernode );
            if ( ordering == 2 ) /* if graph ordering */
            {
               for (jnode=mat_indx[inode];jnode<mat_indx[inode+1];jnode++)
               {
                  index = mat_indx[jnode];
                  if ( aggr_stat[index] == ML_AGGR_READY )
                  {
                     new_node = (ML_Node *) ML_allocate(sizeof(ML_Node));
                     new_node->node_id = index;
                     new_node->next = NULL;
                     if ( node_head == NULL )
                     {
                        node_head = new_node;
                        node_tail = new_node;
                     } else {
                        node_tail->next = new_node;
                        node_tail = new_node;
                     }
                  }
               }
            }
         }
         else
         {
            for ( j = 0; j < supernode->length; j++ )
            {
               jnode = supernode->list[j];
               aggr_stat[jnode] = ML_AGGR_SELECTED;
               aggr_index[jnode] = aggr_count;
               if ( ordering == 2 ) /* if graph ordering */
               {
                  for (kk=mat_indx[jnode];kk<mat_indx[jnode+1];kk++)
                  {
                     if ( aggr_stat[mat_indx[kk]] == ML_AGGR_READY )
                     {
                        new_node = (ML_Node *) ML_allocate(sizeof(ML_Node));
                        new_node->node_id = mat_indx[kk];
                        new_node->next = NULL;
                        if ( node_head == NULL )
                        {
                           node_head = new_node;
                           node_tail = new_node;
                        } else {
                           node_tail->next = new_node;
                           node_tail = new_node;
                        }
                     }
                  }
               }
            }
            supernode->next = NULL;
            supernode->index = aggr_count;
            if ( aggr_count == 0 )
            {
               aggr_head = supernode;
               aggr_curr = supernode;
            }
            else
            {
               aggr_curr->next = supernode;
               aggr_curr = supernode;
            }
            aggr_cnt_array[aggr_count++] = supernode->length;
            if ( aggr_count >= aggr_cnt_leng )
            {
               itmp_array = aggr_cnt_array;
               aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
               nbytes = aggr_cnt_leng * sizeof( int );
               ML_memory_alloc((void**) &aggr_cnt_array, (unsigned int) nbytes,"AMG");
               for ( k = 0; k < aggr_count; k++ )
                  aggr_cnt_array[k] = itmp_array[k];
               ML_memory_free((void**) &itmp_array);
            }
         }
      }
   }
   if ( ordering == 1 ) ML_memory_free((void**) &randomVector);
   else if ( ordering == 2 )
   {
      while ( node_head != NULL )
      {
         new_node = node_head;
         node_head = new_node->next;
         ML_free( new_node );
      }
   }

   m = 0;
   for ( i = 0; i < Nrows; i++ )
      if ( aggr_stat[i] == ML_AGGR_READY ) m++;
   k = ML_Comm_GsumInt( comm, m);
   if ( k > 0 && mypid == 0 && printflag  < ML_Get_PrintLevel())
      printf("Aggregation(UC) : Phase 1 (WARNING) - %d READY nodes left\n",k);
   m = 0;
   for ( i = 0; i < Nrows; i++ )
      if ( aggr_stat[i] == ML_AGGR_SELECTED ) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows);
   j = ML_Comm_GsumInt( comm, aggr_count );

   if ( mypid == 0 && printflag  < ML_Get_PrintLevel())
   {
      printf("Aggregation(UC) : Phase 1 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(UC) : Phase 1 - total aggregates = %d \n",j);
   }
#ifdef newstuff

   ML_Aggregate_Phase2_3_Cleanup(ml_ag, Amat, &aggr_count, Nrows, aggr_index,
   				   Nrows, comm, true_bdry, "UC_Phase2_3",NULL);
#else
   /* ============================================================= */
   /* Phase 2 : aggregate the rest of the nodes into one of the     */
   /*           existing LOCAL aggregates. (attach_scheme)          */
   /* ============================================================= */

   count = 0;
   for ( inode = 0; inode < Nrows; inode++ )
   {
      /* ---------------------------------------------------------- */
      /* for all nodes that have not been aggregated                */
      /* ---------------------------------------------------------- */

      if ( aggr_stat[inode] == ML_AGGR_NOTSEL ||
           aggr_stat[inode] == ML_AGGR_READY )
      {
         if ( attach_scheme == ML_AGGR_MINRANK )
         {
            /* ---------------------------------------------------- */
            /* search for a neighboring aggregate that has the      */
            /* fewest number of nodes                               */
            /* ---------------------------------------------------- */

            search_flag = 0;
            mincount = 100000;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1];
                 jnode++)
            {
               index = mat_indx[jnode];
               if ( index < Nrows )
               {
                  if ( aggr_stat[index] == ML_AGGR_SELECTED )
                  {
                     search_flag = 1;
                     m = aggr_index[index];
                     if ( aggr_cnt_array[m] < mincount )
                     {
                        mincount = aggr_cnt_array[m];
                        k = index;
                     }
                  }
               }
            }
            if ( search_flag == 1 )
            {
               index = k;
               m = aggr_index[index];
            }

         }
         else if ( attach_scheme == ML_AGGR_MAXLINK )
         {
            /* ---------------------------------------------------- */
            /* search for a neighboring aggregate that has the most */
            /* connection to my node                                */
            /* ---------------------------------------------------- */

            search_flag = 0;
            length = mat_indx[inode+1] - mat_indx[inode];
            nbytes = length * sizeof( int );
            if ( nbytes > 0 )
               ML_memory_alloc((void**) &int_buf, (unsigned int) nbytes, "AGR");
            length = 0;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1];
                 jnode++)
            {
               index = mat_indx[jnode];
               if ( aggr_index[index] >= 0 )
                  int_buf[length++] = aggr_index[index];
            }
            ML_sort(length, int_buf);
            m = -1;
            maxcount = 0;
            if ( length > 0 ) {k = int_buf[0]; j = 1; m = k;}
            for ( jnode = 1; jnode < length; jnode++ )
            {
               if ( int_buf[jnode] == k ) j++;
               else
               {
                  if ( j > maxcount )
                  {
                     maxcount = j;
                     m = k;
                  }
                  k = int_buf[jnode];
                  j = 1;
               }
            }
            if ( m >= 0 ) search_flag = 1;
            if ( nbytes > 0 ) ML_memory_free((void**) &int_buf);
         } else {
            printf("ML_Aggregate_CoarsenUncoupled error : invalid scheme.\n");
            exit(1);
         }

         /* ------------------------------------------------------- */
         /* if found, add the node to the existing aggregate        */
         /* ------------------------------------------------------- */

         if ( search_flag == 1 )
         {
            aggr_cnt_array[m]++;
            aggr_index[inode] = m;
            aggr_stat[inode] = ML_AGGR_SELECTED2;
            count++;
         }
      }
   }
   for ( i = 0; i < Nrows; i++ )
      if (aggr_stat[i] == ML_AGGR_SELECTED2) aggr_stat[i] = ML_AGGR_SELECTED;

   m = 0;
   for ( i = 0; i < Nrows; i++ )
      if ( aggr_stat[i] == ML_AGGR_SELECTED ) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows);
   j = ML_Comm_GsumInt( comm, aggr_count );

   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("Aggregation(UC) : Phase 2 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(UC) : Phase 2 - total aggregates = %d \n",j);
   }

   /* ============================================================= */
   /* Phase 3 : for the un-aggregated nodes, form a new aggregate   */
   /* ============================================================= */

   for ( inode = 0; inode < Nrows; inode++ )
   {
      if (aggr_stat[inode] == ML_AGGR_READY ||
          aggr_stat[inode] == ML_AGGR_NOTSEL )
      {
         count = 1;
         for (jnode = mat_indx[inode]; jnode < mat_indx[inode+1]; jnode++)
         {
            index = mat_indx[jnode];
            if ( index < Nrows && aggr_stat[index] != ML_AGGR_SELECTED )
               count++;
         }
         length = mat_indx[inode+1] - mat_indx[inode];

         /* ------------------------------------------------------- */
         /* if enough neighbors have not been aggregated, form one  */
         /* ------------------------------------------------------- */

         supernode = (ML_SuperNode *) ML_allocate(sizeof(ML_SuperNode));
         supernode->list = (int*) ML_allocate(count*sizeof(int));
         if ((supernode->list) == NULL)
         {
            printf("ML_Aggregate_Coarsen - couldn't allocate memory.\n");
            exit(1);
         }

         supernode->maxlength = count;
         supernode->length = 1;
         supernode->list[0] = inode;

         for (jnode = mat_indx[inode]; jnode < mat_indx[inode+1]; jnode++)
         {
            index = mat_indx[jnode];
            if ( index < Nrows&& aggr_stat[index] != ML_AGGR_SELECTED &&
                 aggr_stat[index] != ML_AGGR_BDRY )
               supernode->list[supernode->length++] = index;
         }
         if ( supernode->length > 1 )
         {
            for ( j = 0; j < supernode->length; j++ )
            {
               jnode = supernode->list[j];
               aggr_stat[jnode] = ML_AGGR_SELECTED;
               aggr_index[jnode] = aggr_count;
            }
            supernode->next = NULL;
            supernode->index = aggr_count;
            if ( aggr_count == 0 )
            {
               aggr_head = supernode;
               aggr_curr = supernode;
            }
            else
            {
               aggr_curr->next = supernode;
               aggr_curr = supernode;
            }
            aggr_cnt_array[aggr_count++] = supernode->length;
            if ( aggr_count >= aggr_cnt_leng )
            {
               itmp_array = aggr_cnt_array;
               aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
               nbytes = aggr_cnt_leng * sizeof( int );
               ML_memory_alloc((void**) &aggr_cnt_array, (unsigned int) nbytes, "AGL");
               for ( k = 0; k < aggr_count; k++ )
                  aggr_cnt_array[k] = itmp_array[k];
               ML_memory_free((void**) &itmp_array);
            }
         }
         else
         {
            for ( j = 0; j < supernode->length; j++ )
            {
               jnode = supernode->list[j];
               aggr_stat[jnode] = ML_AGGR_BDRY;
            }
            if ( supernode->maxlength > 0 ) ML_free( supernode->list );
            ML_free( supernode );
         }
      }
   }

   m = 0;
   for ( i = 0; i < Nrows; i++ )
      if ( aggr_stat[i] == ML_AGGR_SELECTED ) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows);
   j = ML_Comm_GsumInt( comm, aggr_count );

   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("Aggregation(UC) : Phase 3 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(UC) : Phase 3 - total aggregates = %d \n",j);
   }

   /* ============================================================= */
   /* check for error                                               */
   /* ============================================================= */

   m = 0;
   for ( i = 0; i < Nrows; i++ )
      if (aggr_stat[i] != ML_AGGR_SELECTED && aggr_stat[i] != ML_AGGR_BDRY) m++;
   k = ML_Comm_GsumInt( comm, m);
   if ( k > 0 && mypid == 0 )
   {
      printf("Aggregation (UC) error : not all nodes processed.\n");
      exit(1);
   }
#endif

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &aggr_stat);
   ML_memory_free((void**) &aggr_cnt_array);
   aggr_curr = aggr_head;
   while ( aggr_curr != NULL )
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->maxlength > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }
   (*aggr_count_in) = aggr_count;
   (*aggr_index_in) = aggr_index;

   return aggr_count;
}

