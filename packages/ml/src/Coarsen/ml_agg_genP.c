/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators                                */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL) and Ray Tuminaro (SNL)                */
/* Date          : August, 1999                                              */
/* ************************************************************************* */
/* ************************************************************************* */

#include <math.h>
#include "ml_struct.h"
#include "ml_op_utils.h"
#include "ml_agg_genP.h"

extern int ML_AGG_Amat_Getrows(void *data, int N_requested_rows, 
               int requested_rows[], int allocated_space, int columns[], 
               double values[], int row_lengths[]);

/* ************************************************************************* */
/* wrapper function as smoother                                              */
/* ------------------------------------------------------------------------- */

int ML_AGG_Smoother_Wrapper(void *obj, int leng1, double *outvec, int leng2,
                            double *invec)
{
   ML *ml;
   ml = (ML *) obj;
   ML_Iterate( ml, outvec, invec );
   return 1;
}

/* ************************************************************************* */
/* generate multilevel hierarchy based on Vanek's method                     */
/* ------------------------------------------------------------------------- */

int ML_Gen_MGHierarchy_UsingAggregation(ML *ml, int start, 
                       int increment_or_decrement, ML_Aggregate *ag)
{
   int    level, idata;
   double dnnz = 0;
   ML_Aggregate *ml_ag;
#ifdef ML_TIMING
   double t0;
#endif

   /* ----------------------------------------------------------------- */
   /* if user does not provide a ML_Aggregate object, create a default  */
   /* ----------------------------------------------------------------- */

   if ( ag == NULL ) ML_Aggregate_Create( &ml_ag );
   else ml_ag=ag;
   ML_Aggregate_Set_MaxLevels( ml_ag, ml->ML_num_levels);
   ML_Aggregate_Set_StartLevel( ml_ag, start );

   /* ----------------------------------------------------------------- */
   /* create multilevel hierarchy                                       */
   /* ----------------------------------------------------------------- */

   idata = 0;
   idata = ML_gmax_int(idata, ml->comm);
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag ) 
      ML_Aggregate_Print( ml_ag );
#ifdef ML_TIMING
   t0 = GetClock();
#endif
   idata = ML_gmax_int(idata, ml->comm);

   if (increment_or_decrement == ML_INCREASING)
   {
      /* -------------------------------------------------------------- */
      /* coarse scheme == 4 ==> domain decomposition                    */
      /* -------------------------------------------------------------- */
      if ( ml_ag->coarsen_scheme == 4 )
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Increment_Two_Level,
                     ML_AGG_Gen_DDProlongator, NULL, ML_INTERNAL, ml_ag);
      }
      else
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Increment_Level,
                     ML_AGG_Gen_Prolongator, NULL, ML_INTERNAL, ml_ag);
      }
   }
   else if (increment_or_decrement == ML_DECREASING)
   {
      if ( ml_ag->coarsen_scheme == 4 )
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Decrement_Two_Level,
                     ML_AGG_Gen_DDProlongator, NULL, ML_INTERNAL, ml_ag);
      }
      else
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Decrement_Level,
                     ML_AGG_Gen_Prolongator, NULL, ML_INTERNAL, ml_ag);
      }
   }
   else 
   {
      if ( ml->comm->ML_mypid == 0 ) 
      {
         printf("ML_Gen_MGHierarchy_UsingAggregation : Unknown ");
         printf("increment_or_decrement choice\n");
      }
      exit(1);
   }
#ifdef ML_TIMING
   t0 = GetClock() - t0;
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag ) 
      printf("Aggregation total setup time = %e seconds\n", t0);
#endif

   /* ----------------------------------------------------------------- */
   /* compute operator complexity                                       */
   /* ----------------------------------------------------------------- */

   if (increment_or_decrement == ML_INCREASING)
      dnnz = ml->Amat[level-start-1].N_nonzeros;
   else if (increment_or_decrement == ML_DECREASING)
      dnnz = ml->Amat[start+1-level].N_nonzeros;
   dnnz = ML_gsum_double( dnnz, ml->comm );
   ml_ag->operator_complexity += dnnz;

   idata = ML_gmax_int(idata, ml->comm);
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag ) 
      ML_Aggregate_Print_Complexity( ml_ag );
   idata = ML_gmax_int(idata, ml->comm);

   if ( ag == NULL ) ML_Aggregate_Destroy( &ml_ag );
   return(level);
}

/* ************************************************************************* */
/* generate multilevel hierarchy given a subroutine for generating           */
/* prolongation operators (version 2 : with ML_Aggregate given)              */
/* ------------------------------------------------------------------------- */

int ML_Gen_MGHierarchy(ML *ml, int fine_level,
        int (*next_level)(ML *, int, ML_Operator *, ML_Aggregate *ag2),
        int (*user_gen_prolongator)(ML *, int, int, void *, ML_Aggregate *),
        void *data, int internal_or_external, ML_Aggregate *ag)
{
   int level, next, flag, count=1;
#ifdef ML_TIMING
   double t0;
#endif

   level = fine_level;
   next  = next_level(ml, level, &(ml->Amat[fine_level]), ag);

   while (next >= 0) 
   {
      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("ML_Gen_MGHierarchy : applying coarsening \n");

      if (internal_or_external == ML_INTERNAL)
      {
         flag = user_gen_prolongator(ml, level, next,
                                     (void*)&(ml->Amat[level]),ag);
      }
      else 
      {
         flag = user_gen_prolongator(ml, level, next, data, ag);
      }
      if (flag < 0) break;
      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("ML_Gen_MGHierarchy : applying coarsening \n");
      ML_Gen_Restrictor_TransP(ml, level, next);

      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("ML_Gen_MGHierarchy : Gen_RAP\n");

#ifdef ML_TIMING
      t0 = GetClock();
#endif
      ML_Gen_AmatrixRAP(ml, level, next);
#ifdef ML_TIMING
      t0 = GetClock() - t0;
      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("RAP time for level %2d = %e\n", level, t0);
#endif

      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("ML_Gen_MGHierarchy : Gen_RAP done\n");

      level = next;
      next  = next_level(ml, next, &(ml->Amat[next]), ag);
      count++;
   }
   return(count);
}

/* ************************************************************************* */
/* generate smooth prolongator                                               */
/* ------------------------------------------------------------------------- */

int ML_AGG_Gen_Prolongator(ML *ml,int level, int clevel, void *data,
                             ML_Aggregate *ag)
{
   int         Ncoarse, Nfine, gNfine, gNcoarse;
   double      max_eigen = -1.;
   ML_Operator *Amat, *Pmatrix, *AGGsmoother = NULL;
   struct      ML_AGG_Matrix_Context widget;
   ML_Krylov   *kdata;

#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   widget.near_bdry = NULL;
   Amat     = (ML_Operator *) data;
   Amat->num_PDEs = ag->num_PDE_eqns;
   /*
   widget.near_bdry = (char *) ML_allocate(sizeof(char)*Amat->outvec_leng);
   ML_AGG_Compute_Near_Bdry(Amat, widget.near_bdry);
   */

   Nfine    = Amat->outvec_leng;
   gNfine   = ML_Comm_GsumInt( ml->comm, Nfine);
   ML_Aggregate_Set_CurrentLevel( ag, level );
   Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&Pmatrix,ml->comm);
   gNcoarse = ML_Comm_GsumInt( ml->comm, Ncoarse);
   if ( gNcoarse == 0 || ((1.0*gNfine)/(1.0*gNcoarse+0.1) < 1.05) )
   {
      if ( Pmatrix != NULL ) ML_Operator_Destroy(Pmatrix);
      return -1;
   }

#ifdef ML_NONSYMMETRIC
   max_eigen = ML_Operator_MaxNorm(Amat, ML_TRUE);
   widget.omega  = ag->smoothP_damping_factor / max_eigen;
   ml->spectral_radius[level] = max_eigen;
#else

   if ( ag->smoothP_damping_factor != 0.0 )
   {
      if ( ag->spectral_radius_scheme == 1 ) /* compute it using CG */
      {
         kdata = ML_Krylov_Create( ml->comm );
         ML_Krylov_Set_PrintFreq( kdata, 0 );
         ML_Krylov_Set_ComputeEigenvalues( kdata );
         ML_Krylov_Set_Amatrix(kdata, Amat);
         ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
         max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);
         ML_Krylov_Destroy( &kdata );
         if ( max_eigen <= 0.0 )
         {
            printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
            max_eigen = 1.0;
         }
         if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
            printf("Gen_Prolongator : max eigen = %e \n", max_eigen);

         widget.omega  = ag->smoothP_damping_factor / max_eigen;
         ml->spectral_radius[level] = max_eigen;
      }
      else   /* using matrix max norm */
      {
         max_eigen = ML_Operator_MaxNorm(Amat, ML_TRUE);
         widget.omega  = ag->smoothP_damping_factor / max_eigen;
         ml->spectral_radius[level] = max_eigen;
      }
   }
   else  /* damping fact = 0 ==> no need to compute spectral radius */
   {
      ml->spectral_radius[level] = 1.0;
      widget.omega  = 0.0;
   }
#endif

   widget.drop_tol = ag->drop_tol_for_smoothing;
   widget.Amat   = &(ml->Amat[level]);
   widget.aggr_info = ag->aggr_info[level];
   AGGsmoother = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(AGGsmoother, widget.Amat->invec_leng,
                        widget.Amat->outvec_leng, ML_EXTERNAL,&widget,
                        widget.Amat->matvec->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(AGGsmoother, ML_EXTERNAL,
                          widget.Amat->getrow->Nrows, 
                          ML_AGG_JacobiSmoother_Getrows);
   ML_CommInfoOP_Clone(&(AGGsmoother->getrow->pre_comm),
                          widget.Amat->getrow->pre_comm);

   ML_2matmult(AGGsmoother, Pmatrix, &(ml->Pmat[clevel]) );

   ML_Operator_Destroy(Pmatrix);
   ML_Operator_Destroy(AGGsmoother);
   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));

   if (widget.near_bdry != NULL) ML_free(widget.near_bdry);
#ifdef ML_TIMING
   ml->Pmat[clevel].build_time =  GetClock() - t0;
   ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
   return 0;
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level        */
/* number larger than the fine levels                                        */
/* ------------------------------------------------------------------------- */

int ML_AGG_Increment_Level(ML *ml, int current_level, ML_Operator *Amat,
                           ML_Aggregate *ag)
{
   int total_size, temp;

   if (current_level == ml->ML_num_levels-1) return(-1);

   total_size = Amat->invec_leng;
   ML_gsum_vec_int(&total_size, &temp, 1, ml->comm);
   if ( total_size <= ag->max_coarse_size ) return(-1);

   return(current_level+1);
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level number */
/* smaller than the fine levels                                              */
/* ------------------------------------------------------------------------- */

int ML_AGG_Decrement_Level(ML *ml, int current_level, ML_Operator *Amat,
                           ML_Aggregate *ag)
{
   int total_size, temp;

   if (current_level == 0 ) return(-1);

   total_size = Amat->invec_leng;
   ML_gsum_vec_int(&total_size, &temp, 1, ml->comm);
   if ( total_size <= ag->max_coarse_size ) return(-1);

   return(current_level-1);
}

/* ************************************************************************* */
/* function for enforcing a 2-level scheme                                   */
/* ------------------------------------------------------------------------- */

int ML_AGG_Increment_Two_Level(ML *ml,int current_level,ML_Operator *Amat,
                               ML_Aggregate *ag)
{
   (void) Amat;
   (void) ml;
   if ( current_level == ag->begin_level ) return (current_level+1);
   return(-1);
}

/* ************************************************************************* */
/* function for enforcing a 2-level scheme                                   */
/* ------------------------------------------------------------------------- */

int ML_AGG_Decrement_Two_Level(ML *ml,int current_level,ML_Operator *Amat,
                               ML_Aggregate *ag)
{
   (void) Amat;
   (void) ml;
   if ( current_level == ag->begin_level ) return (current_level-1);
   return(-1);
}

/* ************************************************************************* */
/* ************************************************************************* */
/* getrow function for the aggregation tentative prolongator                 */
/* ------------------------------------------------------------------------- */

int ML_AGG_JacobiSmoother_Getrows(void *data, int N_requested_rows, 
   int requested_rows[], int allocated_space, int columns[], 
   double values[], int row_lengths[])
{
   struct ML_AGG_Matrix_Context *widget;
   ML_GetrowFunc  *getrow_obj;
   int            info, diag = -1, i, j /*, *aggr_info */;
   double         diag_val = 1.0, dropped, threshold = 0.0;

   widget = (struct ML_AGG_Matrix_Context *) data;
   if (widget->near_bdry != NULL) {
     if (widget->near_bdry[requested_rows[0]] == 'T') {
       if (allocated_space < 1) return(0);
       columns[0] = requested_rows[0];
       values[0]  = 1.0;
       row_lengths[0] = 1;
       return(1);
     }
   }

   /* ----------------------------------------------------------------- */
   /* error checking                                                    */
   /* ----------------------------------------------------------------- */

   getrow_obj = widget->Amat->getrow;
   if (N_requested_rows > 1) 
   {
      printf("Too bad. This routine only works with 1 row at a time\n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* if omega = 0, just return identity                                */
   /* ----------------------------------------------------------------- */

   if ( widget->omega == 0.0 )
   {
      row_lengths[0] = 1;
      values[0] = 1.0;
      columns[0] = requested_rows[0];
      return 1;
   }

   /* ----------------------------------------------------------------- */
   /* fetch row                                                         */
   /* ----------------------------------------------------------------- */

   if ( getrow_obj->ML_id == ML_EXTERNAL) 
      info = getrow_obj->external(widget->Amat->data, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
   else if ( getrow_obj->ML_id == ML_INTERNAL) 
      info = getrow_obj->internal(widget->Amat, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
   else 
   {
      printf("Invalid getrow id (%d)\n",getrow_obj->ML_id);
      exit(1);
   }
   if (info == 0) return(0);

   /* ----------------------------------------------------------------- */
   /* compute threshold for dropping                                    */
   /* ----------------------------------------------------------------- */

   if ( widget->drop_tol > 0.0 )
   {
      for (i = 0; i < row_lengths[0]; i++) 
      {
         if (columns[i] == requested_rows[0]) 
         {
            threshold = fabs(values[i])*widget->drop_tol;
            break;
         }
      }
      j = 0;
      dropped = 0.0;
      for (i = 0; i < row_lengths[0]; i++) 
      {
         if ( fabs(values[i]) >= threshold) 
         {
            columns[j] = columns[i];
            values[j]  = values[i];
            if (columns[j] == requested_rows[0]) { diag = j; }
            j++;
         }
         else dropped += values[i];
      }
      row_lengths[0] = j;
   }
   else
   {
      dropped = 0.0;
      for (i = 0; i < row_lengths[0]; i++) 
         if (columns[i] == requested_rows[0]) { diag = i; break;}
   }

   /* ----------------------------------------------------------------- */
   /* if diagonal is not found, append one                              */
   /* ----------------------------------------------------------------- */

   if (diag == -1) 
   {
      if (row_lengths[0] >= allocated_space) return(0);
      columns[row_lengths[0]] = requested_rows[0];
      values[row_lengths[0]]  = 0.0;
      diag = row_lengths[0];
      row_lengths[0]++;
   }
   else diag_val = values[diag];

   values[diag] += dropped;

   /* ----------------------------------------------------------------- */
   /* The following segment is for filtering (aggregate - not used)     */
   /* ----------------------------------------------------------------- */

/*
   aggr_info = widget->aggr_info;
   N = widget->Amat->outvec_leng;
   for (i = 0; i < row_lengths[0]; i++) 
   {
      if (columns[i] < N &&
          aggr_info[columns[i]] != aggr_info[requested_rows[0]])
      {
         values[diag] += values[i];
         values[i] = 0.0;
      }
   }
   N = 0;
   for (i = 0; i < row_lengths[0]; i++)
   {
      if ( values[i] != 0.0 ) 
      {
         values[N] = values[i]; 
         columns[N++] = columns[i];}
      }
   }
   row_lengths[0] = N;
   diag_val = values[diag];
*/

   /* ----------------------------------------------------------------- */
   /* compute I - omega D^{-1} A                                        */
   /* ----------------------------------------------------------------- */

   for (i = 0; i < row_lengths[0]; i++) 
      values[i] *= (-widget->omega)/diag_val;
   values[diag] += 1.;

   return(1);
}

/* ************************************************************************* */
/* getrow function for the aggregation tentative prolongator                 */
/* ------------------------------------------------------------------------- */

int ML_AGG_Amat_Getrows(void *data, int N_requested_rows, 
   int requested_rows[], int allocated_space, int columns[], 
   double values[], int row_lengths[])
{
   struct ML_AGG_Matrix_Context *widget;
   ML_GetrowFunc  *getrow_obj;
   int            info;

   widget = (struct ML_AGG_Matrix_Context *) data;
   getrow_obj = widget->Amat->getrow;
   if (N_requested_rows > 1) 
   {
      printf("Too bad. This routine only works with 1 row at a time\n");
      exit(1);
   }

   if ( getrow_obj->ML_id == ML_EXTERNAL) 
      info = getrow_obj->external(widget->Amat->data, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
   else if ( getrow_obj->ML_id == ML_INTERNAL) 
      info = getrow_obj->internal(widget->Amat, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
   else 
   {
      printf("Invalid getrow id (%d)\n",getrow_obj->ML_id);
      exit(1);
   }
   if (info == 0) return(0);

   return(1);
}

/* ************************************************************************* */
/* generate smooth prolongator for 2-level DD method                         */
/* ------------------------------------------------------------------------- */

int ML_AGG_Gen_DDProlongator(ML *ml,int level, int clevel, void *data,
                             ML_Aggregate *ag)
{
   int          i, j, Nfine, nbytes, newNlevels, nnz, *col_ind;
   int          k, newClevel, lengc, lengf, ap_ncols, *ap_cols;
   int          *new_ia, *new_ja, p_ncols, *p_cols, max_nz_per_row;
   double       max_eigen, norm, *darray, *darray2, **ap_aa;
   double       *diagonal, *new_val, **p_aa, *col_val;
   ML_Operator  *Amat, *tentP, *APMat;
   ML_Krylov    *kdata;
   struct ML_AGG_Matrix_Context widget;
   struct ML_AGG_Matrix_Context *context;
   ML           *newml;
   ML_Aggregate *newag;
   struct  ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm *aggr_comm;
   ML_GetrowFunc *getrow_obj;
   int           (*getrowfunc)(void *,int,int*,int,int*,double*,int*);

#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   /* ----------------------------------------------------------------- */
   /* coarsen local smoothed aggregation method                         */
   /* ----------------------------------------------------------------- */

   if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
      printf("Aggregation : building multilevel hierarchy at level %d\n",level);
   widget.near_bdry = NULL; 
   Amat     = (ML_Operator *) data;
   Nfine    = Amat->outvec_leng;
   getrow_obj = Amat->getrow;
   if (getrow_obj->ML_id == ML_EXTERNAL) getrowfunc = getrow_obj->external;
   else                                  getrowfunc = getrow_obj->internal;
   max_nz_per_row = 100;
   col_ind = (int *)    malloc( max_nz_per_row * sizeof(int) );
   col_val = (double *) malloc( max_nz_per_row * sizeof(double) );
   nnz = 0;
   for ( i = 0; i < Nfine; i++ )
   {
      while (getrowfunc(Amat->data,1,&i,max_nz_per_row,col_ind,col_val,&k)== 0)
      {
         free( col_ind );
         free( col_val );
         max_nz_per_row = max_nz_per_row * 2 + 1;
         col_ind = (int *)    malloc( max_nz_per_row * sizeof(int) );
         col_val = (double *) malloc( max_nz_per_row * sizeof(double) );
      }
      nnz += k;
   }
   free( col_ind );
   free( col_val );
   nnz = ML_Comm_GsumInt( ml->comm, nnz);
   if ( ag->operator_complexity == 0.0 )
   {
      ag->fine_complexity = 1.0 * nnz;
      ag->operator_complexity = 1.0 * nnz;
   }
   else ag->operator_complexity += 1.0 * nnz;

   /* ----------------------------------------------------------------- */
   /* setup local smoothed aggregation method                           */
   /* ----------------------------------------------------------------- */

   if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
      printf("Aggregation : setting up diagonal block at level %d\n",level);

   newNlevels = 15;
   ML_Create(&newml, newNlevels);
   ML_Set_OutputLevel(newml, 0);
   ML_Set_ResidualOutputFrequency(newml, 0);
   ML_Set_Comm_MyRank(newml, 0);
   ML_Set_Comm_Nprocs(newml, 1);
   nbytes = sizeof(struct ML_AGG_Matrix_Context);
   context = (struct ML_AGG_Matrix_Context *) malloc( nbytes );
   context->Amat = Amat;
   context->near_bdry = NULL;
   ML_Init_Amatrix(newml, newNlevels-1, Nfine,  Nfine, (void *) context);
   ML_Set_Amatrix_Matvec(newml,  newNlevels-1, ML_AGG_DD_Matvec);
   newml->Amat[newNlevels-1].data_destroy = ML_AGG_Matrix_Context_Clean;
   newml->Amat[newNlevels-1].N_nonzeros = 5 * Nfine;
   ML_Set_Amatrix_Getrow(newml, newNlevels-1, ML_AGG_DD_Getrow, NULL, Nfine);
   diagonal = (double *) malloc(Nfine * sizeof(double));
   ML_AGG_Extract_Diag(Amat, diagonal);
   ML_Set_Amatrix_Diag( newml, newNlevels-1, Nfine, diagonal);
   free( diagonal );
   ML_Aggregate_Create( &newag );
   ML_Aggregate_Set_OutputLevel( newag, 0 );
   ML_Aggregate_Set_CoarsenScheme_Uncoupled( newag );
   ML_Aggregate_Set_Threshold( newag, 0.08 );
ML_Aggregate_Set_DampingFactor( newag, 0.0/3.0 );
   ML_Aggregate_Set_MaxCoarseSize( newag, 1 );
   ML_Aggregate_Set_PSmootherType( newag, 0 );
   newClevel = ML_Gen_MGHierarchy_UsingAggregation(newml, newNlevels-1,
                                  ML_DECREASING, newag);
   newClevel = newNlevels - newClevel;
   for (k = newNlevels-1; k > newClevel; k--) 
   {
      ML_Gen_Smoother_SymGaussSeidel(newml, k, ML_PRESMOOTHER, 1, 1.);
      ML_Gen_Smoother_SymGaussSeidel(newml, k, ML_POSTSMOOTHER, 1, 1.);
   }
   ML_Gen_CoarseSolverSuperLU( newml, newClevel );
   ML_Gen_Solver(newml, ML_MGV, newNlevels-1, newClevel);
   ML_Aggregate_Destroy( &newag );

   /* ----------------------------------------------------------------- */
   /* set up Krylov solver to compute eigenvalues                       */
   /* ----------------------------------------------------------------- */

   if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
      printf("Aggregation : computing max eigenvalues at level %d\n",level);

/*
   if ( ag->smoothP_damping_factor != 0.0 )
*/
   if ( 1 )
   {
      kdata = ML_Krylov_Create( ml->comm );
      ML_Krylov_Set_PrintFreq( kdata, 0 );
      ML_Krylov_Set_ComputeEigenvalues( kdata );
      ML_Krylov_Set_Amatrix(kdata, Amat);
      ML_Krylov_Set_Precon(kdata, (void *) newml);
      ML_Krylov_Set_PreconFunc(kdata, ML_AGG_DD_Solve);
      ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
      max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);
      ML_Krylov_Destroy( &kdata );
      if ( max_eigen <= 0.0 )
      {
         printf("Gen_DDProlongator warning : max eigen <= 0.0 \n");
         max_eigen = 1.0;
      }
      if ( ml->comm->ML_mypid == 0 ) 
         printf("Gen_DDProlongator : max eigen = %e \n", max_eigen);

      widget.omega  = ag->smoothP_damping_factor / max_eigen;
      ml->spectral_radius[level] = max_eigen;
   }
   else
   {
      widget.omega = 0.0;
      ml->spectral_radius[level] = 1.0;
   }

   /* ================================================================= */
   /* set up smoothed prolongator (I - alpha D^(-1) A) P                */
   /* 1. compute P (local ml P applied to 1)                            */
   /* 2. A * P                                                          */
   /* 3. D^{-1} * A * P                                                 */
   /* 4. P - alpha D^{-1} * A * P                                       */
   /* ================================================================= */

   i = 1;
   j = ML_gmax_int(i, ml->comm );
   if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
      printf("Aggregation : computing tentative prolongators at level %d\n",level);

   /* ----------------------------------------------------------------- */
   /* 1. compute tentP (local ml P applied to 1)                        */
   /* ----------------------------------------------------------------- */

   darray  = (double *) malloc( Nfine * sizeof(double) );
   darray2 = (double *) malloc( Nfine * sizeof(double) );

   for ( i = 0; i < newml->Amat[newClevel].outvec_leng; i++ )  
      darray[i] = 1.0;

   for ( i = newClevel; i < newNlevels-1; i++ )  
   {
      lengc = newml->Amat[i].outvec_leng;
      lengf = newml->Amat[i+1].outvec_leng;
      ML_Operator_ApplyAndResetBdryPts(&(newml->Pmat[i]),lengc,darray,lengf,
                                       darray2);
      for ( j = 0; j < lengf; j++ ) darray[j] = darray2[j];
   }  
   free( darray2 );
   norm = 0.0;
   for ( j = 0; j < Nfine; j++ ) norm += (darray[j] * darray[j]);
   norm = sqrt(norm);
   for (j = 0; j < Nfine; j++) darray[j] = darray[j] / norm;

   nbytes = ( Nfine + 1 ) * sizeof(int);
   ML_memory_alloc((void**)&(new_ia), nbytes, "AD1");
   nbytes = Nfine * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), nbytes, "AD2");
   nbytes = Nfine * sizeof(double);
   ML_memory_alloc((void**)&(new_val), nbytes, "AD3");
   for (i = 0; i <= Nfine; i++) new_ia[i] = i;
   for (i = 0; i < Nfine; i++) new_ja[i] = 0;
/* */
if ( ml->comm->ML_mypid == 0 )
   printf("Tentative prolongator set to 1.\n");
for (i = 0; i < Nfine; i++) darray[i] = 1.0/sqrt((double) Nfine);
/* */ 
   for (i = 0; i < Nfine; i++) new_val[i] = darray[i];

   p_ncols = 1;
   p_cols = (int *) malloc(sizeof(int));
   p_cols[0] = 0;
   p_aa = (double **) malloc(sizeof(double*));
   p_aa[0] = darray;
 
   ML_memory_alloc((void**) &csr_data,sizeof(struct ML_CSR_MSRdata),"AVP");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   tentP = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(tentP,1,Nfine,ML_EMPTY,csr_data,Nfine,NULL,0);
   tentP->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm), "AD4");
   aggr_comm->comm = ml->comm;
   aggr_comm->N_send_neighbors = 0;
   aggr_comm->N_recv_neighbors = 0;
   aggr_comm->send_neighbors = NULL;
   aggr_comm->recv_neighbors = NULL;
   aggr_comm->send_leng = NULL;
   aggr_comm->recv_leng = NULL;
   aggr_comm->send_list = NULL;
   aggr_comm->local_nrows = 1;
   ML_CommInfoOP_Generate( &(tentP->getrow->pre_comm),
                           ML_Aggregate_ExchangeBdry, aggr_comm, ml->comm, 1, 0);
   ML_Operator_Set_Getrow(tentP, ML_EXTERNAL, Nfine, CSR_getrows);
   ML_Operator_Set_ApplyFunc(tentP, ML_INTERNAL, CSR_matvec);

   /* ----------------------------------------------------------------- */
   /* 2. compute AP = A * tentP                                         */
   /* 3. compute P = tentP - alpha * D^{-1} A tentP                     */
   /* ----------------------------------------------------------------- */

/*
   if ( ag->smoothP_damping_factor != 0.0 )
*/
   if ( 1 )
   {
      i = 1;
      j = ML_gmax_int(i, ml->comm );
      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("Aggregation : computing comm pattern of A*tentP at level %d\n",
              level);

      APMat = &(ml->Pmat[clevel]);
      ML_2matmult(Amat, tentP, APMat);
      ML_AGG_Extract_Matrix(APMat, &ap_ncols, &ap_cols, &ap_aa);
 
      i = 1;
      j = ML_gmax_int(i, ml->comm );
      if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
         printf("Aggregation : computing prolongators at level %d\n",level);

      ML_Set_MaxIterations(newml, 10);
      darray  = (double *) malloc( Nfine * sizeof(double) );
      for ( i = 0; i < ap_ncols; i++ )
      {
         for ( j = 0; j < Nfine; j++ ) darray[j] = 0.0;
         /*ML_Solve_MGV( newml, ap_aa[i], darray );*/
         ML_Iterate( newml, darray, ap_aa[i] );
         if ( i == 0 )
         {
            for ( j = 0; j < Nfine; j++ ) 
               ap_aa[i][j] = p_aa[0][j] - widget.omega * darray[j];
         }
         else
         {
            for ( j = 0; j < Nfine; j++ ) 
               ap_aa[i][j] = - widget.omega * darray[j];
         }
      }
      free( darray );
   }
   else
   {
      APMat = &(ml->Pmat[clevel]);
      ap_ncols = p_ncols;
      ap_cols  = p_cols;
      ap_aa = p_aa;
      p_cols = NULL;
      p_aa = NULL;
      p_ncols = 0;
   }
   if ( p_cols != NULL ) free( p_cols );
   for ( i = 0; i < p_ncols; i++ ) free( p_aa[i] );
   if ( p_aa != NULL ) free( p_aa );

   nnz = 0;
   for ( i = 0; i < ap_ncols; i++ )
      for ( j = 0; j < Nfine; j++ ) if ( ap_aa[i][j] != 0.0 ) nnz++;

   nbytes = ( Nfine + 1 ) * sizeof(int);
   ML_memory_alloc((void**)&(new_ia), nbytes, "ADA");
   nbytes = nnz * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), nbytes, "ADB");
   nbytes = nnz * sizeof(double);
   ML_memory_alloc((void**)&(new_val), nbytes, "ADC");

   nnz = 0;
   new_ia[0] = nnz;
   for ( i = 0; i < Nfine; i++ )
   {
      for ( j = 0; j < ap_ncols; j++ ) 
         if ( ap_aa[j][i] != 0.0 ) 
         {
            new_ja[nnz] = ap_cols[j]; 
            new_val[nnz++] = ap_aa[j][i];
         }
      new_ia[i+1] = nnz;
   }
   max_nz_per_row = 0;
   for ( i = 0; i < Nfine; i++ )
   {
      nnz = 0;
      for ( j = 0; j < ap_ncols; j++ )
         if ( ap_aa[j][i] != 0.0 ) nnz++;
      if ( nnz > max_nz_per_row ) max_nz_per_row = nnz;
   }
   
   ML_memory_alloc((void**)&csr_data,sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
   ML_Operator_Set_ApplyFuncData(APMat,1,Nfine,ML_EMPTY,csr_data,
                                 Nfine,NULL,ap_ncols-1);
   APMat->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_Operator_Set_Getrow(APMat, ML_EXTERNAL, Nfine, CSR_getrows);
   ML_Operator_Set_ApplyFunc(APMat, ML_INTERNAL, CSR_matvec);
   APMat->max_nz_per_row = max_nz_per_row;
/*
   if ( ag->smoothP_damping_factor == 0.0 )
   {
      ML_CommInfoOP_Generate( &(APMat->getrow->pre_comm), 
                           ML_Aggregate_ExchangeBdry, aggr_comm, ml->comm, 1, 0);
   }
*/

   free( ap_cols );
   for ( i = 0; i < ap_ncols; i++ ) free( ap_aa[i] );
   free( ap_aa );
   ML_Destroy(&newml);
   ML_Operator_Destroy(tentP);

   i = 1;
   j = ML_gmax_int(i, ml->comm );
   if ( ml->comm->ML_mypid == 0 && ag->print_flag ) 
      printf("Aggregation : building P complete at level %d\n",level);

/*
   ML_Set_Smoother(ml, level, ML_PRESMOOTHER, newml, ML_AGG_Smoother_Wrapper,NULL);
*/
#ifdef ML_TIMING
   ml->Pmat[clevel].build_time =  GetClock() - t0;
   ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
   return 0;
}

/* ************************************************************************* */
/* local matvec                                                              */
/* ------------------------------------------------------------------------- */

int ML_AGG_DD_Matvec(void *obj,int leng1,double p[],int leng2,double ap[])
{
   int         i, j, m, max_row_nnz=100, nRows, index, *col_ind;
   double      dtmp, *col_val;
   ML_Operator *Amat;
   int         (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   struct ML_AGG_Matrix_Context *context;
   ML_GetrowFunc                 *getrow_obj;

   context = (struct ML_AGG_Matrix_Context *) obj;
   Amat    = (ML_Operator *) context->Amat;
   nRows   = Amat->outvec_leng;
   if ( nRows != leng1 || leng1 != leng2 )
   {
      printf("ML_AGG_DD_Matvec ERROR : inleng != outleng.\n");
      exit(-1);
   }
   getrow_obj = Amat->getrow;
   if (getrow_obj->ML_id == ML_EXTERNAL) getrowfunc = getrow_obj->external;
   else                                  getrowfunc = getrow_obj->internal;
   if ( getrowfunc == NULL )
   {
      printf("ML_AGG_DD_Matvec ERROR : null getrowfunc.\n");
      exit(-1);
   }
   col_ind = (int *)    malloc( max_row_nnz * sizeof(int) );
   col_val = (double *) malloc( max_row_nnz * sizeof(double) );

   for ( i = 0; i < nRows; i++ )
   {
      while (getrowfunc(Amat->data,1,&i,max_row_nnz,col_ind,col_val,&m)== 0)
      {
         free( col_ind );
         free( col_val );
         max_row_nnz = max_row_nnz * 2 + 1;
         col_ind = (int *)    malloc( max_row_nnz * sizeof(int) );
         col_val = (double *) malloc( max_row_nnz * sizeof(double) );
      }
      dtmp = 0.0;
      
      for ( j = 0; j < m; j++ )
      {
         index = col_ind[j];
         if ( index < nRows ) dtmp += ( col_val[j] * p[index] );
      }
      ap[i] = dtmp;
   }
   free( col_ind );
   free( col_val );

   return 1;
}

/* ************************************************************************* */
/* local getrow                                                              */
/* ------------------------------------------------------------------------- */

int ML_AGG_DD_Getrow(void *obj,int inNrows, int *rowlist,int alloc_space, 
                     int *col_ind, double *col_val, int *rowcnt)
{
   int         i, count, status, nRows, *local_ind = NULL;
   double      *local_val = NULL;
   ML_Operator *Amat;
   int         (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   struct ML_AGG_Matrix_Context *context;
   ML_GetrowFunc                *getrow_obj;

   if ( inNrows != 1 )
   {
      printf("ML_AGG_DD_Getrow ERROR : inNrows > 1 not supported.\n");
      exit(-1);
   }
   context = (struct ML_AGG_Matrix_Context *) obj;
   Amat    = (ML_Operator *) context->Amat;
   nRows   = Amat->outvec_leng;
   getrow_obj = Amat->getrow;
   if (getrow_obj->ML_id == ML_EXTERNAL) getrowfunc = getrow_obj->external;
   else                                  getrowfunc = getrow_obj->internal;
   if ( getrowfunc == NULL )
   {
      printf("ML_AGG_DD_Getrow ERROR : null getrowfunc.\n");
      exit(-1);
   }

   if ( alloc_space > 0 )
   {
      local_ind = (int *)    malloc( alloc_space * sizeof(int));
      local_val = (double *) malloc( alloc_space * sizeof(double));
   }
   status = getrowfunc(Amat->data, 1, rowlist, alloc_space, local_ind,
                       local_val, rowcnt);
   if ( status == 0 ) 
   {
      free( local_ind );
      free( local_val );
      return 0;
   }
   count = 0;
   for ( i = 0; i < (*rowcnt); i++ )
   {
      if ( local_ind[i] < nRows ) 
      {
         col_ind[count] = local_ind[i];
         col_val[count++] = local_val[i];
      }
   }
   (*rowcnt) = count;
   free( local_ind );
   free( local_val );
   return 1;
}  

/* ************************************************************************* */
/* extract diagonal                                                          */
/* ------------------------------------------------------------------------- */

int ML_AGG_Extract_Diag(ML_Operator *Amat, double *diagonal)
{
   int           i, j, m, max_row_nnz=100, nRows, *col_ind;
   double        *col_val;
   int           (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   ML_GetrowFunc *getrow_obj;

   nRows   = Amat->outvec_leng;
   getrow_obj = Amat->getrow;
   if (getrow_obj->ML_id == ML_EXTERNAL) getrowfunc = getrow_obj->external;
   else                                  getrowfunc = getrow_obj->internal;
   if ( getrowfunc == NULL )
   {
      printf("ML_AGG_Extract_Diag ERROR : null getrowfunc.\n");
      exit(-1);
   }
   col_ind = (int *)    malloc( max_row_nnz * sizeof(int) );
   col_val = (double *) malloc( max_row_nnz * sizeof(double) );

   for ( i = 0; i < nRows; i++ )
   {
      while (getrowfunc(Amat->data,1,&i,max_row_nnz,col_ind,col_val,&m)== 0)
      {
         free( col_ind );
         free( col_val );
         max_row_nnz = max_row_nnz * 2 + 1;
         col_ind = (int *)    malloc( max_row_nnz * sizeof(int) );
         col_val = (double *) malloc( max_row_nnz * sizeof(double) );
      }
      for (j = 0; j < m; j++) if (col_ind[j] == i) diagonal[i] = col_val[j];
   }
   free( col_ind );
   free( col_val );

   return 1;
}

/* ************************************************************************* */
/* destroy aggregate matrix context                                          */
/* ------------------------------------------------------------------------- */

void ML_AGG_Matrix_Context_Clean(void *data)
{
   struct ML_AGG_Matrix_Context *context;

   context = (struct ML_AGG_Matrix_Context *) data;
   free (context);
}

/* ************************************************************************* */
/* solve local subproblem using smoothed aggregation                         */
/* ------------------------------------------------------------------------- */

int ML_AGG_DD_Solve(void *data, int leng1, double *outvec, int leng2, 
                    double *invec)
{
   ML  *ml = (ML *) data;
   ML_Solve_MGV( ml, invec, outvec );
   return 1;
}

/* ************************************************************************* */
/* solve local subproblem using smoothed aggregation                         */
/* ------------------------------------------------------------------------- */

int ML_AGG_Extract_Matrix(ML_Operator *mat, int *ncols, int **cols,
                          double ***vals)
{
   int           i, j, nnz, local_nrows, *col_ind, row_size, max_size;
   int           index, local_ncols, *local_cols;
   double        *col_val, **local_vals;
   int           (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   ML_GetrowFunc *getrow_obj;

   local_nrows = mat->outvec_leng;
   getrow_obj = mat->getrow;
   if (getrow_obj->ML_id==ML_EXTERNAL) getrowfunc = getrow_obj->external;
   else                                getrowfunc = getrow_obj->internal;

   /* ----------------------------------------------------------------- */
   /* compute number of nonzeros                                        */
   /* ----------------------------------------------------------------- */

   max_size = 3;
   col_ind = (int *)    malloc( max_size * sizeof(int) );
   col_val = (double *) malloc( max_size * sizeof(double) );
   nnz = 0;
   for ( i = 0; i < local_nrows; i++ )
   {
      while (getrowfunc(mat->data,1,&i,max_size,col_ind,col_val,&row_size)==0)
      {
         free( col_ind );
         free( col_val );
         max_size = max_size *2 + 1;
         col_ind = (int *)    malloc( max_size * sizeof(int) );
         col_val = (double *) malloc( max_size * sizeof(double) );
      }
      nnz += row_size;
      if ( row_size > max_size ) max_size = row_size;
   }
   free( col_ind );
   free( col_val );

   /* ----------------------------------------------------------------- */
   /* extract matrix                                                    */
   /* ----------------------------------------------------------------- */

   col_ind = (int *) malloc( nnz * sizeof(int));
   col_val = (double *) malloc( nnz * sizeof(double));
   nnz = 0;
   for ( i = 0; i < local_nrows; i++ )
   {
      getrowfunc(mat->data,1,&i,max_size,&col_ind[nnz],&col_val[nnz],&row_size);
      nnz += row_size;
   }

   /* ----------------------------------------------------------------- */
   /* find number of distinct nonzero columns                           */
   /* ----------------------------------------------------------------- */

   ML_az_sort( col_ind, nnz, NULL, NULL );
   local_ncols = 0;
   for ( i = 1; i < nnz; i++ )
   {
      if ( col_ind[i] != col_ind[local_ncols] ) 
         col_ind[++local_ncols] = col_ind[i];
   }
   local_ncols++;
   local_cols = (int *) malloc(local_ncols * sizeof(int));
   for ( i = 0; i < local_ncols; i++ ) local_cols[i] = col_ind[i];
   free( col_ind );
   free( col_val );

   /* ----------------------------------------------------------------- */
   /* fill in the matrix                                                */
   /* ----------------------------------------------------------------- */

   local_vals = (double **) malloc(local_ncols * sizeof(double*));
   for ( i = 0; i < local_ncols; i++ )
   { 
      local_vals[i] = (double *) malloc(local_nrows * sizeof(double));
      for ( j = 0; j < local_nrows; j++ ) local_vals[i][j] = 0.0;
   }

   col_ind = (int *)    malloc( max_size * sizeof(int));
   col_val = (double *) malloc( max_size * sizeof(double));
   for ( i = 0; i < local_nrows; i++ )
   {
      getrowfunc(mat->data,1,&i,max_size,col_ind,col_val,&row_size);
      for ( j = 0; j < row_size; j++ )
      {
         index = ML_sorted_search( col_ind[j], local_ncols, local_cols);
         if ( index >= 0 ) local_vals[index][i] = col_val[j];
      }
   }
   free( col_ind );
   free( col_val );

   (*ncols) = local_ncols;
   (*cols)  = local_cols;
   (*vals)  = local_vals;
   return 1;
}

/* ************************************************************************* */
/* generate smooth prolongator for 2-level DD method                         */
/* ------------------------------------------------------------------------- */

int ML_AGG_Gen_DDProlongator2(ML *ml,int level, int clevel, void *data,
                             ML_Aggregate *ag)
{
   int          i, k, Nfine, nbytes, newNlevels, newClevel;
   int          *new_ia, *new_ja;
   double       *new_val, omega, *diagonal;
   ML           *newml;
   ML_Operator  *Amat, *tentP, *AGGsmoother;
   ML_Aggregate *newag;
   ML_Aggregate_Comm            *aggr_comm;
   struct ML_CSR_MSRdata        *csr_data;
   struct ML_AGG_Matrix_Context widget, *context;

#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   /* ----------------------------------------------------------------- */
   /* coarsen local smoothed aggregation method                         */
   /* ----------------------------------------------------------------- */
   widget.near_bdry = NULL;
   Amat  = (ML_Operator *) data;
   Nfine = Amat->outvec_leng;
   omega = ag->smoothP_damping_factor;

   /* ----------------------------------------------------------------- */
   /* setup local smoothed aggregation method                           */
   /* ----------------------------------------------------------------- */

   if ( omega != 0.0 )
   {
      newNlevels = 15;
      ML_Create(&newml, newNlevels);
      ML_Set_OutputLevel(newml, 0);
      ML_Set_ResidualOutputFrequency(newml, 0);
      ML_Set_Comm_MyRank(newml, 0);
      ML_Set_Comm_Nprocs(newml, 1);
      nbytes = sizeof(struct ML_AGG_Matrix_Context);
      context = (struct ML_AGG_Matrix_Context *) malloc( nbytes );
      context->Amat = Amat;
      context->near_bdry = NULL;
      ML_Init_Amatrix(newml, newNlevels-1, Nfine,  Nfine, (void *) context);
      ML_Set_Amatrix_Matvec(newml,  newNlevels-1, ML_AGG_DD_Matvec);
      newml->Amat[newNlevels-1].data_destroy = ML_AGG_Matrix_Context_Clean;
      newml->Amat[newNlevels-1].N_nonzeros = 5 * Nfine;
      ML_Set_Amatrix_Getrow(newml, newNlevels-1, ML_AGG_DD_Getrow, NULL, Nfine);
      diagonal = (double *) malloc(Nfine * sizeof(double));
      ML_AGG_Extract_Diag(Amat, diagonal);
      ML_Set_Amatrix_Diag( newml, newNlevels-1, Nfine, diagonal);
      free( diagonal );
      ML_Aggregate_Create( &newag );
      ML_Aggregate_Set_OutputLevel( newag, 0 );
      ML_Aggregate_Set_CoarsenScheme_Uncoupled( newag );
      ML_Aggregate_Set_MaxCoarseSize( newag, 50 );
      ML_Aggregate_Set_PSmootherType( newag, 0 );
      newClevel = ML_Gen_MGHierarchy_UsingAggregation(newml, newNlevels-1,
                                     ML_DECREASING, newag);
      newClevel = newNlevels - newClevel;

      for (k = newNlevels-1; k > newClevel; k--)
      {
         ML_Gen_Smoother_SymGaussSeidel(newml, k, ML_PRESMOOTHER, 1, 1.);
         ML_Gen_Smoother_SymGaussSeidel(newml, k, ML_POSTSMOOTHER, 1, 1.);
      }
      ML_Gen_CoarseSolverSuperLU( newml, newClevel );
      ML_Gen_Solver(newml, ML_MGV, newNlevels-1, newClevel);
      ML_Aggregate_Destroy( &newag );
   }

   /* ----------------------------------------------------------------- */
   /* compute tentP (local ml P applied to 1)                           */
   /* ----------------------------------------------------------------- */

   nbytes = ( Nfine + 1 ) * sizeof(int);
   ML_memory_alloc((void**)&(new_ia), nbytes, "AD1");
   nbytes = Nfine * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), nbytes, "AD2");
   nbytes = Nfine * sizeof(double);
   ML_memory_alloc((void**)&(new_val), nbytes, "AD3");
   for (i = 0; i <= Nfine; i++) new_ia[i] = i;
   for (i = 0; i < Nfine; i++) new_ja[i] = 0;
/*
   norm = sqrt((double) Nfine);
norm = 1.0;
   for (i = 0; i < Nfine; i++) new_val[i] = 1.0 / norm;
*/
   ML_memory_alloc((void**) &csr_data,sizeof(struct ML_CSR_MSRdata),"AVP");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
/*
   tentP = &(ml->Pmat[clevel]);
*/
tentP = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(tentP,1,Nfine,ML_EMPTY,csr_data,Nfine,NULL,0);
   tentP->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm), "AD4");
   aggr_comm->comm = ml->comm;
   aggr_comm->N_send_neighbors = 0;
   aggr_comm->N_recv_neighbors = 0;
   aggr_comm->send_neighbors = NULL;
   aggr_comm->recv_neighbors = NULL;
   aggr_comm->send_leng = NULL;
   aggr_comm->recv_leng = NULL;
   aggr_comm->send_list = NULL;
   aggr_comm->local_nrows = 1;
   ML_CommInfoOP_Generate( &(tentP->getrow->pre_comm),
                           ML_Aggregate_ExchangeBdry, aggr_comm, ml->comm, 1, 0);
   ML_Operator_Set_Getrow(tentP, ML_EXTERNAL, Nfine, CSR_getrows);
   ML_Operator_Set_ApplyFunc(tentP, ML_INTERNAL, CSR_matvec);
   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));

/*###*/
   widget.Amat = Amat;
   widget.aggr_info = ag->aggr_info[level];
   AGGsmoother = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(AGGsmoother, widget.Amat->invec_leng,
                        widget.Amat->outvec_leng, ML_EXTERNAL,&widget,
                        widget.Amat->matvec->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(AGGsmoother, ML_EXTERNAL,
                          widget.Amat->getrow->Nrows, 
                          ML_AGG_Amat_Getrows);
   ML_CommInfoOP_Clone(&(AGGsmoother->getrow->pre_comm),
                          widget.Amat->getrow->pre_comm);

   ML_2matmult(AGGsmoother, tentP, &(ml->Pmat[clevel]) );
if ( ml->comm->ML_mypid == 1 )
ML_Operator_Print(&(ml->Pmat[clevel]), "Pmat2");

   ML_Operator_Destroy(tentP);
   ML_Operator_Destroy(AGGsmoother);
/*###*/

#ifdef ML_TIMING
   ml->Pmat[clevel].build_time =  GetClock() - t0;
   ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
   return 0;
}

/* ************************************************************************* */
/* Compute the DOFs that are on the boundary and those that are right next   */
/* to the boundary.                                                          */
/* ------------------------------------------------------------------------- */
int ML_AGG_Compute_Near_Bdry(ML_Operator *Amatrix, char *near_bdry)
{
  int Nrows, Nghost = 0, allocated = 0, *rowi_col = NULL, rowi_N, count2;
  int i, j, bsize, flag;
  double *dtemp, *rowi_val = NULL, sum;


  Nrows = Amatrix->outvec_leng;

   /* ============================================================= */
   /* Figure out where the Dirichlet points are on the fine grid.   */
   /* ============================================================= */

  if (Amatrix->getrow->pre_comm != NULL)
    Nghost = Amatrix->getrow->pre_comm->total_rcv_length;

  /* near_bdry = (char *) ML_allocate(sizeof(char)*(Nrows+Nghost+1)); */
   dtemp = (double *) ML_allocate(sizeof(double)*(Nrows+Nghost+1));
   if (dtemp == NULL) pr_error("ml_agg_MIS: out of space.\n");

   for (i = 0; i < Nrows+Nghost; i++) dtemp[i] = 0.;

   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      count2 = 0;
      for (j = 0; j < rowi_N; j++) if (rowi_val[j] != 0.) count2++;
      if (count2 <= 1) dtemp[i] = 1.;
   }
  
   /* if one DOF within a node is fixed, mark all the DOFs within node */

   bsize = Amatrix->num_PDEs;
   for (i = 0; i < Nrows/bsize; i++) {
     sum = 0.;
     for (j = 0; j < bsize; j++) {
       sum += dtemp[i*bsize+j];
     }
     if (sum != 0.) {
       for (j = 0; j < bsize; j++) dtemp[i*bsize+j] = 1.;
     }
   }
      

   
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,Amatrix->outvec_leng,
                    Amatrix->comm, ML_OVERWRITE);
   for (i = 0; i < Nrows+Nghost; i++) {
      if (dtemp[i] == 1.) near_bdry[i] = 'T';
      else near_bdry[i] = 'F';
   }

   /* Figure out who touches a Dirichlet point */

   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      count2 = 0;
      for (j = 0; j < rowi_N; j++) if (dtemp[rowi_col[j]] != 0.) count2++;
      if (count2 != 0) near_bdry[i] = 'T';
   }

   for (i = 0; i < Nrows/bsize; i++) {
     flag = 0;
     for (j = 0; j < bsize; j++) {
       if (near_bdry[i*bsize+j] == 'T') flag = 1;
     }
     if (flag == 1) {
        for (j = 0; j < bsize; j++) {
	  near_bdry[i*bsize+j] = 'T';
        }
     }
   }

   
   
   free(rowi_col); free(rowi_val);
   rowi_col = NULL; rowi_val = NULL;
   allocated = 0; 

   ML_free(dtemp);

   return 0;
}
