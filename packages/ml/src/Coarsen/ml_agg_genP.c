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

#ifdef GEOMETRIC_2D
extern double lambda_max,lambda_min, sten_a,sten_b,sten_c, sten_d,sten_e,sten_f,
  sten_g,sten_h,sten_i;
extern int mls_or_gs, mls_order;
#endif

#include <math.h>
#include <stdlib.h>
#include "ml_struct.h"
#include "ml_op_utils.h"
#include "ml_agg_genP.h"
#include "ml_memory.h"

extern int ML_Anasazi_Get_FieldOfValuesBox_Interface(ML_Operator * Amat,
						     struct ML_Field_Of_Values * fov );
extern int ML_Anasazi_Get_FieldOfValuesBoxNonScaled_Interface(ML_Operator * Amat,
						     struct ML_Field_Of_Values * fov );
extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
					       int MaxIters, double Tolerance,
					       int IsProblemSymmetric,
					       int UseDiagonalScaling,
					       double * LambdaMax );

/* ************************************************************************* */
/* wrapper function as smoother                                              */
/* ------------------------------------------------------------------------- */

int ML_AGG_Smoother_Wrapper(void *obj, int leng1, double *outvec, int leng2,
                            double *invec)
{
   ML *ml;
   ml = (ML *) obj;
   ML_Iterate( ml, outvec, invec );
   ML_avoid_unused_param( (void *) &leng1);
   ML_avoid_unused_param( (void *) &leng2);
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

   if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel())
   {
     printf("Entering ML_Gen_MGHierarchy_UsingAggregation\n");
     fflush(stdout);
   }
   ML_memory_check("L%d:gen_hier start",start);

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
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag < ML_Get_PrintLevel()) 
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
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag < ML_Get_PrintLevel()) 
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
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag < ML_Get_PrintLevel()) 
      ML_Aggregate_Print_Complexity( ml_ag );
   idata = ML_gmax_int(idata, ml->comm);

   if ( ag == NULL ) ML_Aggregate_Destroy( &ml_ag );
   ML_memory_check("gen hierarchy end");

   return(level);
}

/* ************************************************************************* */
/* generate multilevel hierarchy given a subroutine for generating           */
/* prolongation operators (version 2 : with ML_Aggregate given)              */
/* ------------------------------------------------------------------------- */

int ML_Gen_MGHierarchy(ML *ml, int fine_level,
        int (*next_level)(ML *, int,  void *),
        int (*user_gen_prolongator)(ML *, int, int, void *),
        void *data, int internal_or_external, ML_Aggregate *ag)
{
   int level, next, flag, count=1;
   int i, j, bail_flag, N_input_vector;
   ML_Operator *Pmat;
   ML_CommInfoOP *getrow_comm; 
#ifdef ML_TIMING
   double t0;
#endif
#ifdef GEOMETRIC_2D
   struct  ML_CSR_MSRdata *csr_data;
   int m,i,j,row,k, bsize;
   double s;
#endif

#ifdef MARZIO
#define USE_ATandT
#endif
#ifdef USE_AT
   char str[80];
   double dtemp;
   struct ML_Field_Of_Values * fov;
#endif

   if (ag->nullspace_corrupted == ML_YES) {
     printf("Can not reuse aggregate object when the fine grid operator\n");
     printf("has a nontrivial null space. It is possible to keep \n");
     printf("tentative prolongator within smoothed aggregation by\n");
     printf("invoking ML_Aggregate_Set_Reuse(...) before hierarchy\n");
     printf("generation and then on subsequent hierarchy generations use\n");
     printf("ML_Gen_MGHierarchy_UsingSmoothedAggr_ReuseExistingAgg().\n");

     exit(-1);
   }

   ml->ML_finest_level = fine_level;
   level = fine_level;
   next  = next_level(ml, level, ag);

   while (next >= 0) 
   {
      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("ML_Gen_MGHierarchy : applying coarsening \n");
#ifdef USE_AT
      /*
        Here is the idea:
        1) build an interpolation operator based on A^T. Save the
           tentative prolongator.
        3) transpose it for the restriction as usual.
        4) clean out the interpolation. Now build the interpolation
           operator using the existing tentative prolongator and
           the symmetry option.
      */
      if (ag->Restriction_smoothagg_transpose == ML_TRUE ) {
        ml->symmetrize_matrix = ML_TRUE;
        ag->keep_P_tentative = ML_YES;
        ag->use_transpose = ML_TRUE;
	if (1 == 1) {  
	  ml->symmetrize_matrix = ML_FALSE;
	  /*	  printf(" what is this %e\n",ml->Amat[level].lambda_max); */
	  fov = (struct ML_Field_Of_Values * )(ag->field_of_values);
	  dtemp = fov->R_coeff[0] + fov->eta * fov->R_coeff[1] + pow(fov->eta,2) * fov->R_coeff[2];
	  ml->Amat[level].lambda_max = dtemp;
	}
      }
#endif

      if (internal_or_external == ML_INTERNAL)
      {
         flag = (*user_gen_prolongator)(ml, level, next,(void *)ag);
      }
      else 
      {
         flag = (*user_gen_prolongator)(ml, level, next, (void *)ag);
      }
      if (flag < 0) break;
      ML_memory_check("L%d: prolongator end",level);

      /* Now check to make sure prolongator has zero columns. */
      Pmat = ml->Pmat+next;
      bail_flag = 0;
      N_input_vector = Pmat->invec_leng;
      getrow_comm = Pmat->getrow->pre_comm;
      if ( getrow_comm != NULL)
      {
         for (i = 0; i < getrow_comm->N_neighbors; i++) {
            for (j = 0; j < getrow_comm->neighbors[i].N_send; j++) {
               if (getrow_comm->neighbors[i].send_list[j] >= N_input_vector) {
                  bail_flag = 1;
               }
            }
         }
      }
      /* If check has failed on any processor, clean up current level & break
         from main loop. */
      ML_gsum_scalar_int(&bail_flag,&j,ml->comm);
      if (bail_flag)
      {
         if (Pmat->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
            printf("(%d) In ML_Gen_MGHierarchy: Bailing from AMG hierarchy build on level %d, where fine level = %d ........\n",
                   Pmat->comm->ML_mypid,level,fine_level);
            fflush(stdout);
         }
         if (ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
            printf("(%d) In ML_Gen_MGHierarchy: Nlevels = %d fine_level = %d  coarsest_level = %d\n",
               ml->comm->ML_mypid,fine_level-count+1,fine_level,count);
            fflush(stdout);
         }
         break; /* from main loop */
 
      }
      /* end of check */
      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("ML_Gen_MGHierarchy : applying coarsening \n");
      ML_Gen_Restrictor_TransP(ml, level, next);
      ML_Operator_ChangeToSinglePrecision(&(ml->Pmat[next]));
      ML_memory_check("L%d: TransP end",level);

#ifdef USE_AT
      if (ag->Restriction_smoothagg_transpose == ML_TRUE ) {
        ag->use_transpose = ML_FALSE;
        ML_Operator_Clean( &(ml->Pmat[next]) );
        ML_Operator_Init( &(ml->Pmat[next]), ml->comm );
	ML_Operator_Set_1Levels(&(ml->Pmat[next]), &(ml->SingleLevel[next]), NULL);
	ML_Operator_Set_BdryPts(&(ml->Pmat[next]), NULL);
	sprintf(str,"Pmat_%d",next); ML_Operator_Set_Label( &(ml->Pmat[next]),str);

        if (internal_or_external == ML_INTERNAL)
          flag = (*user_gen_prolongator)(ml, level, next,(void *)ag);
        else
          flag = (*user_gen_prolongator)(ml, level, next, (void *)ag);
      }

#endif

      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("ML_Gen_MGHierarchy : Gen_RAP\n");

#ifdef ML_TIMING
      t0 = GetClock();
#endif
      ML_Gen_AmatrixRAP(ml, level, next);

      ML_Operator_ImplicitTranspose(&(ml->Rmat[level]),
      			    &(ml->Pmat[next]), ML_TRUE);

      ML_memory_check("L%d: RAP end",level);

#ifdef GEOMETRIC_2D
   csr_data = ml->Amat[next].data;
   m = (int) sqrt( ((double)ml->Amat[next].invec_leng) + .00001); k = m*m+1;
   csr_data->columns[0] = k;
   bsize = csr_data->columns[m*m];
   for (j = 0; j < m; j++) {
     for (i = 0; i < m; i++) {
       row = i +j*m;
       csr_data->columns[k]= row+1;
       if ((row)%m !=     m-1                       ) csr_data->values[k++] = sten_f;
       csr_data->columns[k]= row-1;
       if ((row)%m !=       0                       ) csr_data->values[k++] = sten_d;
       csr_data->columns[k]= row+m-1; 
       if (((row/(m))%m != m-1) &&((row)%m != 0)    ) csr_data->values[k++] = sten_a;
       csr_data->columns[k]= row+m;
       if ((row/(m))%m != m-1                       ) csr_data->values[k++] = sten_b;
       csr_data->columns[k]= row+m+1;
       if (((row/(m))%m != m-1) &&((row)%m !=   m-1)) csr_data->values[k++] = sten_c;
       csr_data->columns[k]= row-m-1;
       if (((row/(m))%m !=   0) &&((row)%m !=     0)) csr_data->values[k++] = sten_g;
       csr_data->columns[k]= row-m;
       if ((row/(m))%m !=   0                       ) csr_data->values[k++] = sten_h;
       csr_data->columns[k]= row-m+1;
       if (((row/(m))%m !=   0) && ((row)%m !=  m-1)) csr_data->values[k++] = sten_i;
       csr_data->values[row]     = sten_e;
       csr_data->columns[row+1] = k;
     }
   }
   if (k > bsize) {
     printf("problemssssssssssssssssssss %d %d\n",k,bsize);
     exit(1);
   }
#endif

#ifdef ML_TIMING
      t0 = GetClock() - t0;
      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("RAP time for level %2d = %e\n", level, t0);
#endif

      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("ML_Gen_MGHierarchy : Gen_RAP done\n");

      level = next;
      next  = next_level(ml, next, ag);
      count++;
   }
   return(count);
}

/* ************************************************************************* */
/* generate smooth prolongator                                               */
/* ------------------------------------------------------------------------- */

int ML_AGG_Gen_Prolongator(ML *ml,int level, int clevel, void *data)
{
   int         Ncoarse, Nfine, gNfine, gNcoarse, jj;
   double      max_eigen = -1.;
   ML_Operator *Amat, *Pmatrix = NULL, *AGGsmoother = NULL;
   ML_Operator **prev_P_tentatives;
   struct      ML_AGG_Matrix_Context widget;
   ML_Krylov   *kdata;
   ML_Operator *t2 = NULL, *t3 = NULL;
   ML_Aggregate * ag = (ML_Aggregate *) data;
   struct MLSthing *mls_widget = NULL;
   ML_Operator *blockMat = NULL, *Ptemp;

#ifdef GEOMETRIC_2D
   int nx, nxcoarse, ii, coarse_me, fine_me, *rowptr, *bindx, k,free_ptr,start;
   int i,j, end;
   double *val;
   struct  ML_CSR_MSRdata *csr_data;
#endif

#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel())
   {
     printf("Entering ML_AGG_Gen_Prolongator\n");
     fflush(stdout);
   }

   Amat = &(ml->Amat[level]);

   if (Amat->num_PDEs < ag->num_PDE_eqns) Amat->num_PDEs = ag->num_PDE_eqns;
   if (ag->block_scaled_SA == 1) {
     /* 
         Create block scaled and compute its eigenvalues
	 a) if the user has requested it, save this Amat into
	    the aggregate data structure.
      */
     mls_widget = ML_Smoother_Create_MLS();
     ML_Gen_BlockScaledMatrix_with_Eigenvalues(Amat, -1, NULL,
					       &blockMat, mls_widget);
     max_eigen = blockMat->lambda_max;
   }
   else    max_eigen = Amat->lambda_max;
   
   widget.near_bdry = NULL;
   Amat->num_PDEs = ag->num_PDE_eqns;
   prev_P_tentatives = ag->P_tentative;

#ifdef EXTREME_DEBUG
   printf("### %e %e\n",  ag->smoothP_damping_factor, max_eigen);
#endif

   /*
   widget.near_bdry = (char *) ML_allocate(sizeof(char)*Amat->outvec_leng);
   ML_AGG_Compute_Near_Bdry(Amat, widget.near_bdry);
   */

   Nfine    = Amat->outvec_leng;
#ifdef GEOMETRIC_2D
   nx = (int) floor( sqrt( (double) Nfine) + .00001 );
   printf("in here %d %d\n",Nfine,nx);
   nxcoarse = (nx-1)/2;
   rowptr = (int *) ML_allocate(sizeof(int)*(Nfine+1));
   bindx  = (int *) ML_allocate(sizeof(int)*(Nfine*4+1));
   val    = (double *) ML_allocate(sizeof(double)*Nfine*4);
   for (ii = 0; ii <= Nfine; ii++) rowptr[ii] = 4*ii;
   for (ii = 0; ii <= 4*Nfine; ii++) bindx[ii] = -1;
   for (ii = 0; ii < nxcoarse; ii++) {
      for (jj = 0; jj < nxcoarse; jj++) {
	coarse_me = jj*nxcoarse+ii;
	fine_me = (2*jj+1)*nx + 2*ii+1;
	if (jj == 0) {
	  if (ii == 0) {
	    bindx[rowptr[fine_me-nx-1]] = coarse_me;
	    val[  rowptr[fine_me-nx-1]] = .25;
	  }
	  bindx[rowptr[fine_me-nx]] = coarse_me;
	  val[  rowptr[fine_me-nx]] = .5;
	  bindx[rowptr[fine_me-nx+1]] = coarse_me;
	  val[  rowptr[fine_me-nx+1]] = .25;
	  if (ii != nxcoarse-1) {
	    bindx[rowptr[fine_me-nx+1]+1] = coarse_me+1;
	    val[  rowptr[fine_me-nx+1]+1] = .25;
	  }
	}
	if (ii == 0) {
	  bindx[rowptr[fine_me-1]] = coarse_me;
	  val[  rowptr[fine_me-1]] = .5;
	  bindx[rowptr[fine_me-1+nx]] = coarse_me;
	  val[  rowptr[fine_me-1+nx]] = .25;
	  if (jj != nxcoarse-1) {
	    bindx[rowptr[fine_me-1+nx]+1] = coarse_me+nxcoarse;
	    val[  rowptr[fine_me-1+nx]+1] = .25;
	  }
	}
	bindx[rowptr[fine_me]] = coarse_me;
	val[  rowptr[fine_me]] = 1.;
	bindx[rowptr[fine_me+1]] = coarse_me;
	val[  rowptr[fine_me+1]] = .5;
        if (ii != nxcoarse-1) {
	  bindx[rowptr[fine_me+1]+1] = coarse_me + 1;
	  val[  rowptr[fine_me+1]+1] = .5;
	}
	bindx[rowptr[fine_me+nx]] = coarse_me;
	val[  rowptr[fine_me+nx]] = .5;
        if (jj != nxcoarse-1) {
	  bindx[rowptr[fine_me+nx]+1] = coarse_me + nxcoarse;
	  val[  rowptr[fine_me+nx]+1] = .5;
	}
	bindx[rowptr[fine_me+1+nx]] = coarse_me;
	val[  rowptr[fine_me+1+nx]] = .25; 
	k = 1;
	if (ii != nxcoarse-1) {
	  bindx[rowptr[fine_me+1+nx]+k] = coarse_me+1;
	  val[  rowptr[fine_me+1+nx]+k] = .25;
	  k++;
	}
	if (jj != nxcoarse-1) {
	  bindx[rowptr[fine_me+1+nx]+k] = coarse_me+nxcoarse;
	  val[  rowptr[fine_me+1+nx]+k] = .25;
	  k++;
	}
	if ((ii != nxcoarse-1) && (jj != nxcoarse-1)) {
	  bindx[rowptr[fine_me+1+nx]+k] = coarse_me+nxcoarse+1;
	  val[  rowptr[fine_me+1+nx]+k] = .25;
	  k++;
	}

      }
   }
   csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata)); 
   ml->Pmat[clevel].data = csr_data;
   csr_data->rowptr = rowptr;
   csr_data->values = val;
   csr_data->columns = bindx;

   /* compress bindx[] = '1  out matrix */

   free_ptr = rowptr[0];
   start = rowptr[0];
   for (i = 0; i < Nfine; i++) {
     rowptr[i] = rowptr[i+1] - rowptr[i];
   }
   for (i = 0; i < Nfine; i++) {
     end = start + rowptr[i] - 1;
     for (j = start; j <= end; j++) {
       if (bindx[j] != -1) {
	 val[free_ptr] = val[j];
	 bindx[free_ptr] = bindx[j];
	 free_ptr++;
       }
       else rowptr[i]--;
     }
     start = end+1;
   }

   start = 0;
   for (i = 0; i <= Nfine; i++) {
     j = rowptr[i];
     rowptr[i] = start;
     start += j;
   }

   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));
   ML_Operator_Set_ApplyFuncData(&(ml->Pmat[clevel]), nxcoarse*nxcoarse,
				 Nfine, ML_EMPTY, csr_data,
				 Nfine, NULL, 0);
   ML_Operator_Set_Getrow(&(ml->Pmat[clevel]), ML_EXTERNAL, 
			  Nfine, CSR_getrows);
   ml->Pmat[clevel].max_nz_per_row = 4;
   ml->Pmat[clevel].N_nonzeros = bindx[Nfine];
   ML_Operator_Set_ApplyFunc (&(ml->Pmat[clevel]), ML_INTERNAL, CSR_matvec);
   if (mls_widget != NULL) ML_Smoother_Destroy_MLS(mls_widget);
   return 1;
#endif
   gNfine   = ML_Comm_GsumInt( ml->comm, Nfine);
   ML_Aggregate_Set_CurrentLevel( ag, level );

   if (ag->smoothP_damping_factor != 0.0)
   {
     if ((ag->keep_P_tentative == ML_YES) && (prev_P_tentatives != NULL) &&
	     (prev_P_tentatives[clevel] != NULL))
     {
       Pmatrix = prev_P_tentatives[clevel];
       Ncoarse = Pmatrix->invec_leng;
     }
	 else if (ML_Aggregate_Get_Flag_SmoothExistingTentativeP(ag) == ML_YES)
     {
       Pmatrix = ML_Operator_halfClone( &(ml->Pmat[clevel]) );
       /* ml->Pmat[clevel] is destroyed first.  Half cloning assumes that
          everything associated with matrix->getrow is destroyed with the
          original operator.  In this case, the getrow object should be
          destroyed with the clone.  Hence, we lie and say Pmatrix is not
          the result of a clone....*/
       Pmatrix->halfclone = ML_FALSE;

       Pmatrix->data_destroy = ml->Pmat[clevel].data_destroy;
       ml->Pmat[clevel].data_destroy = NULL;
       ml->Pmat[clevel].data = NULL;

       ML_memory_free( (void**)&(Pmatrix->matvec) );
       Pmatrix->matvec = ml->Pmat[clevel].matvec;
       ml->Pmat[clevel].matvec = NULL;

       ML_memory_free( (void**)&(Pmatrix->getrow) );
       Pmatrix->getrow = ml->Pmat[clevel].getrow;
       ml->Pmat[clevel].getrow = NULL;
       Pmatrix->label = ml->Pmat[clevel].label;
       ml->Pmat[clevel].label = NULL;

       ML_Operator_Clean(&(ml->Pmat[clevel]));
       ML_memory_alloc((void**)&(ml->Pmat[clevel].getrow),
                       sizeof(ML_GetrowFunc),"OF2");

       ml->Pmat[clevel].matvec = Pmatrix->matvec;
       Pmatrix->matvec = NULL;
       ml->Pmat[clevel].label = Pmatrix->label;
       Pmatrix->label = NULL;
       Ncoarse = Pmatrix->invec_leng;
     }
     else {
       Pmatrix = ML_Operator_Create(ml->comm);
       Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&Pmatrix,ml->comm);
     }
   }
   else
   {
     Pmatrix = &(ml->Pmat[clevel]);
     if (Pmatrix->invec_leng == 0)
       Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&Pmatrix,ml->comm);
     else Ncoarse = Pmatrix->invec_leng;
   }
   gNcoarse = ML_Comm_GsumInt( ml->comm, Ncoarse);
   if ( gNcoarse == 0 || ((1.0*gNfine)/(1.0*gNcoarse+0.1) < 1.05) )
   {
     if (( Pmatrix != NULL ) && (ag->smoothP_damping_factor != 0.0))
     {
		if (ML_Aggregate_Get_Flag_SmoothExistingTentativeP(ag) == ML_YES)
        {
           ml->Pmat[clevel].data_destroy = Pmatrix->data_destroy;
           Pmatrix->data_destroy = NULL;

           ML_memory_free( (void**)&(ml->Pmat[clevel].getrow) );
           ml->Pmat[clevel].getrow = Pmatrix->getrow;
           Pmatrix->getrow = NULL;
           
           ml->Pmat[clevel].data = Pmatrix->data;
           Pmatrix->data = NULL;
        }
        ML_Operator_Destroy(&Pmatrix);
     }
     if (mls_widget != NULL) ML_Smoother_Destroy_MLS(mls_widget);
     return -1;
   }

#ifdef ML_NONSYMMETRIC
   max_eigen = ML_Operator_MaxNorm(Amat, ML_TRUE);
   widget.omega  = ag->smoothP_damping_factor / max_eigen;
   ml->spectral_radius[level] = max_eigen;
#else

   if ( ag->smoothP_damping_factor != 0.0 )
   {
#ifdef SYMMETRIZE
     t2 = ML_Operator_Create(Amat->comm);
     ML_Operator_Transpose_byrow(Amat,t2);
     t3 = ML_Operator_Create(Amat->comm);
     ML_Operator_Add(Amat,t2,t3,ML_CSR_MATRIX,1.);
     max_eigen = t3->lambda_max;
#else
#ifdef USE_AT
     if (ag->use_transpose == ML_TRUE) {
       t3 = ML_Operator_Create(Amat->comm);
       ML_Operator_Transpose_byrow(Amat,t3);
       t3->lambda_max = 1.5;
       max_eigen = t3->lambda_max;
       /*       printf("A^T lambda_max is settttttt\n");*/
     }
     else
#endif
     if (ml->symmetrize_matrix == ML_TRUE) {
       t2 = ML_Operator_Create(Amat->comm);
       ML_Operator_Transpose_byrow(Amat,t2);
       t3 = ML_Operator_Create(Amat->comm);
       ML_Operator_Add(Amat,t2,t3,ML_CSR_MATRIX,1.);
       max_eigen = t3->lambda_max;
     }
#endif

     if ((max_eigen < -666.) && (max_eigen > -667)) {

       switch( ag->spectral_radius_scheme ) {

       case 1:  /* compute it using CG */
	 
         kdata = ML_Krylov_Create( ml->comm );
         ML_Krylov_Set_PrintFreq( kdata, 0 );
         ML_Krylov_Set_ComputeEigenvalues( kdata );
#ifdef SYMMETRIZE
         ML_Krylov_Set_Amatrix(kdata, t3);
#else
#ifdef USE_AT
         if (ag->use_transpose ==ML_TRUE) ML_Krylov_Set_Amatrix(kdata, t3);
         else
#endif
	 if (ml->symmetrize_matrix ==ML_TRUE) ML_Krylov_Set_Amatrix(kdata, t3);
         else ML_Krylov_Set_Amatrix(kdata, Amat);
#endif
         ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
         max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);
	 Amat->lambda_max = max_eigen; 
	 Amat->lambda_min = kdata->ML_eigen_min; 
         ML_Krylov_Destroy( &kdata );
         if ( max_eigen <= 0.0 )
         {
            printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
            max_eigen = 1.0;
         }
	 /* I use a print statement after computations
         if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
            printf("Gen_Prolongator : max eigen = %e \n", max_eigen);
	 */
         widget.omega  = ag->smoothP_damping_factor / max_eigen;
         ml->spectral_radius[level] = max_eigen;

	 break;

       case 2: /* Use Anasazi */
	 ML_Anasazi_Get_SpectralNorm_Anasazi( Amat, 10, 1e-5,
					      ML_FALSE, ML_TRUE, &max_eigen);
	 Amat->lambda_max = max_eigen; 
	 Amat->lambda_min = -12345.6789;
	 if ( max_eigen <= 0.0 ) {
	   printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
	   max_eigen = 1.0;
	 }
	 /*
	   if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
	   printf("Gen_Prolongator : max eigen = %e \n", max_eigen);
	 */
	 widget.omega  = ag->smoothP_damping_factor / max_eigen;
	 ml->spectral_radius[level] = max_eigen;

	 break;

       case 3: /* use ML's power method */
	 kdata = ML_Krylov_Create( ml->comm );
	 ML_Krylov_Set_PrintFreq( kdata, 0 );
	 ML_Krylov_Set_ComputeNonSymEigenvalues( kdata );
	 ML_Krylov_Set_Amatrix(kdata, Amat);
	 ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
	 max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);
	 Amat->lambda_max = max_eigen; 
	 Amat->lambda_min = kdata->ML_eigen_min; 
	 ML_Krylov_Destroy( &kdata );
	 if ( max_eigen <= 0.0 )
	   {
	     printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
	     max_eigen = 1.0;
	   }
         widget.omega  = ag->smoothP_damping_factor / max_eigen;
         ml->spectral_radius[level] = max_eigen;
	 
	 break;

       default: /* using matrix max norm */
         max_eigen = ML_Operator_MaxNorm(Amat, ML_TRUE);
         widget.omega  = ag->smoothP_damping_factor / max_eigen;
         ml->spectral_radius[level] = max_eigen;
	 break;
	 
       }

       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nProlongator/Restriction smoother (level %d) : damping factor = %e\n"
		"Prolongator/Restriction smoother (level %d) : ( = %e / %e)\n\n",
		level, ag->smoothP_damping_factor/ max_eigen,
		level,
		ag->smoothP_damping_factor, max_eigen );
       }
       
     }
     else { /* no need to compute eigenvalue .... we already have it */
       widget.omega  = ag->smoothP_damping_factor / max_eigen;
       ml->spectral_radius[level] = max_eigen;
     }
     
     if ( ml->comm->ML_mypid == 0 && 7 < ML_Get_PrintLevel()) {
       printf("Gen_Prolongator (level %d) : Max eigenvalue = %e\n",
	      ag->cur_level,
	      max_eigen);
     }
     
   }
   else  /* damping fact = 0 ==> no need to compute spectral radius */
   {
      ml->spectral_radius[level] = 1.0;
      widget.omega  = 0.0;

      if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	printf("\nProlongator/Restriction smoother (level %d) : damping factor = 0.0\n\n",
	      level );
      }

   }
   
   
#endif

   if ( ag->smoothP_damping_factor != 0.0 ) {
     widget.drop_tol = ag->drop_tol_for_smoothing;
#ifdef SYMMETRIZE
     widget.Amat   = t3;
#else
#ifdef USE_AT
     if (ag->use_transpose == ML_TRUE) widget.Amat   = t3;
     else
#endif
     if (ml->symmetrize_matrix == ML_TRUE) widget.Amat   = t3;
     else widget.Amat   = &(ml->Amat[level]);
#endif
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

#ifdef EXTREME_DEBUG
     printf("\n###the damping parameter is %e\n\n",widget.omega);
#endif


   if (ag->block_scaled_SA == 1) {
     /* Computed the following:
      *	a) turn off the usual 2 mat mult.
      *	b) Ptemp = A*P 
      *	c) Ptemp = Dinv*Ptemp;
      *	d) do an ML_Operator_Add() with the original P.
      */

     Ptemp = ML_Operator_Create(Amat->comm);
     ML_2matmult(Amat, Pmatrix, Ptemp, ML_CSR_MATRIX );
     ML_AGG_DinvP(Ptemp, mls_widget, Amat->num_PDEs);
     ML_Operator_Add(Pmatrix,Ptemp, &(ml->Pmat[clevel]),ML_CSR_MATRIX,
		     -ag->smoothP_damping_factor / max_eigen);
     ML_Operator_Destroy(&Ptemp);
   }
   else 
     ML_2matmult(AGGsmoother, Pmatrix, &(ml->Pmat[clevel]), ML_CSR_MATRIX );
     
#ifdef SYMMETRIZE
     if (t3 != NULL) ML_Operator_Destroy(&t3);
     if (t2 != NULL) ML_Operator_Destroy(&t2);
#else
     if ((ml->symmetrize_matrix == ML_TRUE) ||
         (ag->use_transpose ==ML_TRUE) ) {
       if (t3 != NULL) ML_Operator_Destroy(&t3);
       if (t2 != NULL) ML_Operator_Destroy(&t2);
     }
#endif
     if (ag->keep_P_tentative == ML_NO)  ML_Operator_Destroy(&Pmatrix);
     else {
       if (prev_P_tentatives == NULL) {
	 ag->P_tentative = ML_Operator_ArrayCreate(ag->max_levels);
         prev_P_tentatives = ag->P_tentative;
	 for (jj = 0; jj < ag->max_levels; jj++) prev_P_tentatives[jj] = NULL;
       }
       prev_P_tentatives[clevel] = Pmatrix;
     }

     ML_Operator_Destroy(&AGGsmoother);
   }
   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));

   if (widget.near_bdry != NULL) ML_free(widget.near_bdry);
   /*
   if (Amat->comm->ML_mypid==0 && ag->print_flag)
   {
      printf("Pe: Total nonzeros = %d (Nrows = %d)\n",
             ml->Pmat[clevel].N_nonzeros, ml->Pmat[clevel].outvec_leng);
   }
   */
   if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel())
   {
     printf("Leaving ML_AGG_Gen_Prolongator\n");
     fflush(stdout);
   }
   if (mls_widget != NULL) ML_Smoother_Destroy_MLS(mls_widget);

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

int ML_AGG_Increment_Level(ML *ml, int current_level,
                           void *data)
{
   int total_size, temp;
   ML_Operator * Amat = &(ml->Amat[current_level]);
   ML_Aggregate * ag = (ML_Aggregate *)data;
     
   if (current_level == ml->ML_num_levels-1) return(-1);

   total_size = Amat->invec_leng;
   ML_gsum_scalar_int(&total_size, &temp, ml->comm);
   if ( total_size <= ag->max_coarse_size ) return(-1);

   return(current_level+1);
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level number */
/* smaller than the fine levels                                              */
/* ------------------------------------------------------------------------- */

int ML_AGG_Decrement_Level(ML *ml, int current_level, void * data)
{
   int total_size, temp;
   ML_Operator * Amat = &(ml->Amat[current_level]);
   ML_Aggregate * ag = (ML_Aggregate *)data;
   
   if (current_level == 0 ) return(-1);

   total_size = Amat->invec_leng;
   ML_gsum_scalar_int(&total_size, &temp, ml->comm);
   if ( total_size <= ag->max_coarse_size ) return(-1);

   return(current_level-1);
}

/* ************************************************************************* */
/* function for enforcing a 2-level scheme                                   */
/* ------------------------------------------------------------------------- */

int ML_AGG_Increment_Two_Level(ML *ml,int current_level,
                               void * data)
{
  ML_Aggregate *ag = (ML_Aggregate *)data;
  (void) ml;
  if ( current_level == ag->begin_level ) return (current_level+1);
  return(-1);
}

/* ************************************************************************* */
/* function for enforcing a 2-level scheme                                   */
/* ------------------------------------------------------------------------- */

int ML_AGG_Decrement_Two_Level(ML *ml,int current_level,void * data)
{
  ML_Aggregate *ag = (ML_Aggregate *)data;
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
      printf("Invalid getrow id (level %d)\n",getrow_obj->ML_id);
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
            threshold = ML_dabs(values[i])*widget->drop_tol;
            break;
         }
      }
      j = 0;
      dropped = 0.0;
      for (i = 0; i < row_lengths[0]; i++) 
      {
         if ( ML_dabs(values[i]) >= threshold) 
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
#ifdef RST_MODIF
   if (diag_val == 0.) { row_lengths[0] = 0; return 1; }
   for (i = 0; i < row_lengths[0]; i++) 
      values[i] *= -widget->omega/diag_val;
   values[diag] += 1.;
#else
#ifndef MB_MODIF
   if (ML_dabs(diag_val) > 0.0)
   {
      for (i = 0; i < row_lengths[0]; i++) 
         values[i] *= (-widget->omega)/diag_val;
      values[diag] += 1.;
   }
#else
   for (i = 0; i < row_lengths[0]; i++) 
      values[i] *= -widget->omega;
   values[diag] += 1.;
#endif
#endif

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

int ML_AGG_Gen_DDProlongator(ML *ml,int level, int clevel, void *data)
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
   ML_Aggregate * ag = (ML_Aggregate *)data;
   
#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   /* ----------------------------------------------------------------- */
   /* coarsen local smoothed aggregation method                         */
   /* ----------------------------------------------------------------- */

   if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
      printf("Aggregation : building multilevel hierarchy at level %d\n",level);
   widget.near_bdry = NULL; 
   Amat     = &(ml->Amat[level]);
   Nfine    = Amat->outvec_leng;
   getrow_obj = Amat->getrow;
   if (getrow_obj->ML_id == ML_EXTERNAL) getrowfunc = getrow_obj->external;
   else                                  getrowfunc = getrow_obj->internal;
   max_nz_per_row = 100;
   col_ind = (int *)    ML_allocate( max_nz_per_row * sizeof(int) );
   col_val = (double *) ML_allocate( max_nz_per_row * sizeof(double) );
   nnz = 0;
   for ( i = 0; i < Nfine; i++ )
   {
      while (getrowfunc(Amat->data,1,&i,max_nz_per_row,col_ind,col_val,&k)== 0)
      {
         ML_free( col_ind );
         ML_free( col_val );
         max_nz_per_row = max_nz_per_row * 2 + 1;
         col_ind = (int *)    ML_allocate( max_nz_per_row * sizeof(int) );
         col_val = (double *) ML_allocate( max_nz_per_row * sizeof(double) );
      }
      nnz += k;
   }
   ML_free( col_ind );
   ML_free( col_val );
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

   if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
      printf("Aggregation : setting up diagonal block at level %d\n",level);

   newNlevels = 15;
   ML_Create(&newml, newNlevels);
   ML_Set_OutputLevel(newml, 0);
   ML_Set_ResidualOutputFrequency(newml, 0);
   ML_Set_Comm_MyRank(newml, 0);
   ML_Set_Comm_Nprocs(newml, 1);
   nbytes = sizeof(struct ML_AGG_Matrix_Context);
   context = (struct ML_AGG_Matrix_Context *) ML_allocate( nbytes );
   context->Amat = Amat;
   context->near_bdry = NULL;
   ML_Init_Amatrix(newml, newNlevels-1, Nfine,  Nfine, (void *) context);
   ML_Set_Amatrix_Matvec(newml,  newNlevels-1, ML_AGG_DD_Matvec);
   newml->Amat[newNlevels-1].data_destroy = ML_AGG_Matrix_Context_Clean;
   newml->Amat[newNlevels-1].N_nonzeros = 5 * Nfine;
   ML_Set_Amatrix_Getrow(newml, newNlevels-1, ML_AGG_DD_Getrow, NULL, Nfine);
   diagonal = (double *) ML_allocate(Nfine * sizeof(double));
   ML_AGG_Extract_Diag(Amat, diagonal);
   ML_Set_Amatrix_Diag( newml, newNlevels-1, Nfine, diagonal);
   ML_free( diagonal );
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

   if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
      printf("Aggregation : computing max eigenvalues at level %d\n",level);

/*
   if ( ag->smoothP_damping_factor != 0.0 )
*/
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
   /*
     }
     else
     {
     widget.omega = 0.0;
     ml->spectral_radius[level] = 1.0;
     }
      */

   /* ================================================================= */
   /* set up smoothed prolongator (I - alpha D^(-1) A) P                */
   /* 1. compute P (local ml P applied to 1)                            */
   /* 2. A * P                                                          */
   /* 3. D^{-1} * A * P                                                 */
   /* 4. P - alpha D^{-1} * A * P                                       */
   /* ================================================================= */

   i = 1;
   j = ML_gmax_int(i, ml->comm );
   if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
      printf("Aggregation : computing tentative prolongators at level %d\n",level);

   /* ----------------------------------------------------------------- */
   /* 1. compute tentP (local ml P applied to 1)                        */
   /* ----------------------------------------------------------------- */

   darray  = (double *) ML_allocate( Nfine * sizeof(double) );
   darray2 = (double *) ML_allocate( Nfine * sizeof(double) );

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
   ML_free( darray2 );
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
   p_cols = (int *) ML_allocate(sizeof(int));
   p_cols[0] = 0;
   p_aa = (double **) ML_allocate(sizeof(double*));
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
   /*
   if ( 1 )
   {
   */
      i = 1;
      j = ML_gmax_int(i, ml->comm );
      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("Aggregation : computing comm pattern of A*tentP at level %d\n",
              level);

      APMat = &(ml->Pmat[clevel]);
      ML_2matmult(Amat, tentP, APMat, ML_CSR_MATRIX );
      ML_AGG_Extract_Matrix(APMat, &ap_ncols, &ap_cols, &ap_aa);
 
      i = 1;
      j = ML_gmax_int(i, ml->comm );
      if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
         printf("Aggregation : computing prolongators at level %d\n",level);

      ML_Set_MaxIterations(newml, 10);
      darray  = (double *) ML_allocate( Nfine * sizeof(double) );
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
      ML_free( darray );
   /*
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
   */
   if ( p_cols != NULL ) ML_free( p_cols );
   for ( i = 0; i < p_ncols; i++ ) ML_free( p_aa[i] );
   if ( p_aa != NULL ) ML_free( p_aa );

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

   ML_free( ap_cols );
   for ( i = 0; i < ap_ncols; i++ ) ML_free( ap_aa[i] );
   ML_free( ap_aa );
   ML_Destroy(&newml);
   ML_Operator_Destroy(&tentP);

   i = 1;
   j = ML_gmax_int(i, ml->comm );
   if ( ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
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
   col_ind = (int *)    ML_allocate( max_row_nnz * sizeof(int) );
   col_val = (double *) ML_allocate( max_row_nnz * sizeof(double) );

   for ( i = 0; i < nRows; i++ )
   {
      while (getrowfunc(Amat->data,1,&i,max_row_nnz,col_ind,col_val,&m)== 0)
      {
         ML_free( col_ind );
         ML_free( col_val );
         max_row_nnz = max_row_nnz * 2 + 1;
         col_ind = (int *)    ML_allocate( max_row_nnz * sizeof(int) );
         col_val = (double *) ML_allocate( max_row_nnz * sizeof(double) );
      }
      dtmp = 0.0;
      
      for ( j = 0; j < m; j++ )
      {
         index = col_ind[j];
         if ( index < nRows ) dtmp += ( col_val[j] * p[index] );
      }
      ap[i] = dtmp;
   }
   ML_free( col_ind );
   ML_free( col_val );

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
      local_ind = (int *)    ML_allocate( alloc_space * sizeof(int));
      local_val = (double *) ML_allocate( alloc_space * sizeof(double));
   }
   status = getrowfunc(Amat->data, 1, rowlist, alloc_space, local_ind,
                       local_val, rowcnt);
   if ( status == 0 ) 
   {
      ML_free( local_ind );
      ML_free( local_val );
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
   ML_free( local_ind );
   ML_free( local_val );
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
   col_ind = (int *)    ML_allocate( max_row_nnz * sizeof(int) );
   col_val = (double *) ML_allocate( max_row_nnz * sizeof(double) );

   for ( i = 0; i < nRows; i++ )
   {
      while (getrowfunc(Amat->data,1,&i,max_row_nnz,col_ind,col_val,&m)== 0)
      {
         ML_free( col_ind );
         ML_free( col_val );
         max_row_nnz = max_row_nnz * 2 + 1;
         col_ind = (int *)    ML_allocate( max_row_nnz * sizeof(int) );
         col_val = (double *) ML_allocate( max_row_nnz * sizeof(double) );
      }
      for (j = 0; j < m; j++) if (col_ind[j] == i) diagonal[i] = col_val[j];
   }
   ML_free( col_ind );
   ML_free( col_val );

   return 1;
}

/* ************************************************************************* */
/* destroy aggregate matrix context                                          */
/* ------------------------------------------------------------------------- */

void ML_AGG_Matrix_Context_Clean(void *data)
{
   struct ML_AGG_Matrix_Context *context;

   context = (struct ML_AGG_Matrix_Context *) data;
   ML_free(context);
}

/* ************************************************************************* */
/* solve local subproblem using smoothed aggregation                         */
/* ------------------------------------------------------------------------- */

int ML_AGG_DD_Solve(void *data, int leng1, double *outvec, int leng2, 
                    double *invec)
{
   ML  *ml = (ML *) data;
   ML_Solve_MGV( ml, invec, outvec );
   ML_avoid_unused_param( (void *) &leng1);
   ML_avoid_unused_param( (void *) &leng2);
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
   col_ind = (int *)    ML_allocate( max_size * sizeof(int) );
   col_val = (double *) ML_allocate( max_size * sizeof(double) );
   nnz = 0;
   for ( i = 0; i < local_nrows; i++ )
   {
      while (getrowfunc(mat->data,1,&i,max_size,col_ind,col_val,&row_size)==0)
      {
         ML_free( col_ind );
         ML_free( col_val );
         max_size = max_size *2 + 1;
         col_ind = (int *)    ML_allocate( max_size * sizeof(int) );
         col_val = (double *) ML_allocate( max_size * sizeof(double) );
      }
      nnz += row_size;
      if ( row_size > max_size ) max_size = row_size;
   }
   ML_free( col_ind );
   ML_free( col_val );

   /* ----------------------------------------------------------------- */
   /* extract matrix                                                    */
   /* ----------------------------------------------------------------- */

   col_ind = (int *) ML_allocate( nnz * sizeof(int));
   col_val = (double *) ML_allocate( nnz * sizeof(double));
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
   local_cols = (int *) ML_allocate(local_ncols * sizeof(int));
   for ( i = 0; i < local_ncols; i++ ) local_cols[i] = col_ind[i];
   ML_free( col_ind );
   ML_free( col_val );

   /* ----------------------------------------------------------------- */
   /* fill in the matrix                                                */
   /* ----------------------------------------------------------------- */

   local_vals = (double **) ML_allocate(local_ncols * sizeof(double*));
   for ( i = 0; i < local_ncols; i++ )
   { 
      local_vals[i] = (double *) ML_allocate(local_nrows * sizeof(double));
      for ( j = 0; j < local_nrows; j++ ) local_vals[i][j] = 0.0;
   }

   col_ind = (int *)    ML_allocate( max_size * sizeof(int));
   col_val = (double *) ML_allocate( max_size * sizeof(double));
   for ( i = 0; i < local_nrows; i++ )
   {
      getrowfunc(mat->data,1,&i,max_size,col_ind,col_val,&row_size);
      for ( j = 0; j < row_size; j++ )
      {
         index = ML_sorted_search( col_ind[j], local_ncols, local_cols);
         if ( index >= 0 ) local_vals[index][i] = col_val[j];
      }
   }
   ML_free( col_ind );
   ML_free( col_val );

   (*ncols) = local_ncols;
   (*cols)  = local_cols;
   (*vals)  = local_vals;
   return 1;
}

/* ************************************************************************* */
/* generate smooth prolongator for 2-level DD method                         */
/* ------------------------------------------------------------------------- */

int ML_AGG_Gen_DDProlongator2(ML *ml,int level, int clevel, void *data)
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
   ML_Aggregate *ag = (ML_Aggregate *)data;
   
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
      context = (struct ML_AGG_Matrix_Context *) ML_allocate( nbytes );
      context->Amat = Amat;
      context->near_bdry = NULL;
      ML_Init_Amatrix(newml, newNlevels-1, Nfine,  Nfine, (void *) context);
      ML_Set_Amatrix_Matvec(newml,  newNlevels-1, ML_AGG_DD_Matvec);
      newml->Amat[newNlevels-1].data_destroy = ML_AGG_Matrix_Context_Clean;
      newml->Amat[newNlevels-1].N_nonzeros = 5 * Nfine;
      ML_Set_Amatrix_Getrow(newml, newNlevels-1, ML_AGG_DD_Getrow, NULL, Nfine);
      diagonal = (double *) ML_allocate(Nfine * sizeof(double));
      ML_AGG_Extract_Diag(Amat, diagonal);
      ML_Set_Amatrix_Diag( newml, newNlevels-1, Nfine, diagonal);
      ML_free( diagonal );
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

   ML_2matmult(AGGsmoother, tentP, &(ml->Pmat[clevel]), ML_CSR_MATRIX );

   ML_Operator_Destroy(&tentP);
   ML_Operator_Destroy(&AGGsmoother);
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
                    Amatrix->comm, ML_OVERWRITE,NULL);
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

   
   
   ML_free(rowi_col); ML_free(rowi_val);
   rowi_col = NULL; rowi_val = NULL;
   allocated = 0; 

   ML_free(dtemp);

   return 0;
}

/******************************************************************************
Regenerate the multigrid hierarchy with the existing restriction and
prolongation operators.
******************************************************************************/

int  ML_Gen_MGHierarchy_ReuseExistingOperators(ML *ml)
{
   int mesh_level, old_mesh_level;
   ML_Operator *mat;


   mesh_level = ml->ML_finest_level;

   while( ml->SingleLevel[mesh_level].Rmat->to != NULL) {
     old_mesh_level = mesh_level;
     mesh_level = ml->SingleLevel[mesh_level].Rmat->to->levelnum;
     mat = &(ml->Amat[mesh_level]);
     ML_Operator_Clean(mat);
     ML_Operator_Init(mat,ml->comm);
     ML_Gen_AmatrixRAP(ml, old_mesh_level, mesh_level);
   }

   return 0;
}
/******************************************************************************
Regenerate the multigrid hierarchy using smoothed aggregation reusing the 
existing aggregates.
******************************************************************************/

int  ML_Gen_MGHierarchy_UsingSmoothedAggr_ReuseExistingAgg(ML *ml,
							   ML_Aggregate *ag)
{
   int mesh_level, old_mesh_level;
   ML_Operator *mat;


   mesh_level = ml->ML_finest_level;
   if (ag->keep_P_tentative != ML_YES) {
     printf("ML_Gen_MGHierarchy_UsingSmoothedAggr_ReuseExistingAgg: must save\n");
     printf("   aggregation information by setting ML_Aggregate_Set_Reuse(...)\n");
     exit(-1);
   }

   while( ml->SingleLevel[mesh_level].Rmat->to != NULL) {
     old_mesh_level = mesh_level;
     mesh_level = ml->SingleLevel[mesh_level].Rmat->to->levelnum;
     /* clean and regenerate P */

     mat = &(ml->Pmat[mesh_level]);
     if (ag->smoothP_damping_factor != 0.0 ) {
       ML_Operator_Clean(mat);
       ML_Operator_Init(mat,ml->comm);
       ML_AGG_Gen_Prolongator(ml, old_mesh_level, mesh_level, (void*) ag);
     }

     /* clean and regenerate R */

     mat = &(ml->Rmat[old_mesh_level]);
     if (ag->smoothP_damping_factor != 0.0 ) {
       ML_Operator_Clean(mat);
       ML_Operator_Init(mat,ml->comm);
       ML_Gen_Restrictor_TransP(ml, old_mesh_level, mesh_level);
     }

     /* clean and regenerate A */

     mat = &(ml->Amat[mesh_level]);
     ML_Operator_Clean(mat);
     ML_Operator_Init(mat,ml->comm);
     ML_Gen_AmatrixRAP(ml, old_mesh_level, mesh_level);
   }

   return 0;
}

/* ************************************************************************* */
/* new version of ML_Gen_MGHierarchy_UsingAggregation                        */
/* I dropped out the domain_decomposition lines                              */
/* ------------------------------------------------------------------------- */

int ML_Gen_MultiLevelHierarchy_UsingAggregation(ML *ml, int start, 
						int increment_or_decrement,
						ML_Aggregate *ag)
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

   /* FIXME: I don't know when P_tentative is freed !!!!! */
   /*   if( ag->smoothP_damping_factor == 0.0 ) ag->Restriction_smoothagg_transpose == ML_FALSE;*/
   if( ag->Restriction_smoothagg_transpose == ML_TRUE ) ag->keep_P_tentative = ML_TRUE;
   
   idata = 0;
   idata = ML_gmax_int(idata, ml->comm);
   /*
   if ( ml->comm->ML_mypid == 0 && 9 < ML_Get_PrintLevel()) 
      ML_Aggregate_Print( ml_ag );
   */
#ifdef ML_TIMING
   t0 = GetClock();
#endif
   idata = ML_gmax_int(idata, ml->comm);

   if (increment_or_decrement == ML_INCREASING)
   {
     level = ML_Gen_MultiLevelHierarchy(ml, start,
					ML_AGG_Increment_Level,
					ML_MultiLevel_Gen_Restriction,
					ML_MultiLevel_Gen_Prolongator,
					(void *)ml_ag);
   }
   else if (increment_or_decrement == ML_DECREASING)
   {
     level = ML_Gen_MultiLevelHierarchy(ml, start,
					ML_AGG_Decrement_Level,
					ML_MultiLevel_Gen_Restriction,
					ML_MultiLevel_Gen_Prolongator,
					(void *)ml_ag);
   }
   else {
     if ( ml->comm->ML_mypid == 0 ) 
       {
	 printf("ML_Gen_MultiLevelHierarchy_UsingAggregation : Unknown ");
         printf("increment_or_decrement choice\n");
       }
     exit(1);
   }
#ifdef ML_TIMING
   t0 = GetClock() - t0;
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag < ML_Get_PrintLevel()) 
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
   if ( ml->comm->ML_mypid == 0 && ml_ag->print_flag < ML_Get_PrintLevel()) 
      ML_Aggregate_Print_Complexity( ml_ag );
   idata = ML_gmax_int(idata, ml->comm);

   if ( ag == NULL ) ML_Aggregate_Destroy( &ml_ag );
   return(level);
}

/* ************************************************************************* */
/* generate multilevel hierarchy given a subroutine for generating           */
/* prolongation operators. New version of ML_Gen_MGHierarchy.                */
/* ------------------------------------------------------------------------- */

int ML_Gen_MultiLevelHierarchy(ML *ml, int fine_level,
        int (*user_next_level)(ML *, int, void *),
        int (*user_gen_restriction)(ML *, int, int, void *),
        int (*user_gen_prolongator)(ML *, int, int, void *),
        void *user_data)
{
   int level, next, flag, count=1;
   int i, j, bail_flag, N_input_vector;
   ML_Operator *Pmat;
   ML_CommInfoOP *getrow_comm;
   
#ifdef ML_TIMING
   double t0;
#endif

   ml->ML_finest_level = fine_level;
   level = fine_level;
   next  = user_next_level(ml, level, user_data);

   while (next >= 0) 
   {
      if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel()) 
         printf("ML_Gen_MultiLevelHierarchy (level %d) : Gen Restriction and Prolongator \n",
		level );

      flag = (*user_gen_prolongator)(ml, level, next, user_data);
      if (flag < 0) break;

      /* Now check to make sure prolongator has zero columns. */
      Pmat = ml->Pmat+next;
      bail_flag = 0;
      N_input_vector = Pmat->invec_leng;
      getrow_comm = Pmat->getrow->pre_comm;
      if ( getrow_comm != NULL)
      {
         for (i = 0; i < getrow_comm->N_neighbors; i++) {
            for (j = 0; j < getrow_comm->neighbors[i].N_send; j++) {
               if (getrow_comm->neighbors[i].send_list[j] >= N_input_vector) {
                  bail_flag = 1;
               }
            }
         }
      }
      /* If check has failed on any processor, clean up current level & break
         from main loop. */
      ML_gsum_scalar_int(&bail_flag,&j,ml->comm);
      if (bail_flag)
      {
         if (Pmat->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
            printf("(%d) In ML_Gen_MultiLevelHierarchy: Bailing from AMG hierarchy build on level %d, where fine level = %d ........\n",
                   Pmat->comm->ML_mypid,level,fine_level);
            fflush(stdout);
         }
         if (ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
            printf("(%d) In ML_Gen_MultiLevelHierarchy: Nlevels = %d fine_level = %d  coarsest_level = %d\n",
               ml->comm->ML_mypid,fine_level-count+1,fine_level,count);
            fflush(stdout);
         }
         break; /* from main loop */
 
      }
      /* end of check */

      (*user_gen_restriction)(ml, level, next,user_data);

#ifdef ML_TIMING
      t0 = GetClock();
#endif

      if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel())
	printf("ML_Gen_MultiLevelHierarchy (level %d) : Gen RAP\n", level);
      ML_Gen_AmatrixRAP(ml, level, next);

#ifdef ML_TIMING
      t0 = GetClock() - t0;
      if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel()) 
         printf("ML_Gen_MultiLevelHierarchy (level %d) : RAP time = %e\n", level, t0);
#endif

      level = next;
      next  = user_next_level(ml, next, user_data);

      count++;
   }
   return(count);
}

/* ************************************************************************* */
/* generate smooth prolongator                                               */
/* This is supposed to create a `real' prolongator (that is, not to define R */
/* in the case of R = (I - w A ) Pt^T)                                       */
/* NOTE: data is suppposed to be a pointer to the ML_Aggregate structure.    */
/*                                                                           */
/* If ag->smoothP_damping_factor != 0.0, the function estimate lambda max.   */
/* ------------------------------------------------------------------------- */

int ML_MultiLevel_Gen_Prolongator(ML *ml,int level, int clevel, void *data)
{
   int         Nfine, gNfine;
   double      max_eigen = -1.;
   ML_Operator *Amat;
   ML_Operator **prev_P_tentatives;
   struct      ML_AGG_Matrix_Context widget;
   ML_Aggregate *ag = (ML_Aggregate *) data;
   struct ML_Field_Of_Values * fov;
   double dtemp, dtemp2, eta;
#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   widget.near_bdry = NULL;
   Amat     = &(ml->Amat[level]);
   Amat->num_PDEs = ag->num_PDE_eqns;
   prev_P_tentatives = ag->P_tentative;

   Nfine    = Amat->outvec_leng;
   gNfine   = ML_Comm_GsumInt( ml->comm, Nfine);
   ML_Aggregate_Set_CurrentLevel( ag, level );

   /* ********************************************************************** */
   /* May require computationof field-of-values for non-diagonally scaled A  */
   /* ********************************************************************** */

   if( ag->field_of_values != NULL ) {
     
     fov = (struct ML_Field_Of_Values * )(ag->field_of_values);

     if( fov->compute_field_of_values_non_scaled == ML_YES ) {

       ML_Anasazi_Get_FieldOfValuesBoxNonScaled_Interface(Amat,fov);
       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nNon-Scaled Field of Values Box (level %d) : Max Real = %e\n",
		level,
		fov->real_max );
	 printf("Non-Scaled Field of Values Box (level %d) : Max Imag = %e\n",
		level,
		fov->imag_max );
	 printf("Non-Scaled Field of Values Box (level %d) : eta = %e\n\n",
		level,
		fov->eta );
       } 
     } 
   }

       
   /* ********************************************************************** */
   /* May require field-of-values computations for classic ML                */
   /* ********************************************************************** */

   if( ag->Restriction_smoothagg_transpose == ML_FALSE &&
       ag->field_of_values != NULL ) {
     
     fov = (struct ML_Field_Of_Values * )(ag->field_of_values);

     if( fov->compute_field_of_values == ML_YES ) {

       ML_Anasazi_Get_FieldOfValuesBox_Interface(Amat,fov);

       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nField of Values Box (level %d) : Max Real = %e\n",
		level,
		fov->real_max );
	 printf("Field of Values Box (level %d) : Max Imag = %e\n",
		level,
		fov->imag_max );
	 printf("Field of Values Box (level %d) : eta = %e\n\n",
		level,
		fov->eta );
       }
       
     }
     
   }
   
   /* ********************************************************************** */
   /* Methods based on field-of-values requires to stick some parameres now  */
   /* This is not an error! Here use R (this P will become R later)          */
   /* ********************************************************************** */
   
   if( ag->Restriction_smoothagg_transpose == ML_TRUE ) {
     
     fov = (struct ML_Field_Of_Values * )(ag->field_of_values);

     /* compute box surrounding field-of-values */

#if defined(HAVE_ML_ANASAZI) && defined(HAVE_ML_TEUCHOS)

     if( fov->compute_field_of_values == ML_YES && fov->choice != 1 ) {
       
       ML_Anasazi_Get_FieldOfValuesBox_Interface(Amat,fov);
       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nField of Values Box (level %d) : Max Real = %e\n",
		level,
		fov->real_max );
	 printf("Field of Values Box (level %d) : Max Imag = %e\n",
		level,
		fov->imag_max );
	 printf("Field of Values Box (level %d) : eta = %e\n\n",
		level,
		fov->eta );
       }
       
     }
     
     if( fov->choice == 0 ) {

       fov->eta = 0;
       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\n(level %d) : Using non-smoothed aggregation\n\n",
		level );
       }
       
       ag->smoothP_damping_factor = 0.00000001;
       Amat->lambda_max = 1.0;

     } else if( fov->choice == 1 ) {
       
       ML_Anasazi_Get_FieldOfValuesBox_Interface(Amat,fov);
       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nField of Values Box (level %d) : Max Real = %e\n",
		level,
		fov->real_max );
	 printf("Field of Values Box (level %d) : Max Imag = %e\n",
		level,
		fov->imag_max );
	 printf("Field of Values Box (level %d) : eta = %e\n\n",
		level,
		fov->eta );
       }

       eta = fov->eta;
       dtemp = fov->R_coeff[0] + eta * (fov->R_coeff[1]) + pow(eta,2) * (fov->R_coeff[2]);
       if( dtemp < 0.0 ) dtemp = 0.000001;
       dtemp2 = fov->real_max;
       
       ag->smoothP_damping_factor = dtemp;
       Amat->lambda_max = dtemp2;
       
     } else if( fov->choice == 2 ) {

       ML_Anasazi_Get_FieldOfValuesBox_Interface(Amat,fov);
       fov->eta = sqrt(pow(fov->real_max,2) + pow(fov->imag_max,2));

       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nLargest eigenvalue (in modulus) = %e\n\n",
		fov->eta );
       }

       eta = fov->eta;
       dtemp = fov->R_coeff[0] + eta * (fov->R_coeff[1]) + pow(eta,2) * (fov->R_coeff[2]);
       if( dtemp < 0.0 ) dtemp = 0.000001;
       dtemp2 = fov->real_max;
       
       ag->smoothP_damping_factor = dtemp;
       Amat->lambda_max = dtemp2;
       
     } else {
       
       fprintf( stderr,
		"ERROR: value of choice not correct (%d)\n"
		"ERROR: (file %s, line %d)\n",
		fov->choice,
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
       
     }
#else
     fprintf( stderr,
	      "ERROR: You must compile with options --enable-anasazi\n"
	      "ERROR: and --enable-teuchos for eigen-analysis\n"
	      "ERROR: (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
#endif

     if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
       printf("Restriction smoother (level %d) : damping factor = %e\n"
	      "Restriction smoother (level %d) : ( = %e / %e)\n\n",
	      level,
	      ag->smoothP_damping_factor/Amat->lambda_max,
	      level,
	      ag->smoothP_damping_factor,
	      Amat->lambda_max);
     }

     ml->symmetrize_matrix = ML_FALSE;
     ag->keep_P_tentative = ML_YES;
     ag->use_transpose = ML_TRUE;

     max_eigen = Amat->lambda_max;
     
   } else {

     /* This is the normal ML way */
     
     max_eigen = Amat->lambda_max;

   }
   
   ML_AGG_Gen_Prolongator(ml,level,clevel,data);   
   
   return 0;
   
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level        */
/* number larger than the fine levels                                        */
/* ------------------------------------------------------------------------- */

int ML_MultiLevel_Gen_Restriction(ML *ml,int level, int next, void *data)
{

  ML_Operator *Amat;
  ML_Operator **prev_P_tentatives;
  ML_Aggregate *ag = (ML_Aggregate *) data;

  struct ML_Field_Of_Values * fov;
  double dtemp, dtemp2, eta;
  char str[80];

  prev_P_tentatives = ag->P_tentative;
  
  Amat = &(ml->Amat[level]);

  if( ag->Restriction_smoothagg_transpose == ML_TRUE ) {

    if( ag->use_transpose != ML_TRUE ) {
      fprintf( stderr,
	       "ERROR: Something went **very** wrong in `ML_MultiLevel_Gen_ProlongatorForRestriction'\n"
	       "ERROR: (file %s, line %d)\n",
	       __FILE__,
	       __LINE__ );
    
      exit( EXIT_FAILURE );
    }  
    
    /* ********************************************************************** */
    /* Previous P has been built based on A^T and not on A. Now, first we     */
    /* transpose P into R (formed with A^T, so that now R is based on A).     */
    /* NOTE: the damping parameter previously used in P is the one for R.     */
    /* Then, we will have to build a new P based on A.                        */
    /* ********************************************************************** */

    ML_Gen_Restrictor_TransP(ml, level, next);

    /* ********************************************************************** */
    /* Now rebuilt P using the saved tentative guy. We need to compute the    */
    /* good parameter for damping.                                            */
    /* We erase the old Pmat too.                                             */
    /* ********************************************************************** */

    ag->use_transpose = ML_FALSE;
    ML_Operator_Clean( &(ml->Pmat[next]) );
    ML_Operator_Init( &(ml->Pmat[next]), ml->comm );
    ML_Operator_Set_1Levels(&(ml->Pmat[next]), &(ml->SingleLevel[next]), NULL);
    ML_Operator_Set_BdryPts(&(ml->Pmat[next]), NULL);
    sprintf(str,"Pmat_%d",next); ML_Operator_Set_Label( &(ml->Pmat[next]),str);

    fov = (struct ML_Field_Of_Values *)(ag->field_of_values);
    eta = fov->eta;
    
    dtemp = fov->P_coeff[0] + eta * (fov->P_coeff[1]) + pow(eta,2) * (fov->P_coeff[2]);
    if( dtemp < 0.0 ) dtemp = 0.000001;
    dtemp2 = fov->real_max;
 	      
    ag->smoothP_damping_factor = dtemp;
    Amat->lambda_max = dtemp2;

    if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
      printf("Prolongator smoother (level %d) : damping parameter = %e\n"
	     "Prolongator smoother (level %d) : ( = %e / %e )\n\n",
	     level,
	     ag->smoothP_damping_factor/dtemp2,
	     level,
	     ag->smoothP_damping_factor,
	     dtemp2 );
    }    
    
    /* use old-fashioned functions to create the actual prolongator based on A */
    
    ML_AGG_Gen_Prolongator(ml,level,next,data);

  }  else {
    
    ML_Gen_Restrictor_TransP(ml, level, next);
    
  }

  return 0;
  
}
/*
 * This routine is used for block diagonal scaling in the prolongator
 * Smoother. Essentially, it takes the matrix Ptemp and multiplies
 * it by the block diagonal scaling which is inside the MLSthing. 
 *
 * The tricky thing is that we want to access Ptemp by Rows ... so
 * what we do is we transpose it. However, we have to be sure that
 * we use the funny transpose that leaves all the same data within
 * processor (and has the funny post communication). After we are
 * finished doing the computation (which is all local) we need
 * to transpose things back. The only catch is that this 
 * tranpose routine does not work with matrices that have 
 * post communication (which we now have due to the first call
 * to the transpose function). So we take the transposed data
 * and put it into the original Ptemp (that has the correct
 * communication data structure).
 */

int ML_AGG_DinvP(ML_Operator *Ptemp, struct MLSthing *mls_widget,
		 int blk_size)
{
  ML_Operator *Rtemp, *P2temp;
  struct ML_CSR_MSRdata *csr_data;
  int Num_blocks = 0, i, j, block, column;
  int *rowptr, *columns;
  int *new_columns, first, last, nz_ptr;
  double *new_values, *values;
  int **perms;
  double **blockdata;
  char N[2];
  unsigned int   itmp=0;
  int info, one = 1;


  Rtemp = ML_Operator_Create(Ptemp->comm);
  ML_Operator_Transpose(Ptemp, Rtemp);

  csr_data = (struct ML_CSR_MSRdata *) Rtemp->data;

  /* count the total number of blocks */

  rowptr  = csr_data->rowptr;
  columns = csr_data->columns;
  values  = csr_data->values;
  for (i = 0; i < Rtemp->outvec_leng; i++) {
    ML_az_sort( &(columns[rowptr[i]]), 
		rowptr[i+1]-rowptr[i+1], NULL, 
		&(values[rowptr[i]])); 
    block = -1;
    for (j = rowptr[i]; j < rowptr[i+1]; j++) {
      if (columns[j]/blk_size != block) {
	Num_blocks++;
	block = columns[j]/blk_size;
      }
    }
  }
  /* Make sure that all columns within a block */
  /* are present (i.e. pad with zeros if some  */
  /* are missing).                             */
  if (blk_size*Num_blocks != rowptr[Rtemp->outvec_leng]) {
    new_columns= (int    *) ML_allocate(sizeof(int)*Num_blocks*blk_size);
    new_values = (double *) ML_allocate(sizeof(double)*Num_blocks*blk_size);
    first = rowptr[0];
       nz_ptr = 0;
       for (i = 0; i < Rtemp->outvec_leng; i++) {
	 j = first;
	 last = rowptr[i+1];

	 while ( j < last ) {
	   block = columns[j]/blk_size;
	   column = block*blk_size;
	   while (column < block*blk_size+blk_size) {
	     new_columns[nz_ptr] = column;
	     if (  (j < rowptr[i+1]) && (columns[j] == column)) {
	       new_values[nz_ptr++] = values[j];
	       j++;
	     }
	     else new_values[nz_ptr++] = 0.;
	     column++;
	   }
	 }
	 first = last;
	 rowptr[i+1] = nz_ptr;
       }
       ML_free(values); ML_free(columns);
       csr_data->values = new_values;
       csr_data->columns = new_columns;
       columns = csr_data->columns;
       values  = csr_data->values;

     }
     blockdata  = mls_widget->block_scaling->blockfacts;
     perms      = mls_widget->block_scaling->perms;
     strcpy(N,"N");

     for (i = 0; i < Rtemp->outvec_leng; i++) {
       for (j = csr_data->rowptr[i];j < csr_data->rowptr[i+1]; j += blk_size) {
	 block  = csr_data->columns[j]/blk_size;
	 MLFORTRAN(dgetrs)(N,&blk_size,&one,blockdata[block],&blk_size,
			   perms[block], &(values[j]),
			   &blk_size, &info, itmp);
	 if ( info != 0 ) {
	   printf("dgetrs returns with %d at block %d\n",info,i); 
	   exit(1);
	 }
       }
     }

     /* We need to transpose the matrix back. Unfortunately, */
     /* the routine ML_Operator_Transpose() takes a matrix   */
     /* with pre-communication and makes a matrix with post- */
     /* communication. This means that Rtemp has post-       */
     /* communication ... so the second call to ML_Operator_ */
     /* Transpose() will not get the communication right.    */
     /* So what we will do is swap the data pointers with    */
     /* the original Ptemp (which has the correct            */
     /* communication pointer).                              */
	
     P2temp = ML_Operator_Create(Ptemp->comm);
     ML_Operator_Transpose(Rtemp, P2temp);
     ML_Operator_Destroy(&Rtemp);

     csr_data = (struct ML_CSR_MSRdata *) P2temp->data;
     P2temp->data = Ptemp->data;
     Ptemp->data = csr_data;
     ML_Operator_Destroy(&P2temp);

     return 0;
}

