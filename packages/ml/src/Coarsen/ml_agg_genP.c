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
#include <stdlib.h>
#include "ml_struct.h"
#include "ml_smoother.h"
#include "ml_op_utils.h"
#include "ml_agg_genP.h"
#include "ml_memory.h"
#include "ml_lapack.h"
#include "ml_viz_stats.h"
#include "ml_agg_min_energy.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)
extern int ML_Anasazi_Get_FieldOfValuesBox_Interface(ML_Operator * Amat,
						     struct ML_Field_Of_Values * fov );
extern int ML_Anasazi_Get_FieldOfValuesBoxNonScaled_Interface(ML_Operator * Amat,
						     struct ML_Field_Of_Values * fov );
extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
                                               ML_Smoother* Smoother,
					       int MaxIters, double Tolerance,
					       int IsProblemSymmetric,
					       int UseDiagonalScaling,
					       double * LambdaMax );

#endif
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

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
   int    level, idata, i;
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
				    ML_AGG_Gen_DDProlongator, NULL, ml_ag);
      }
      else
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Increment_Level,
				    ML_AGG_Gen_Prolongator, NULL, ml_ag);
      }
   }
   else if (increment_or_decrement == ML_DECREASING)
   {
      if ( ml_ag->coarsen_scheme == 4 )
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Decrement_Two_Level,
				    ML_AGG_Gen_DDProlongator, NULL, ml_ag);
      }
      else
      {
         level = ML_Gen_MGHierarchy(ml, start, ML_AGG_Decrement_Level,
				    ML_AGG_Gen_Prolongator, NULL, ml_ag);
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
   if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 9 )
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

   if (ML_Get_PrintLevel() > 10) {
     for (i=0; i<level; i++) {
       int thisLevel = ml->LevelID[i];
         ML_Operator_Profile(ml->Amat+thisLevel,NULL);
         if (i != level-1) ML_Operator_Profile(ml->Rmat+thisLevel,NULL);
         if (i != 0)       ML_Operator_Profile(ml->Pmat+thisLevel,NULL);
     }
   }

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
		       void *data, ML_Aggregate *ag)
{
   int level, next, flag, count=1;
   int i, j, k, bail_flag, N_input_vector;
   ML_Operator *Pmat;
   ML_CommInfoOP *getrow_comm;
   ML_Operator *Ptent;
   ML_Aggregate_Viz_Stats *grid_info;
#ifdef ML_TIMING
   double t0;
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
   k = 0;
   ml->LevelID[k++] = fine_level;

   while (next >= 0)
   {

      /* This if-else supports an ALEGRA capability. */
      if (data == NULL)
        flag = (*user_gen_prolongator)(ml, level, next,(void *)ag);
      else
        flag = (*user_gen_prolongator)(ml, level, next,(void *)data);

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
      ML_Gen_Restrictor_TransP(ml, level, next, NULL);
      ML_Operator_ChangeToSinglePrecision(&(ml->Pmat[next]));
      ML_memory_check("L%d: TransP end",level);

#ifdef ML_TIMING
      t0 = GetClock();
#endif
      ML_Gen_AmatrixRAP(ml, level, next);
      ML_memory_check("L%d: RAP end",level);

      /* project the coordinates (if any) to the next coarser level */
      if (ag->P_tentative != NULL)
        Ptent = ag->P_tentative[next];
      else
        Ptent = &(ml->Pmat[next]);

      grid_info =(ML_Aggregate_Viz_Stats *) ml->Amat[level].to->Grid->Grid;
      if (grid_info) {
        if (grid_info->x != NULL) {
          if (ML_Get_PrintLevel() > 4 && ml->comm->ML_mypid == 0)
            printf("ML_Gen_MGHierarchy: Projecting node coordinates from level %d to level %d\n",
                   level,next);
          ML_Project_Coordinates(ml->Amat+level, Ptent, ml->Amat+next);
        }
      }

      ML_Comm_Barrier(ml->comm);

      ML_repartition_Acoarse(ml, level, next, ag, ML_TRUE, ML_FALSE);

      ML_Comm_Barrier(ml->comm);

      ML_Operator_ImplicitTranspose(&(ml->Rmat[level]),
      			    &(ml->Pmat[next]), ML_TRUE);

      ML_memory_check("L%d: repartition end",level);

#ifdef ML_TIMING
      t0 = GetClock() - t0;
      if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 10 )
         printf("RAP time for level %2d = %e\n", level, t0);
#endif

      level = next;
      ml->LevelID[k++] = next;
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
   int         Ncoarse, Nfine, gNfine, gNcoarse, ii, jj;
   double      max_eigen = -1.;
   ML_Operator *Amat, *Pmatrix = NULL, *AGGsmoother = NULL;
   ML_Operator **prev_P_tentatives;
   struct      ML_AGG_Matrix_Context widget;
   ML_Krylov   *kdata;
   ML_Operator *t2 = NULL, *t3 = NULL;
   ML_Aggregate * ag = (ML_Aggregate *) data;
   struct MLSthing *mls_widget = NULL;
   ML_Operator *blockMat = NULL, *Ptemp;
   int numSmSweeps;   /* polynomial degree of prolongator smoother */
   double *dampingFactors; /* coefficients of prolongator smoother */
   ML_Operator *tmpmat1=NULL,*tmpmat2=NULL;

   double t1=0, createPtentTime=0, smooPTime=0, eigenTime=0;
#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif
   widget.Adiag = NULL;

   StartTimer(&t1);

   Amat = &(ml->Amat[level]);
   numSmSweeps = ML_Aggregate_Get_DampingSweeps(ag,level);

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
#ifdef USE_MOREACCURATE
if (Amat->diagonal == NULL) ML_Operator_Getrow_Diag(Amat, &(widget.Adiag));
else ML_DVector_GetDataPtr( Amat->diagonal, &(widget.Adiag) );
#endif


   /*
   widget.near_bdry = (char *) ML_allocate(sizeof(char)*Amat->outvec_leng);
   ML_AGG_Compute_Near_Bdry(Amat, widget.near_bdry);
   */

   Nfine    = Amat->outvec_leng;
   gNfine   = ML_Comm_GsumInt( ml->comm, Nfine);
   ML_Aggregate_Set_CurrentLevel( ag, level );

   if ( (ag->smoothP_damping_factor!=0.0) && (numSmSweeps>0) )
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

       if (Pmatrix->getrow != NULL) {
         if (Pmatrix->getrow->pre_comm != NULL)
           ML_CommInfoOP_Destroy(&(Pmatrix->getrow->pre_comm));
       }
       ML_memory_free( (void**)&(Pmatrix->getrow) );
       Pmatrix->getrow = ml->Pmat[clevel].getrow;
       ml->Pmat[clevel].getrow = NULL;
       if (Pmatrix->label != NULL) ML_free(Pmatrix->label);
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
/* the next two lines added to prevent breakdown in adaptive SA */
   gNcoarse = gNcoarse / Pmatrix->num_PDEs;
   gNfine = gNfine / Amat->num_PDEs;
   if ( gNcoarse == 0 || ((1.0*gNfine)/(1.0*gNcoarse+0.1) < 1.05) )
   {
     if (( Pmatrix != NULL ) && (ag->smoothP_damping_factor != 0.0)
                             && (numSmSweeps>0) )
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

   StopTimer(&t1,&createPtentTime);
   StartTimer(&t1);

   if ( ag->smoothP_damping_factor != 0.0 && numSmSweeps > 0 )
   {
     /*********************************************************
     * If we symmetrize we need the symmetrized matrix so we
     * don't call ML_Gimmie_Eigenvalues() here.
     * ML_Gimmie_Eigenvalues(Amat, ML_DIAGSCALE,
     *             Amat->spectral_radius_scheme, ml->symmetrize_matrix);
     * max_eigen = Amat->lambda_max;
     *********************************************************/

     if (ml->symmetrize_matrix == ML_TRUE) {
       t2 = ML_Operator_Create(Amat->comm);
       ML_Operator_Transpose_byrow(Amat,t2);
       t3 = ML_Operator_Create(Amat->comm);
       ML_Operator_Add(Amat,t2,t3,ML_CSR_MATRIX,1.);
       max_eigen = t3->lambda_max;
     }

     if ((max_eigen < -666.) && (max_eigen > -667)) {
       if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5)
         printf("Calculating eigenvalue estimate using ");

       switch( Amat->spectral_radius_scheme ) {

       case ML_USE_CG:  /* compute it using CG */
         if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5)
           printf("CG method\n");

         kdata = ML_Krylov_Create( ml->comm );
         ML_Krylov_Set_PrintFreq( kdata, 0 );
         ML_Krylov_Set_MaxIterations(kdata, Amat->spectral_radius_max_iters);
         ML_Krylov_Set_ComputeEigenvalues( kdata );
         if (ml->symmetrize_matrix ==ML_TRUE) ML_Krylov_Set_Amatrix(kdata, t3);
         else ML_Krylov_Set_Amatrix(kdata, Amat);
         ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
         /* This is a bit screwy in that max_eigen corresponds to Dinv A. */
         /* This is good for Cheby smoothers but bad for smoothed agg     */
         /* if we actually filter A (significantly) as we should compute  */
         /* the eigenvalue of Dinv Afilt.                                 */
         max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);

         Amat->lambda_max = max_eigen;
         Amat->lambda_min = kdata->ML_eigen_min;
         ML_Krylov_Destroy( &kdata );
         if ( max_eigen <= 0.0 ) {
            printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
            max_eigen = 1.0;
         }

         break;

       case ML_USE_ANASAZI: /* Use Anasazi */
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)
         ML_Anasazi_Get_SpectralNorm_Anasazi(Amat, 0, 10, 1e-5,
                         ML_FALSE, ML_TRUE, &max_eigen);
         if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5)
           printf("Anasazi\n");
#else
         fprintf(stderr,
             "--enable-epetra --enable-anasazi --enable-teuchos required\n"
             "(file %s, line %d)\n",
             __FILE__,
             __LINE__);
         exit(EXIT_FAILURE);
#endif
         Amat->lambda_max = max_eigen;
         Amat->lambda_min = -12345.6789;
         if ( max_eigen <= 0.0 ) {
           printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
           max_eigen = 1.0;
         }

         break;

       case ML_USE_POWER: /* use ML's power method */
         if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5)
           printf("power method\n");
         kdata = ML_Krylov_Create( ml->comm );
         ML_Krylov_Set_PrintFreq( kdata, 0 );
         ML_Krylov_Set_MaxIterations(kdata, Amat->spectral_radius_max_iters);
         ML_Krylov_Set_ComputeNonSymEigenvalues( kdata );
         ML_Krylov_Set_Amatrix(kdata, Amat);
         ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
         max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);
         Amat->lambda_max = max_eigen;
         Amat->lambda_min = kdata->ML_eigen_min;
         ML_Krylov_Destroy( &kdata );
         if ( max_eigen <= 0.0 ) {
           printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
           max_eigen = 1.0;
         }

         break;

       default: /* using matrix max norm */
         if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5)
           printf("matrix max norm\n");
         max_eigen = ML_Operator_MaxNorm(Amat, ML_TRUE);
         break;

       } /* switch( Amat->spectral_radius_scheme ) */

     } else {
       if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5)
         printf("Using stashed eigenvalue estimate");
     } /* if ((max_eigen < -666.) && (max_eigen > -667)) */


     widget.omega  = ag->smoothP_damping_factor / max_eigen;
     ml->spectral_radius[level] = max_eigen;

     if ( ml->comm->ML_mypid == 0 && 7 < ML_Get_PrintLevel())
       printf("Gen_Prolongator (level %d) : Max eigenvalue = %2.4e\n",
          ag->cur_level, max_eigen);

   }
   else  /* damping fact = 0 ==> no need to compute spectral radius */
   {
      ml->spectral_radius[level] = 1.0;
      widget.omega  = 0.0;
      if ( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() )
        printf("\nProlongator/Restriction smoother (level %d) : damping = %2.3e , sweeps = %d\n", level, ag->smoothP_damping_factor, numSmSweeps );

   } /* if ( ag->smoothP_damping_factor != 0.0 ) */

   StopTimer(&t1,&eigenTime);
   StartTimer(&t1);

   /* Smooth tentative prolongator. */
   if ( ag->smoothP_damping_factor != 0.0 && numSmSweeps > 0 )
   {
     dampingFactors = (double *) ML_allocate( sizeof(double) * numSmSweeps );
     if (numSmSweeps == 1)
       dampingFactors[0] = ag->smoothP_damping_factor;
     else
       /* Calculate the proper Chebyshev polynomial coefficients. */
       ML_AGG_Calculate_Smoothing_Factors(numSmSweeps, dampingFactors);

     if ( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
       printf("\n");
       for (ii=0; ii<numSmSweeps; ii++)
       {
         if (ag->minimizing_energy == -1)
         {
           printf("\nProlongator smoother (level %d) : damping factor #%d = %2.4e\nProlongator smoother (level %d) : ( = %2.4e / %2.4e)\nNon-smoothed restriction is used.\n",
                  level, ii+1, dampingFactors[ii]/ max_eigen, level,
                  dampingFactors[ii], max_eigen );
         }
         else
         {
           printf("Prolongator/Restriction smoother (level %d) : damping factor #%d = %2.4e\nProlongator/Restriction smoother (level %d) : ( = %2.4e / %2.4e)\n",
                  level, ii+1, dampingFactors[ii]/ max_eigen, level,
                  dampingFactors[ii], max_eigen );
         }
         if (ml->symmetrize_matrix == ML_TRUE) printf("Smoothing with A + A^T\n");
         printf("\n");
       }
     }

     /* Create the prolongator smoother operator, I-omega*inv(D)*A. */
     AGGsmoother = ML_Operator_Create(ml->comm);
     ML_Operator_Set_Label(AGGsmoother,"Prolongator smoother");
     widget.drop_tol = ag->drop_tol_for_smoothing;
     if (ml->symmetrize_matrix == ML_TRUE) widget.Amat   = t3;
     else widget.Amat   = &(ml->Amat[level]);
     widget.aggr_info = ag->aggr_info[level];
     ML_Operator_Set_ApplyFuncData(AGGsmoother, widget.Amat->invec_leng,
                                   widget.Amat->outvec_leng, &widget,
                                   widget.Amat->matvec->Nrows, NULL, 0);
     ML_Operator_Set_Getrow(AGGsmoother,
                            widget.Amat->getrow->Nrows,
#ifdef USE_MOREACCURATE
                            ML_AGG_JacobiMoreAccurate_Getrows);
#else
                            ML_AGG_JacobiSmoother_Getrows);
#endif
     ML_CommInfoOP_Clone(&(AGGsmoother->getrow->pre_comm),
                           widget.Amat->getrow->pre_comm);

     tmpmat2 = Pmatrix;
     for (ii=0; ii < numSmSweeps; ii++)
     {
       /* Set the appropriate prolongator smoother damping factor. */
       widget.omega  = dampingFactors[ii] / max_eigen;

       if (ii < numSmSweeps-1) {
         tmpmat1 = tmpmat2;
         tmpmat2 = ML_Operator_Create(Amat->comm);
       }
       else {
         tmpmat1 = tmpmat2;
         tmpmat2 = &(ml->Pmat[clevel]);
       }

       if (ag->block_scaled_SA == 1) {
         /* Computed the following:
          *    a) turn off the usual 2 mat mult.
          *    b) Ptemp = A*P
          *    c) Ptemp = Dinv*Ptemp;
          *    d) do an ML_Operator_Add() with the original P.
          */
         Ptemp = ML_Operator_Create(Amat->comm);
         ML_2matmult(Amat, tmpmat1, Ptemp, ML_CSR_MATRIX );
         ML_AGG_DinvP(Ptemp, mls_widget, Amat->num_PDEs, ROW_SCALE_WITH_D);

         ML_Operator_Add(tmpmat1, Ptemp, tmpmat2, ML_CSR_MATRIX,
                         -dampingFactors[ii] / max_eigen);
         ML_Operator_Destroy(&Ptemp);
       }
       else
         ML_2matmult(AGGsmoother, tmpmat1, tmpmat2, ML_CSR_MATRIX );

       /*Free intermediate matrix I-omega*inv(D)*A, as long as it's not Ptent.*/
       if (tmpmat1 && tmpmat1 != Pmatrix) ML_Operator_Destroy(&tmpmat1);

     } /* for (ii=0; ii < numSmSweeps; ii++) */

     ML_Operator_Destroy(&AGGsmoother);
     ML_free(dampingFactors);

     if ((ml->symmetrize_matrix == ML_TRUE) || (ag->use_transpose == ML_TRUE) )
     {
       if (t3 != NULL) ML_Operator_Destroy(&t3);
       if (t2 != NULL) ML_Operator_Destroy(&t2);
     }
     if (ag->keep_P_tentative == ML_NO)
       ML_Operator_Destroy(&Pmatrix);
     else {
       if (prev_P_tentatives == NULL) {
         ag->P_tentative = ML_Operator_ArrayCreate(ag->max_levels);
         prev_P_tentatives = ag->P_tentative;
         for (jj = 0; jj < ag->max_levels; jj++) prev_P_tentatives[jj] = NULL;
       }
       prev_P_tentatives[clevel] = Pmatrix;
     }

   } /* if ( ag->smoothP_damping_factor != 0.0 ) */
   else if ( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() )
     printf("Gen_Prolongator (level %d) : not smoothing prolongator\n", level);
   StopTimer(&t1,&smooPTime);
   if ( ML_Get_PrintLevel() > 9) {
     char eigmethod[32];
#    ifdef ML_TIMING
     if (ml->comm->ML_mypid == 0)
       printf("Detailed timing for forming smoothed prolongator (level %d)\n", level);
#    endif
     ReportTimer(createPtentTime,"Create tentative P (total)",ml->comm);
     switch( Amat->spectral_radius_scheme ) {
       case ML_USE_CG:
         sprintf(eigmethod,"Eigen estimate (cg)       ");
         break;
       case ML_USE_POWER:
         sprintf(eigmethod,"Eigen estimate (power)    ");
         break;
       default:
         sprintf(eigmethod,"Eigen estimate (max norm) ");
         break;
     }
     ReportTimer(eigenTime,      eigmethod,ml->comm);
     ReportTimer(smooPTime,      "Smooth tentative P        ",ml->comm);
#    ifdef ML_TIMING
     if (ml->comm->ML_mypid == 0)
       printf("\n");
#    endif
   }
   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));

   if (widget.near_bdry != NULL) ML_free(widget.near_bdry);
   if (mls_widget != NULL) ML_Smoother_Destroy_MLS(mls_widget);

#ifdef ML_TIMING
   ml->Pmat[clevel].build_time =  GetClock() - t0;
   ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
   return 0;
} /* ML_AGG_Gen_Prolongator() */

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

int ML_AGG_JacobiSmoother_Getrows(ML_Operator *data, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[],
   double values[], int row_lengths[])
{
   struct ML_AGG_Matrix_Context *widget;
   ML_GetrowFunc  *getrow_obj;
   int            info, diag = -1, i, j /*, *aggr_info */;
   double         diag_val = 1.0, dropped, threshold = 0.0;
   ML_Operator    *mat_in;

   mat_in = (ML_Operator *) data;
   widget = (struct ML_AGG_Matrix_Context *) ML_Get_MyGetrowData(mat_in);
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

   info = getrow_obj->func_ptr(widget->Amat, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
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
/* ************************************************************************* */
/* getrow function I - w Dinv A. This one is more accurate than the standard */
/* ML one as it uses the same style of dropping as was done for the          */
/* aggregation. Additionally, the Dinv is that of filtered A.                */
/* ------------------------------------------------------------------------- */

int ML_AGG_JacobiMoreAccurate_Getrows(ML_Operator *data, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[],
   double values[], int row_lengths[])
{
   struct ML_AGG_Matrix_Context *widget;
   ML_GetrowFunc  *getrow_obj;
   int            info, diag = -1, i, j /*, *aggr_info */;
   double         diag_val = 1.0, dropped, threshold = 0.0, *thediag;
   ML_Operator    *mat_in;

   mat_in = (ML_Operator *) data;
   widget = (struct ML_AGG_Matrix_Context *) ML_Get_MyGetrowData(mat_in);
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

   info = getrow_obj->func_ptr(widget->Amat, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
   if (info == 0) return(0);

   /* ----------------------------------------------------------------- */
   /* compute threshold for dropping                                    */
   /* ----------------------------------------------------------------- */

thediag = widget->Adiag;
   if ( widget->drop_tol > 0.0 )
   {
      for (i = 0; i < -row_lengths[0]; i++)
      {
         if (columns[i] == requested_rows[0])
         {
            threshold = ML_dabs(values[i])*widget->drop_tol;
            break;
         }
      }

threshold = sqrt(ML_dabs(thediag[requested_rows[0]]))*widget->drop_tol;
      j = 0;
      dropped = 0.0;
      for (i = 0; i < row_lengths[0]; i++)
      {
         if ( ML_dabs(values[i]) >= threshold*sqrt(ML_dabs(thediag[columns[i]])))
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
diag_val += dropped;  /* rst new */

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

int ML_AGG_Amat_Getrows(ML_Operator *data, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[],
   double values[], int row_lengths[])
{
   struct ML_AGG_Matrix_Context *widget;
   ML_GetrowFunc  *getrow_obj;
   int            info;
   ML_Operator    *mat_in;

   mat_in = (ML_Operator *) data;
   widget = (struct ML_AGG_Matrix_Context *) ML_Get_MyGetrowData(mat_in);
   getrow_obj = widget->Amat->getrow;
   if (N_requested_rows > 1)
   {
      printf("Too bad. This routine only works with 1 row at a time\n");
      exit(1);
   }

   info = getrow_obj->func_ptr(widget->Amat, N_requested_rows,
			    requested_rows, allocated_space, columns,
			    values, row_lengths);
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
   int          *new_ia, *new_ja, p_ncols, *p_cols, max_nz_per_row, min_nz_per_row;
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
   int           (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void          *getrowdata;
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
   getrowfunc = getrow_obj->func_ptr;
   getrowdata = (void *) Amat;
   max_nz_per_row = 100;
   col_ind = (int *)    ML_allocate( max_nz_per_row * sizeof(int) );
   col_val = (double *) ML_allocate( max_nz_per_row * sizeof(double) );
   nnz = 0;
   for ( i = 0; i < Nfine; i++ )
   {
      while (getrowfunc((ML_Operator *) getrowdata,1,&i,max_nz_per_row,col_ind,col_val,&k)== 0)
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
   ML_Aggregate_Set_OutputLevel( newag, 0.);
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
   ML_Krylov_Set_MaxIterations(kdata, Amat->spectral_radius_max_iters);
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
   ML_memory_alloc((void**)&(new_ia), (unsigned int) nbytes, "AD1");
   nbytes = Nfine * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), (unsigned int) nbytes, "AD2");
   nbytes = Nfine * sizeof(double);
   ML_memory_alloc((void**)&(new_val), (unsigned int) nbytes, "AD3");
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
   ML_Operator_Set_ApplyFuncData(tentP,1,Nfine,csr_data,Nfine,NULL,0);
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
   ML_Operator_Set_Getrow(tentP, Nfine, CSR_getrow);
   ML_Operator_Set_ApplyFunc(tentP, CSR_matvec);

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
   ML_memory_alloc((void**)&(new_ia), (unsigned int) nbytes, "ADA");
   nbytes = nnz * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), (unsigned int) nbytes, "ADB");
   nbytes = nnz * sizeof(double);
   ML_memory_alloc((void**)&(new_val), (unsigned int) nbytes, "ADC");

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
   min_nz_per_row = 1e6;
   for ( i = 0; i < Nfine; i++ )
   {
      nnz = 0;
      for ( j = 0; j < ap_ncols; j++ )
         if ( ap_aa[j][i] != 0.0 ) nnz++;
      if ( nnz > max_nz_per_row ) max_nz_per_row = nnz;
      if ( nnz < min_nz_per_row && nnz > 0 ) min_nz_per_row = nnz;
   }

   ML_memory_alloc((void**)&csr_data,sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
   ML_Operator_Set_ApplyFuncData(APMat,1,Nfine,csr_data,
                                 Nfine,NULL,ap_ncols-1);
   APMat->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_Operator_Set_Getrow(APMat, Nfine, CSR_getrow);
   ML_Operator_Set_ApplyFunc(APMat, CSR_matvec);
   APMat->max_nz_per_row = max_nz_per_row;
   APMat->min_nz_per_row = min_nz_per_row;
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

int ML_AGG_DD_Matvec(ML_Operator *obj,int leng1,double p[],int leng2,double ap[])
{
   int         i, j, m, max_row_nnz=100, nRows, index, *col_ind;
   double      dtmp, *col_val;
   ML_Operator *Amat;
   int         (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void          *getrowdata;
   struct ML_AGG_Matrix_Context *context;
   ML_GetrowFunc                 *getrow_obj;
   ML_Operator  *mat_in;

   mat_in = (ML_Operator *) obj;
   context = (struct ML_AGG_Matrix_Context *) ML_Get_MyMatvecData(mat_in);
   Amat    = (ML_Operator *) context->Amat;
   nRows   = Amat->outvec_leng;
   if ( nRows != leng1 || leng1 != leng2 )
   {
      printf("ML_AGG_DD_Matvec ERROR : inleng != outleng.\n");
      exit(-1);
   }
   getrow_obj = Amat->getrow;
   getrowfunc = getrow_obj->func_ptr;
   getrowdata = (void *) Amat;
   if ( getrowfunc == NULL )
   {
      printf("ML_AGG_DD_Matvec ERROR : null getrowfunc.\n");
      exit(-1);
   }
   col_ind = (int *)    ML_allocate( max_row_nnz * sizeof(int) );
   col_val = (double *) ML_allocate( max_row_nnz * sizeof(double) );

   for ( i = 0; i < nRows; i++ )
   {
      while (getrowfunc((ML_Operator *) getrowdata,1,&i,max_row_nnz,col_ind,col_val,&m)== 0)
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

int ML_AGG_DD_Getrow(ML_Operator *obj,int inNrows, int *rowlist,int alloc_space,
                     int *col_ind, double *col_val, int *rowcnt)
{
   int         i, count, status, nRows, *local_ind = NULL;
   double      *local_val = NULL;
   ML_Operator *Amat;
   int         (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void        *getrowdata;
   struct ML_AGG_Matrix_Context *context;
   ML_GetrowFunc                *getrow_obj;
   ML_Operator *mat_in;

   mat_in = (ML_Operator *) obj;
   if ( inNrows != 1 )
   {
      printf("ML_AGG_DD_Getrow ERROR : inNrows > 1 not supported.\n");
      exit(-1);
   }
   context = (struct ML_AGG_Matrix_Context *) ML_Get_MyGetrowData(mat_in);
   Amat    = (ML_Operator *) context->Amat;
   nRows   = Amat->outvec_leng;
   getrow_obj = Amat->getrow;
   getrowfunc = getrow_obj->func_ptr;
   getrowdata = (void *) Amat;
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
   status = getrowfunc((ML_Operator *) getrowdata, 1, rowlist, alloc_space, local_ind,
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
   int           (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void        *getrowdata;
   ML_GetrowFunc *getrow_obj;

   nRows   = Amat->outvec_leng;
   getrow_obj = Amat->getrow;
   getrowfunc = getrow_obj->func_ptr;
   getrowdata = (void *) Amat;

   if ( getrowfunc == NULL )
   {
      printf("ML_AGG_Extract_Diag ERROR : null getrowfunc.\n");
      exit(-1);
   }
   col_ind = (int *)    ML_allocate( max_row_nnz * sizeof(int) );
   col_val = (double *) ML_allocate( max_row_nnz * sizeof(double) );

   for ( i = 0; i < nRows; i++ )
   {
      while (getrowfunc((ML_Operator *) getrowdata,1,&i,max_row_nnz,col_ind,col_val,&m)== 0)
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
   int           (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void        *getrowdata;
   ML_GetrowFunc *getrow_obj;

   local_nrows = mat->outvec_leng;
   getrow_obj = mat->getrow;
   getrowfunc = getrow_obj->func_ptr;
   getrowdata = (void *) mat;

   /* ----------------------------------------------------------------- */
   /* compute number of nonzeros                                        */
   /* ----------------------------------------------------------------- */

   max_size = 3;
   col_ind = (int *)    ML_allocate( max_size * sizeof(int) );
   col_val = (double *) ML_allocate( max_size * sizeof(double) );
   nnz = 0;
   for ( i = 0; i < local_nrows; i++ )
   {
      while (getrowfunc((ML_Operator *) getrowdata,1,&i,max_size,col_ind,col_val,&row_size)==0)
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
      getrowfunc((ML_Operator *) getrowdata,1,&i,max_size,&col_ind[nnz],&col_val[nnz],&row_size);
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
      getrowfunc((ML_Operator *) getrowdata,1,&i,max_size,col_ind,col_val,&row_size);
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
      ML_Aggregate_Set_OutputLevel( newag, 0.);
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
   ML_memory_alloc((void**)&(new_ia), (unsigned int) nbytes, "AD1");
   nbytes = Nfine * sizeof(int);
   ML_memory_alloc((void**)&(new_ja), (unsigned int) nbytes, "AD2");
   nbytes = Nfine * sizeof(double);
   ML_memory_alloc((void**)&(new_val), (unsigned int) nbytes, "AD3");
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
   ML_Operator_Set_ApplyFuncData(tentP,1,Nfine,csr_data,Nfine,NULL,0);
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
   ML_Operator_Set_Getrow(tentP, Nfine, CSR_getrow);
   ML_Operator_Set_ApplyFunc(tentP,
CSR_matvec);
   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));

/*###*/
   widget.Amat = Amat;
   widget.aggr_info = ag->aggr_info[level];
   AGGsmoother = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(AGGsmoother, widget.Amat->invec_leng,
                        widget.Amat->outvec_leng, &widget,
                        widget.Amat->matvec->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(AGGsmoother,
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
       ML_Gen_Restrictor_TransP(ml, old_mesh_level, mesh_level, NULL);
     }

     /* clean and regenerate A */

     mat = &(ml->Amat[mesh_level]);
     ML_Operator_Clean(mat);
     ML_Operator_Init(mat,ml->comm);
     ML_Gen_AmatrixRAP(ml, old_mesh_level, mesh_level);
   }

   return 0;
}

/*****************************************************************************/
/* I believe this routine works in the following way. The two dimensional    */
/* array Filter has already been computed. Filter[i][0] indicates how many   */
/* column entries in the ith row WILL be removed. Filter[i][1:Filter[i][0]]  */
/* indicates which column values should be removed.  In addition, the matrix */
/* diagonal is also modified so that piecewise constants are still in the    */
/* null space (assuming that they were originally). This code appears to work*/
/* even when we have a block PDE system.                                     */
/*****************************************************************************/

static int ML_Aux_Getrow(ML_Operator *data, int N_requested_rows, int requested_rows[],
                         int allocated_space, int columns[], double values[],
                         int row_lengths[])
{
  int ierr;
  int i, j, count, mod;
  int BlockCol, BlockRow, RowMod;
  double DiagValue = 0.0;
  int DiagID;
  int* Filter;

  ierr = (*(data->aux_data->aux_func_ptr))(data, N_requested_rows, requested_rows,
                                      allocated_space, columns, values, row_lengths);
  if (ierr == 0)
    return(0);

  if (N_requested_rows != 1) {
    fprintf(stderr, "ML_Aux_Getrow() works only if N_requested_rows == 1\n"
            "(file %s, line %d)\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  /* new part */
  mod       = data->invec_leng / data->aux_data->filter_size;
  BlockRow  = requested_rows[0] / mod;
  RowMod    = requested_rows[0] % mod;
  Filter    = data->aux_data->filter[BlockRow];
  count     = 0;
  DiagID    = -1;
  DiagValue = 0.0;

  for (i = 0 ; i < row_lengths[0] ; ++i)
  {
    BlockCol = columns[i] / mod;

    if (BlockCol ==  BlockRow)
    {
      columns[count] = columns[i];
      values[count]  = values[i];
      if (requested_rows[0] == columns[i]) DiagID = count;
      ++count;
      goto after;
    }

    /* different block col, same equation */
    for (j = 0 ; j < Filter[0] ; ++j)
    {
      /* look for elements to discard */
      if (Filter[j + 1] == BlockCol)
      {
        /* We modify the matrix diagonal so that the constant vector */
        /* is still in the null space (assuming that the constant    */
        /* was in the original matrix null space). For block pde     */
        /* systems we need to check that both the column and the row */
        /* correspond to the same DOF within the node.               */
        /* Note: This code will not preserve the null space if it is */
        /* not given by piecewise constants.                         */

        if (columns[i] % mod == RowMod) {
          DiagValue += values[i];
        }
        goto after;
      }
    }
    columns[count] = columns[i];
    values[count]  = values[i];
    ++count;
after:
    continue;
  }

  if (DiagID == -1)
  {
    fprintf(stderr, "Diagonal not defined for row %d\n", requested_rows[0]);
  }

  values[DiagID] += DiagValue;
  row_lengths[0] = count;

  return(ierr);
}

#define ML_FUNCTION_NAME "ML_Project_Coordinates"
void ML_Project_Coordinates(ML_Operator* Amat, ML_Operator* Pmat,
                            ML_Operator* Cmat)
{
  double* new_x_coord = NULL;
  double* new_y_coord = NULL;
  double* new_z_coord = NULL;
  int PDEs = Cmat->num_PDEs;
  int Nghost, i;
  ML_Operator* Rmat;
  int size_old, size_new;
  double* tmp_old,* tmp_new;
  double* aggr_sizes;
  int (*getrow)(ML_Operator*, int, int [], int, int [], double [], int []) = NULL;
  int (*matvec)(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[]) = NULL;
  ML_Aggregate_Viz_Stats *Agrid_info;
  ML_Aggregate_Viz_Stats *Cgrid_info;
  int oldPDEs=Amat->num_PDEs;

  if (PDEs != 1)
  {
    getrow = Pmat->getrow->func_ptr;
    matvec = Pmat->matvec->func_ptr;

    if (getrow != CSR_getrow && getrow != sCSR_getrows)
    {
      fprintf(stderr, "ERROR: only CSR_getrow() and sCSR_getrows() are currently supported\n"
              "ERROR: (file %s, line %d)\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    if (matvec != CSR_matvec && matvec != sCSR_matvec)
    {
      fprintf(stderr, "ERROR: only CSR_matvec() and sCSR_matvec() are currently supported\n"
              "ERROR: (file %s, line %d)\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    Pmat->getrow->func_ptr = CSR_get_one_row;
    Pmat->matvec->func_ptr = CSR_ones_matvec;
  }
  Agrid_info = (ML_Aggregate_Viz_Stats *) Amat->to->Grid->Grid;

  if (Agrid_info == NULL)
    pr_error("Amat->to->Grid->Grid == NULL\nERROR: (file %s, line %d)\n",
            __FILE__, __LINE__);

  Rmat = ML_Operator_Create(Pmat->comm);
#ifdef ML_TESTING_ONLY
  ML_Operator_Transpose(Pmat,Rmat);
  matvec = Rmat->matvec->func_ptr;
  /* Comment out the next line to check whether ML catches dividing by zero
    sized aggregate. */
  Rmat->matvec->func_ptr = CSR_ones_matvec;
#else
  ML_CommInfoOP_TransComm(Pmat->getrow->pre_comm,&(Rmat->getrow->post_comm),
			    Pmat->invec_leng);

  ML_Operator_Set_ApplyFuncData(Rmat, Pmat->outvec_leng,
                                Pmat->invec_leng,
                                Pmat->data, -1, CSR_trans_ones_matvec, 0);

  Rmat->getrow->func_ptr = NULL;
  Rmat->data_destroy = NULL;
#endif /*ifdef ML_TESTING_ONLY */

  if (Cmat->getrow->pre_comm == NULL)
    Nghost = 0;
  else {
    if (Cmat->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Cmat->getrow->pre_comm);
    Nghost = Cmat->getrow->pre_comm->total_rcv_length;
  }

  size_old = Rmat->invec_leng;
  size_new = Rmat->outvec_leng + Nghost;
  tmp_old = (double*) ML_allocate(sizeof(double) * (size_old + 1));
  tmp_new = (double*) ML_allocate(sizeof(double) * (size_new + 1));
  for (i=0; i<size_new+1; i++) tmp_new[i] = 0.0;
  aggr_sizes = (double*) ML_allocate(sizeof(double) * (size_new + 1));

  /* computes how many nodes are included in each aggregate */

  for (i = 0 ; i < size_old ; ++i)
    tmp_old[i] = 0.0;

  /*  for (i = 0 ; i < size_old ; i += PDEs)
      tmp_old[i] = 1.0;*/
  for (i = 0 ; i < size_old ; i += oldPDEs)
      tmp_old[i] = 1.0;

  ML_Operator_Apply(Rmat, Rmat->invec_leng, tmp_old, Rmat->outvec_leng, aggr_sizes);

  ML_exchange_bdry(aggr_sizes,Cmat->getrow->pre_comm,Cmat->outvec_leng,
                   Cmat->comm, ML_OVERWRITE,NULL);

  /* project the coordinates */

  Cgrid_info = (ML_Aggregate_Viz_Stats *) Cmat->to->Grid->Grid;
  if (Agrid_info->x!= NULL)
  {
    for (i = 0 ; i < size_old ; i+=oldPDEs)
      tmp_old[i] = Agrid_info->x[i / oldPDEs];
      /*      tmp_old[i] = Agrid_info->x[i / PDEs];*/

    ML_Operator_Apply(Rmat, size_old, tmp_old, Rmat->outvec_leng, tmp_new);

    new_x_coord = (double *) ML_allocate(sizeof(double) * (size_new / PDEs+1));

    ML_exchange_bdry(tmp_new,Cmat->getrow->pre_comm,Cmat->outvec_leng,
                     Cmat->comm, ML_OVERWRITE,NULL);

    for (i = 0 ; i < size_new ; i+=PDEs) {
      if (aggr_sizes[i] != 0.0)
        new_x_coord[i / PDEs] = tmp_new[i] / aggr_sizes[i];
      else {
        if (tmp_new[i] == 0.0)
          new_x_coord[i / PDEs] = 0.0;
        else {
          char msg[240];
          sprintf(msg,"(pid %d) agg %d size = %f but nonzero coordinate = %f",
                  Cmat->comm->ML_mypid, i, aggr_sizes[i],tmp_new[i]);
          pr_error("*ML_ERR* %s\n*ML_ERR* function %s\n*ML_ERR* file %s\n*ML_ERR* line %d\n",
                    msg,ML_FUNCTION_NAME,__FILE__, __LINE__);
        }
      }
    }

    if (Cgrid_info->x != NULL) ML_free(Cgrid_info->x);
    Cgrid_info->x = new_x_coord;
  }

  if (Agrid_info->y != NULL)
  {
    for (i = 0 ; i < size_old ; i+=oldPDEs)
      tmp_old[i] = Agrid_info->y[i / oldPDEs];

    /*    for (i = 0 ; i < Rmat->invec_leng ; i+=PDEs)
          tmp_old[i] = Agrid_info->y[i / PDEs];*/

    ML_Operator_Apply(Rmat, Rmat->invec_leng, tmp_old, Rmat->outvec_leng, tmp_new);

    new_y_coord = (double *) ML_allocate(sizeof(double) * (size_new / PDEs+1));

    ML_exchange_bdry(tmp_new,Cmat->getrow->pre_comm,Cmat->outvec_leng,
                     Cmat->comm, ML_OVERWRITE,NULL);

    for (i = 0 ; i < size_new ; i+=PDEs) {
      if (aggr_sizes[i] != 0.0)
        new_y_coord[i / PDEs] = tmp_new[i] / aggr_sizes[i];
      else {
        if (tmp_new[i] == 0.0)
          new_y_coord[i / PDEs] = 0.0;
        else {
          char msg[240];
          sprintf(msg,"(pid %d) agg %d size = %f but nonzero coordinate = %f",
                  Cmat->comm->ML_mypid, i, aggr_sizes[i],tmp_new[i]);
          pr_error("*ML_ERR* %s\n*ML_ERR* function %s\n*ML_ERR* file %s\n*ML_ERR* line %d\n",                    msg,ML_FUNCTION_NAME,__FILE__, __LINE__);
        }
      }
    }

    if (Cgrid_info->y != NULL) ML_free(Cgrid_info->y);
    Cgrid_info->y = new_y_coord;
  }

  if (Agrid_info->z != NULL)
  {
    for (i = 0 ; i < size_old ; i+=oldPDEs)
      tmp_old[i] = Agrid_info->z[i / oldPDEs];

    /*    for (i = 0 ; i < Rmat->invec_leng ; i+=PDEs)
	  tmp_old[i] = Agrid_info->z[i / PDEs];*/

    ML_Operator_Apply(Rmat, Rmat->invec_leng, tmp_old, Rmat->outvec_leng, tmp_new);

    new_z_coord = (double *) ML_allocate(sizeof(double) * (size_new / PDEs+1));

    ML_exchange_bdry(tmp_new,Cmat->getrow->pre_comm,Cmat->outvec_leng,
                     Cmat->comm, ML_OVERWRITE,NULL);

    for (i = 0 ; i < size_new ; i+=PDEs) {
      if (aggr_sizes[i] != 0.0)
        new_z_coord[i / PDEs] = tmp_new[i] / aggr_sizes[i];
      else {
        if (tmp_new[i] == 0.0)
          new_z_coord[i / PDEs] = 0.0;
        else {
          char msg[240];
          sprintf(msg,"(pid %d) agg %d size = %f but nonzero coordinate = %f",
                  Cmat->comm->ML_mypid, i, aggr_sizes[i],tmp_new[i]);
          pr_error("*ML_ERR* %s\n*ML_ERR* function %s\n*ML_ERR* file %s\n*ML_ERR* line %d\n",
                    msg,ML_FUNCTION_NAME,__FILE__, __LINE__);
        }
      }
    }

    if (Cgrid_info->z != NULL) ML_free(Cgrid_info->z);
    Cgrid_info->z = new_z_coord;
  }

  ML_free(tmp_old);
  ML_free(tmp_new);
  ML_free(aggr_sizes);

  /*FIXME is this right??*/
  Cgrid_info->Ndim = Agrid_info->Ndim;

  if (PDEs != 1)
  {
    Pmat->getrow->func_ptr = getrow;
    Pmat->matvec->func_ptr = matvec;
  }

#ifdef ML_TESTING_ONLY
  Rmat->matvec->func_ptr = matvec;
#endif

  ML_Operator_Destroy(&Rmat);
}
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

/*
 * This function allocates the filter field of the aux_data
 * structure. Filter[i] contains the thrown away nonzero
 * connections in amalgamated row i. Filter[i][0] is the number
 * of removed connections; connection `j' is stored in
 * Filter[0][j + 1]
 */
static void ML_Init_Aux(ML* ml, int level)
{
  int i, j, n, count, num_PDEs, BlockRow, BlockCol;
  double threshold;
  int* columns;
  double* values;
  int allocated, entries = 0;
  int N_dimensions;
  int DiagID;
  double DiagValue;
  int** filter;
  double dist;
  double* x_coord,* y_coord,* z_coord;
  ML_Aggregate_Viz_Stats *grid_info;
  double *LaplacianDiag;
  int     Nghost;

  ML_Operator* A = &(ml->Amat[level]);
  grid_info = (ML_Aggregate_Viz_Stats *) A->to->Grid->Grid;
  num_PDEs = A->num_PDEs;
  N_dimensions = grid_info->Ndim;
  x_coord = grid_info->x;
  y_coord = grid_info->y;
  z_coord = grid_info->z;
  if (x_coord == NULL) {
      printf("ML_Init_Aux: Cannot use aux options without supplying coordinates\n");
      exit(1);
  }

  threshold = A->aux_data->threshold;
  ML_Operator_AmalgamateAndDropWeak(A, num_PDEs, 0.0);

  n = A->invec_leng;
  Nghost = ML_CommInfoOP_Compute_TotalRcvLength(A->getrow->pre_comm);

  LaplacianDiag = (double *) ML_allocate((A->getrow->Nrows+Nghost+1)*
                                         sizeof(double));

  filter = (int**) ML_allocate(sizeof(int*) * n);

  allocated = 128;
  columns = (int *)    ML_allocate(allocated * sizeof(int));
  values  = (double *) ML_allocate(allocated * sizeof(double));

  for (i = 0 ; i < n ; ++i) {
    BlockRow = i;
    DiagID = -1;
    DiagValue = 0.0;

    ML_get_matrix_row(A,1,&i,&allocated,&columns,&values, &entries,0);

    for (j = 0; j < entries; j++) {
      BlockCol = columns[j];
      if (BlockRow != BlockCol) {
        dist = 0.0;
        switch (N_dimensions) {
        case 3:
          dist += (z_coord[BlockRow] - z_coord[BlockCol]) * (z_coord[BlockRow] - z_coord[BlockCol]);
        case 2:
          dist += (y_coord[BlockRow] - y_coord[BlockCol]) * (y_coord[BlockRow] - y_coord[BlockCol]);
        case 1:
          dist += (x_coord[BlockRow] - x_coord[BlockCol]) * (x_coord[BlockRow] - x_coord[BlockCol]);
        }

        if (dist == 0.0) {
          printf("node %d = %e ", i, x_coord[BlockRow]);
          if (N_dimensions > 1) printf(" %e ", y_coord[BlockRow]);
          if (N_dimensions > 2) printf(" %e ", z_coord[BlockRow]);
          printf("\n");
          printf("node %d = %e ", j, x_coord[BlockCol]);
          if (N_dimensions > 1) printf(" %e ", y_coord[BlockCol]);
          if (N_dimensions > 2) printf(" %e ", z_coord[BlockCol]);
          printf("\n");
          printf("Operator has inlen = %d and outlen = %d\n",
                 A->invec_leng, A->outvec_leng);
        }

        dist = 1.0 / dist;
        DiagValue += dist;
      }
      else if (columns[j] == i) {
        DiagID = j;
      }
    }

    if (DiagID == -1) {
      fprintf(stderr, "ERROR: matrix has no diagonal!\n"
              "ERROR: (file %s, line %d)\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    LaplacianDiag[BlockRow] = DiagValue;
  }
  if ( A->getrow->pre_comm != NULL )
     ML_exchange_bdry(LaplacianDiag,A->getrow->pre_comm,A->getrow->Nrows,
                      A->comm, ML_OVERWRITE,NULL);


  for (i = 0 ; i < n ; ++i) {
    BlockRow = i;

    ML_get_matrix_row(A,1,&i,&allocated,&columns,&values, &entries,0);

    for (j = 0; j < entries; j++) {
      BlockCol = columns[j];
      if (BlockRow != BlockCol) {
        dist = 0.0;
        switch (N_dimensions) {
        case 3:
          dist += (z_coord[BlockRow] - z_coord[BlockCol]) * (z_coord[BlockRow] - z_coord[BlockCol]);
        case 2:
          dist += (y_coord[BlockRow] - y_coord[BlockCol]) * (y_coord[BlockRow] - y_coord[BlockCol]);
        case 1:
          dist += (x_coord[BlockRow] - x_coord[BlockCol]) * (x_coord[BlockRow] - x_coord[BlockCol]);
        }

        dist = 1.0 / dist;
        values[j] = dist;
      }
    }

    count = 0;
    for (j = 0 ; j < entries ; ++j) {
      if (  (i != columns[j]) &&
            (values[j]*values[j] <
       LaplacianDiag[BlockRow]*LaplacianDiag[columns[j]]*threshold*threshold)){
        columns[count++] = columns[j];
      }
    }

    /* insert the rows */
    filter[BlockRow] = (int*) ML_allocate(sizeof(int) * (count + 1));
    filter[BlockRow][0] = count;

    for (j = 0 ; j < count ; ++j)
      filter[BlockRow][j + 1] = columns[j];

  }

  ML_free(columns);
  ML_free(values);

  ML_free(LaplacianDiag);

  ML_Operator_UnAmalgamateAndDropWeak(A, num_PDEs, 0.0);

  A->aux_data->aux_func_ptr  = A->getrow->func_ptr;
  A->getrow->func_ptr = ML_Aux_Getrow;
  A->aux_data->filter = filter;
  A->aux_data->filter_size = n;

}

static void ML_Finalize_Aux(ML* ml, const int level)
{
  int i;
  ML_Operator* A = &(ml->Amat[level]);

  A->getrow->func_ptr = A->aux_data->aux_func_ptr;
  A->aux_data->aux_func_ptr = 0;
  for (i = 0 ; i < A->aux_data->filter_size ; ++i)
    ML_free((A->aux_data->filter[i]));
  ML_free(A->aux_data->filter);
}

/* ************************************************************************* */
/* new version of ML_Gen_MGHierarchy_UsingAggregation                        */
/* I dropped out the domain_decomposition lines                              */
/* ------------------------------------------------------------------------- */

int ML_Gen_MultiLevelHierarchy_UsingAggregation(ML *ml, int start,
						int increment_or_decrement,
						ML_Aggregate *ag)
{
   int    level, idata, i;
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
   if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 10 )
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

   if (ML_Get_PrintLevel() > 10)
     for (i=0; i<level; i++) {
       int thisLevel = ml->LevelID[i];
         ML_Operator_Profile(ml->Amat+thisLevel,NULL);
         if (i != level-1) ML_Operator_Profile(ml->Rmat+thisLevel,NULL);
         if (i != 0)       ML_Operator_Profile(ml->Pmat+thisLevel,NULL);
     }

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
   int i, j, k, bail_flag, N_input_vector;
   ML_Operator *Amat, *Pmat, *Ptent;
   ML_CommInfoOP *getrow_comm;
   int aux_flag;
   ML_Aggregate_Viz_Stats *grid_info;
   ML_Aggregate *ag;
   char str[80];
   int R_is_Ptranspose;

#ifdef ML_TIMING
   double trap,    /*RAP*/
          tpart,   /*repartitioning: Zoltan, migrating Acoarse, prolongator, restrictor*/
          tprol,   /*prolongator generation (e.g., agg+QR)*/
          trest,   /*restriction generation (e.g., agg+QR)*/
          tmisc,   /*any operations not accounted for elsewhere*/
          ttotal;  /*entire level*/
#endif

   ml->ML_finest_level = fine_level;
   level = fine_level;
   next  = user_next_level(ml, level, user_data);
   k = 0;
   ml->LevelID[k++] = fine_level;

   while (next >= 0)
   {
#    ifdef ML_TIMING
     ttotal = GetClock();
#    endif
     Amat = &(ml->Amat[level]);
     if (Amat->to != NULL && Amat->to->Grid != NULL)
       grid_info =(ML_Aggregate_Viz_Stats *) Amat->to->Grid->Grid;
     else
       grid_info = NULL;

     if ( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
       printf("-----------------------------------------------------------------------\n");
       printf("ML_Gen_MultiLevelHierarchy (level %d) : Gen Restriction and Prolongator\n", level );
       printf("-----------------------------------------------------------------------\n");
     }

     aux_flag = (ml->Amat[fine_level].aux_data->enable &&
                 level <= ml->Amat[fine_level].aux_data->max_level);

     if (aux_flag)
     {
       if (ml->comm->ML_mypid == 0 && 6 < ML_Get_PrintLevel())
       {
         printf("ML_Gen_MultiLevelHierarchy (level %d) : Using auxiliary matrix\n",
                level);
         printf("ML_Gen_MultiLevelHierarchy (level %d) : threshold = %e\n",
                level, Amat->aux_data->threshold);
       }

       ML_Init_Aux(ml, level);
     }

#     ifdef ML_TIMING
      tprol = GetClock();
#     endif
      flag = (*user_gen_prolongator)(ml, level, next, user_data);
#     ifdef ML_TIMING
      tprol = GetClock() - tprol;
#     endif

      if (flag < 0) break;

      /* Now check to make sure prolongator has zero columns. */
      Pmat = ml->Pmat+next;
      bail_flag = 0;
      N_input_vector = Pmat->invec_leng;
      getrow_comm = Pmat->getrow->pre_comm;
      if (getrow_comm != NULL)
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
         if (Pmat->comm->ML_mypid == 0 && 6 < ML_Get_PrintLevel()) {
            printf("(%d) In ML_Gen_MultiLevelHierarchy: Bailing from AMG hierarchy build on level %d, where fine level = %d ........\n",
                   Pmat->comm->ML_mypid,level,fine_level);
            fflush(stdout);
         }
         if (ml->comm->ML_mypid == 0 && 6 < ML_Get_PrintLevel()) {
            printf("(%d) In ML_Gen_MultiLevelHierarchy: Nlevels = %d fine_level = %d  coarsest_level = %d\n",
               ml->comm->ML_mypid,fine_level-count+1,fine_level,count);
            fflush(stdout);
         }
         break; /* from main loop */

      }
      /* end of check */

      /* reduce memory usage */
      ML_Operator_ChangeToSinglePrecision(&(ml->Pmat[next]));

#     ifdef ML_TIMING
      trest = GetClock();
#     endif
      (*user_gen_restriction)(ml, level, next,user_data);
#     ifdef ML_TIMING
      trest = GetClock() - trest;
#     endif

      /* reduce memory usage */
      /* ML_Operator_ChangeToSinglePrecision(&(ml->Rmat[level])); */

     if (aux_flag)
     {
       ML_Finalize_Aux(ml, level);
     }

#ifdef ML_TIMING
      trap = GetClock();
#endif
      ML_Gen_AmatrixRAP(ml, level, next);
#ifdef ML_TIMING
      trap = GetClock() - trap;
#endif

      ag = (ML_Aggregate*) user_data;
      /* project the coordinates (if any) to the next coarser level */
      if (ag->P_tentative != NULL && ag->P_tentative[next] != NULL)
        Ptent = ag->P_tentative[next];
      else
        Ptent = &(ml->Pmat[next]);

      grid_info =(ML_Aggregate_Viz_Stats *) ml->Amat[level].to->Grid->Grid;
      if (grid_info) {
        if (grid_info->x != NULL) {
          if (ML_Get_PrintLevel() > 6 && ml->comm->ML_mypid == 0)
            printf("ML_Gen_MultilevelHierarchy: Projecting node coordinates from level %d to level %d\n",
                   level,next);
          ML_Project_Coordinates(ml->Amat+level, Ptent, ml->Amat+next);
        }
      }
      ML_Comm_Barrier(ml->comm);
      if (ML_Get_PrintLevel() > 10) {
        sprintf(str,"before_repartition");
        ML_Operator_Profile(ml->Amat+next,str);
      }

      R_is_Ptranspose=ML_TRUE;
      /* These values are taken from ML_MultiLevel_Gen_Prolongator().  Note that Jacob Schroder's
         energy minimization isn't checked for here, as it only runs in serial (as of 8/4/09). */
      if (ag->minimizing_energy == 1 || ag->minimizing_energy == 2 || ag->minimizing_energy == 3)
        R_is_Ptranspose = ML_FALSE;
#     ifdef ML_TIMING
      tpart = GetClock();
#     endif
      if (ML_Repartition_Get_StartLevel(ml) <= k)
        ML_repartition_Acoarse(ml, level, next, (ML_Aggregate*)user_data, R_is_Ptranspose, ML_FALSE);
      else if (ML_Get_PrintLevel() > 6 && ml->comm->ML_mypid == 0)
        printf("ML_Gen_MultilevelHierarchy (level %d): repartitioning suppressed until level %d\n",next,ML_Repartition_Get_StartLevel(ml));
#     ifdef ML_TIMING
      tpart = GetClock() - tpart;
#     endif

      if (ML_Get_PrintLevel() > 10) {
        sprintf(str,"after_repartition");
        ML_Operator_Profile(ml->Amat+next,str);
      }
      ML_Comm_Barrier(ml->comm);

#if 0
      ML_Operator_ImplicitTranspose(&(ml->Rmat[level]), &(ml->Pmat[next]), ML_TRUE);
#endif
      if (ag->P_tentative != NULL && ag->P_tentative[next] != NULL)
      {
        Ptent = ((ML_Aggregate*)user_data)->P_tentative[next];
        if (Ptent != NULL)
          ML_Operator_ChangeToSinglePrecision(Ptent);
      }
      else
        Ptent = &(ml->Pmat[next]);

      /* project the coordinates (if any) to the next
       * coarser level */
      if (grid_info) {
        if (grid_info->x != NULL || grid_info->y != NULL ||
            grid_info->z != NULL) {
          ML_Project_Coordinates(Amat, Ptent, &(ml->Amat[next]));
        }
      }

      if (aux_flag)
        ml->Amat[next].aux_data->threshold = ml->Amat[level].aux_data->threshold * 1.0;

#     ifdef ML_TIMING
      ttotal = GetClock() - ttotal;
      if ( ml->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 5 ) {
        printf("Timing summary (in seconds) for set up of level %d\n",level);
        printf("     (level %d) P creation time  = %3.2e\n", level, tprol);
        printf("     (level %d) R creation time  = %3.2e\n", level, trest);
        printf("     (level %d) RAP time         = %3.2e\n", level, trap);
        printf("     (level %d) repartition time = %3.2e\n", level, tpart);
        tmisc = ttotal - tprol - trest - trap - tpart;
        printf("     (level %d) misc. time       = %3.2e\n", level, tmisc);
        printf("     -------------------------------------\n");
        printf("     (level %d) total time       = %3.2e\n\n", level, ttotal);
      }
#     endif

      level = next;
      ml->LevelID[k++] = next;
      next  = user_next_level(ml, next, user_data);

      count++;
   } /*while (next >= 0)*/
   return(count);
} /*ML_Gen_MultiLevelHierarchy()*/

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
   ML_Operator *Amat;
   ML_Aggregate *ag = (ML_Aggregate *) data;
   struct ML_Field_Of_Values * fov;
   int flag=0; /* For the return value */
   int RelativeLevel, NumZDir, Zorientation, *LayerId = NULL, *VertLineId = NULL;
   int  Nnodes;
   struct SemiCoarsen_Struct   widget;
   void *old_field;

#if defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)
   double dtemp, dtemp2, eta;
#endif
#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   Amat     = &(ml->Amat[level]);
   Amat->num_PDEs = ag->num_PDE_eqns;

   ML_Aggregate_Set_CurrentLevel( ag, level );

   RelativeLevel = ml->ML_finest_level - level;
   if (RelativeLevel < 0) RelativeLevel = -RelativeLevel;

   /* ********************************************************************** */
   /* May require computationof field-of-values for non-diagonally scaled A  */
   /* ********************************************************************** */

   if( ag->field_of_values != NULL ) {

     fov = (struct ML_Field_Of_Values * )(ag->field_of_values);

     if( fov->compute_field_of_values_non_scaled == ML_YES ) {

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)
       ML_Anasazi_Get_FieldOfValuesBoxNonScaled_Interface(Amat,fov);
#else
       fprintf(stderr,
	       "--enable-epetra --enable-anasazi --enable-teuchos required\n"
	       "(file %s, line %d)\n",
	       __FILE__,
	       __LINE__);
       exit(EXIT_FAILURE);
#endif

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

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)
       ML_Anasazi_Get_FieldOfValuesBox_Interface(Amat,fov);
#else
	 fprintf(stderr,
		 "--enable-epetra --enable-anasazi --enable-teuchos required\n"
		 "(file %s, line %d)\n",
		 __FILE__,
		 __LINE__);
	 exit(EXIT_FAILURE);
#endif

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

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)

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
       dtemp = fov->R_coeff[0] + eta * (fov->R_coeff[1]) + pow(eta,2.) * (fov->R_coeff[2]);
       if( dtemp < 0.0 ) dtemp = 0.000001;
       dtemp2 = fov->real_max;

       ag->smoothP_damping_factor = dtemp;
       Amat->lambda_max = dtemp2;

     } else if( fov->choice == 2 ) {

       ML_Anasazi_Get_FieldOfValuesBox_Interface(Amat,fov);
       fov->eta = sqrt(pow(fov->real_max,2.) + pow(fov->imag_max,2.));

       if( ml->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	 printf("\nLargest eigenvalue (in modulus) = %e\n\n",
		fov->eta );
       }

       eta = fov->eta;
       dtemp = fov->R_coeff[0] + eta * (fov->R_coeff[1]) + pow(eta,2.) * (fov->R_coeff[2]);
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

   }

  NumZDir     = Amat->NumZDir;
  Zorientation= Amat->Zorientation;
  if ( (ml->Pmat[level]).NumZDir      != -1) NumZDir     = (ml->Pmat[level]).NumZDir;
  if ( (ml->Pmat[level]).Zorientation != -1) Zorientation= (ml->Pmat[level]).Zorientation;

  if (NumZDir == -7 )  {
     if (ml->comm->ML_mypid == 0  && ag->semicoarsen_levels != -1) {
       printf("It appears that repartitioning has been performed and so further semicoarsening is aborted.\n");
       printf("Any further line smoothing is going to numerically act as point smoothing .\n");
     }
     (ml->Pmat[level]).NumZDir = 1;
     ml->Pmat[clevel].Zorientation = 1;
     NumZDir = 1;
     Zorientation= 1;
  }
  if ( (RelativeLevel < ag->semicoarsen_levels ) && (NumZDir != 1) ) {

     Nnodes = Amat->invec_leng/Amat->num_PDEs;
     LayerId    = (int *) ML_allocate(sizeof(int)*(Nnodes+1));
     VertLineId = (int *) ML_allocate(sizeof(int)*(Nnodes+1));

     NumZDir = ML_compute_line_info(LayerId, VertLineId,Amat->invec_leng,
            Amat->num_PDEs,  ag->semicoarsen_coordinate,Zorientation, NumZDir, ml->Grid[level].Grid, ml->comm);
  }
  if ( (RelativeLevel < ag->semicoarsen_levels ) && (NumZDir > 1) ) {
     widget.nz = NumZDir;
     widget.CoarsenRate = ag->coarsen_rate;
     widget.LayerId = LayerId;
     widget.VertLineId = VertLineId;

     old_field = ag->field_of_values;
     ag->field_of_values = (void *) &widget;
     flag=ML_AGG_SemiCoarseP(ml,level, clevel, data);
     ag->field_of_values = old_field;
   }
   else {
   /* this is a bit of a hack to get the line smoothing to default to */
   /* point smoothing (block size of 1) when we are using a combination */
   /* of semi-coarsening and regular coarsening. The user should really */
   /* set the level dependent smoothers so that something sane is used  */
   /* on levels associated with regular coarsening, but if they are     */
   /* lazy we can just do line smoothing with tiny (length 1) lines.    */
   if ( (0 < ag->semicoarsen_levels ) && (NumZDir > 0) ) {
     ml->Pmat[clevel].NumZDir      = 1;
     ml->Pmat[clevel].Zorientation = 1;
   }

   /* Added on Jul-05 */
   switch (ag->minimizing_energy) {
   case -1:
   case 0:
     flag=ML_AGG_Gen_Prolongator(ml,level,clevel,data);
     break;
   case 1: /* Z_1 */
   case 2: /* Z_2 */
   case 3: /* Z_1 */
     flag=ML_AGG_Gen_Prolongator_MinEnergy(ml,level,clevel,data);
     break;
   case 4:
     flag = ML_AGG_Gen_Prolongator_MandelMinEnergy(ml,level,clevel,data);
     break;
   default:
     printf("Value of ag->minimizing_energy not correct (%d)\n"
            "(file %s, line %d)\n",
            ag->minimizing_energy,
            __FILE__,
            __LINE__);
     exit(EXIT_FAILURE);
   }
   }
   if (VertLineId != NULL)  ML_free(VertLineId);
   if (LayerId    != NULL)  ML_free(LayerId);

   return flag;
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level        */
/* number larger than the fine levels                                        */
/* ------------------------------------------------------------------------- */
int ML_MultiLevel_Gen_Restriction(ML *ml,int level, int next, void *data)
{
  ML_Operator *Amat;
  ML_Aggregate *ag = (ML_Aggregate *) data;

  struct ML_Field_Of_Values * fov;
  double dtemp, dtemp2, eta;
  char str[80];

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

    ML_Gen_Restrictor_TransP(ml, level, next, NULL);

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

    dtemp = fov->P_coeff[0] + eta * (fov->P_coeff[1]) + pow(eta,2.) * (fov->P_coeff[2]);
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

  } else {

    /* Added on Jul-05 */
    switch (ag->minimizing_energy) {
    case -1:
      if (ag->P_tentative == NULL)
      {
        fprintf(stderr, "Asked for ag->minimizing_energy == -1, but\n");
        fprintf(stderr, "the tentative prolongator has not been saved.\n");
        fprintf(stderr, "(file %s, line %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
      ML_Gen_Restrictor_TransP(ml, level, next, ag->P_tentative[next]);
      break;
    case 0:
      ML_Gen_Restrictor_TransP(ml, level, next, NULL);
      /* To do NSR we need to turn on NSR in ml_agg_min_energy.cpp */
      /* and change the above line to                              */
      /* ML_AGG_Gen_Restriction_MinEnergy(ml, level, next, data);  */
      break;
    case 1: /* Z_1 */
    case 2: /* Z_2 */
    case 3: /* Z_3 */
      ML_AGG_Gen_Restriction_MinEnergy(ml, level, next, data);
      break;
    case 4:
      ML_Gen_Restrictor_TransP(ml, level, next, NULL);
      break;
    default:
      printf("Value of ag->minimizing_energy not correct (%d)\n"
             "(file %s, line %d)\n",
             ag->minimizing_energy,
             __FILE__,
             __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  return 0;
}

/******************************************************************************
Regenerate the multigrid hierarchy using smoothed aggregation reusing the
existing aggregates.
******************************************************************************/

int  ML_Gen_MultiLevelHierarchy_UsingSmoothedAggr_ReuseExistingAgg(ML *ml,
								   ML_Aggregate *ag)
{
   int mesh_level, old_mesh_level;
   ML_Operator *mat;

   mesh_level = ml->ML_finest_level;
   if (ag->keep_P_tentative != ML_YES) {
     printf("ML_Gen_MultiLevelHierarchy_UsingSmoothedAggr_ReuseExistingAgg: must save\n");
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
       ML_MultiLevel_Gen_Prolongator(ml, old_mesh_level, mesh_level, (void*) ag);
     }

     /* clean and regenerate R */

     mat = &(ml->Rmat[old_mesh_level]);
     if (ML_Use_LowMemory() == ML_TRUE || ag->smoothP_damping_factor != 0.0) {
       ML_Operator_Clean(mat);
       ML_Operator_Init(mat,ml->comm);
       ML_MultiLevel_Gen_Restriction(ml, old_mesh_level, mesh_level, (void *)ag);
     }

     /* clean and regenerate A */

     mat = &(ml->Amat[mesh_level]);
     ML_Operator_Clean(mat);
     ML_Operator_Init(mat,ml->comm);
     ML_Gen_AmatrixRAP(ml, old_mesh_level, mesh_level);

     /* save some memory */

     ML_Operator_ChangeToSinglePrecision(&(ml->Pmat[mesh_level]));
   }

   return 0;
}

/******************************************************************************
Calculate the damping parameters for the case that the tentative prolongator is
smoothed by more than one iteration of damped Jacobi.  This corresponds to a
higher degree Chebyshev polynomial in factored form.
******************************************************************************/

void ML_AGG_Calculate_Smoothing_Factors(int numSweeps, double *factors)
{
  int i;
  int deg = 2*numSweeps+1;       /* degree of the Chebyshev poly */
  double pi=4.e0 * atan(1.e0); /* 3.141592653589793115998e0; */
  double root;

  for (i=0; i<numSweeps; i++) {
    root = cos( (2*i+1) * pi / (2*deg) );
    factors[i] = 1.0 / (root*root);
  }
}

/* This function block scales matrices. It can be invoked in differnt ways.   */
/* Let X be the matrix to be scaled and D be the block diagonal scaling       */
/* matrix. Here are the possibilities:                                        */
/*                                                                            */
/*      a)   D X (ROW_SCALE_WITH_D)                                           */
/*              Implemented by first transposing, scaling, and transposing    */
/*              back: (  X^T D^T )^T.  This was probably done so things were  */
/*              contiguous for the scaling. One REAL tricky thing is that     */
/*              ML_Operator_Transpose() only locally transposes data and then */
/*              sets up a strange post communication object. This is IMPORTANT*/
/*              because if data is shuffled (e.g. a real transpose), it would */
/*              be necessary to communicate D.  However, this is not needed as*/
/*              all data remains local. The catch is that the second transpose*/
/*              has messed up communication as ML_Operator_Transpose() was not*/
/*              designed to work with post communication on input. So we      */
/*              ignore the communication object produced by the transpose and */
/*              use X's communication object in the final result. By the way, */
/*              ML_Operator_ExplicitDinvA() looks like it has the same        */
/*              functionality (without all this transpose stuff). Probably,   */
/*              these functions were written without knowledge of each other. */
/*      b)   X D^T  (COL_SCALE_WITH_DT)                                       */
/*              X is not tranposed here. Communication is used on D.          */
/*      c)   D^T X  (ROW_SCALE_WITH_DT)                                       */
/*              Probably never tested. Should function similar to a) and so   */
/*              hopefully just works.                                         */
/*      d)   X D   (COL_SCALE_WITH_D)                                         */
/*              Probably never tested. Should function similar to b) and so   */
/*              hopefully just works.                                         */
/*                                                                            */
/* One nasty issue is that nnz(D X) >= nnz(X). That is, we might have some    */
/* new nonzeros. This means that both the data object and the communication   */
/* object need to be changed first to accommodate the new nonzeros before any */
/* scaling is applied.                                                        */

int ML_AGG_DinvP(ML_Operator *X, struct MLSthing *mls_widget,
                 int BlkSize,   int ScalingType)
/*
 * X           (On input ) Matrix to be scaled.
 *             (On output) Scaled matrix.
 *
 * mls_widget  This effectively holds the block diagonal scaling matrix.
 *             Do something like y = ML_Operator_ImplicitlyBlockDinvScale(Amat);
 *             to compute it (in y->data).
 */
{
  ML_Operator     *TmpMat, *X_or_XT;
  ML_CommInfoOP   *Blkd_comm = NULL;
  struct ML_CSR_MSRdata *csr_data;

  int      info, one = 1, **IntPtr, NBlocks, Nghost = 0,  column, last;
  int      nz_ptr, first, *new_columns = NULL, NumNnz = 0, i, j;
  int      *columns, *rowptr, **perms, block;
  double   *buffer, **DblPtr, *new_values = NULL, *values = NULL, **blockdata;
  char       N[2];
  double   *WhichDOF = NULL;
  int      *NewLocalIds =NULL, max_per_proc, CurBlock, Nneighbors, *neigh_list = NULL;
  int      Nrcv, *rcv_list, prev, *GlobalBlockIds = NULL;

  blockdata  = mls_widget->block_scaling->blockfacts;
  perms      = mls_widget->block_scaling->perms;

  if (X->getrow->pre_comm != NULL)
     Nghost = ML_CommInfoOP_Compute_TotalRcvLength(X->getrow->pre_comm);

  if ( (ScalingType == ROW_SCALE_WITH_D) || (ScalingType == COL_SCALE_WITH_DT))
     strcpy(N,"T");
  else strcpy(N,"N");

  if ((ScalingType == ROW_SCALE_WITH_D) || (ScalingType == ROW_SCALE_WITH_D)){
    X_or_XT = ML_Operator_Create(X->comm);
    ML_Operator_Transpose(X, X_or_XT); /* Do not change to by_row! */
  }
  else X_or_XT = X;

  csr_data = (struct ML_CSR_MSRdata *) X_or_XT->data;
  rowptr   = csr_data->rowptr;
  columns  = csr_data->columns;
  values   = csr_data->values;

  /* sort the blocks */
  for (i = 0; i < X_or_XT->outvec_leng; i++) {
    ML_az_sort( &(columns[rowptr[i]]), rowptr[i+1]-rowptr[i], NULL,
              &(values[rowptr[i]]));
  }

  /* See if the matrix should be padded with zeros so that the sparsity*/
  /* matches that of the final scaled matrix. Here are the 2 cases:    */
  /*     col scaling :  This is tricky because communication and       */
  /*                    column numbering needs to be changed to account*/
  /*                    for new ghost unknowns. So we do the following:*/
  /*                     a) compute new numbering for ghost nodes      */
  /*                        and put this into the unpadded matrix.     */
  /*                     b) count nonzeros in new matrix               */
  /*                     c) pad matrix with zeros                      */
  /*                     d) update the communication structure.        */
  /*     row scaling :  This is the easier case because we don't have  */
  /*                    to consider ghost variables. That is, the      */
  /*                    communication associated with X and D X is the */
  /*                    same. This also means that the node numbering  */
  /*                    is also okay. Here are the steps:              */
  /*                     a) count nonzeros in new matrix               */
  /*                     b) pad matrix with zeros                      */

  if (((ScalingType == COL_SCALE_WITH_D)||(ScalingType == COL_SCALE_WITH_DT))
       && (X->comm->ML_nprocs > 1)) {

     ML_create_unique_BlockCol_id(X->invec_leng, &GlobalBlockIds,
                BlkSize, X->getrow->pre_comm, &max_per_proc, X->comm);

     WhichDOF = (double *) ML_allocate((X->invec_leng+Nghost+1)*sizeof(double));

     for (i = 0; i < X->invec_leng; i++) WhichDOF[i] = (double) (i%BlkSize);

     ML_exchange_bdry(WhichDOF, X->getrow->pre_comm, X->invec_leng,
                        X->comm, ML_OVERWRITE, NULL);

     CurBlock = X->invec_leng/BlkSize - 1;

     NewLocalIds = (int *) ML_allocate((Nghost+1)*sizeof(int));

     Nneighbors = ML_CommInfoOP_Get_Nneighbors(X->getrow->pre_comm);
     neigh_list = ML_CommInfoOP_Get_neighbors(X->getrow->pre_comm);
     for (i = 0; i < Nneighbors; i++) {
         Nrcv = ML_CommInfoOP_Get_Nrcvlist(X->getrow->pre_comm,neigh_list[i]);
         rcv_list =ML_CommInfoOP_Get_rcvlist(X->getrow->pre_comm,neigh_list[i]);
         prev = -1;
         for (j = 0; j < Nrcv; j++) {
            if ( GlobalBlockIds[rcv_list[j]] > prev )
               CurBlock++;    /* new block */
            else if ( GlobalBlockIds[rcv_list[j]] < prev ) {
               printf("Ugh: I think this means that the send list is out of\n");
               printf("order. That is, the i^th degree of freedom within a\n");
               printf("block is sent before the j^th DOF where j<i.\n");
               exit(1);
            }

            prev = GlobalBlockIds[rcv_list[j]];

            NewLocalIds[rcv_list[j] - X->invec_leng]=
                    CurBlock*BlkSize+((int) WhichDOF[rcv_list[j]]);
         }
         if (rcv_list != NULL) ML_free(rcv_list);
     }
     if (neigh_list != NULL) ML_free(neigh_list);
     for (i = 0; i < rowptr[X->outvec_leng]; i++) {
        if ( columns[i] >= X->invec_leng) {
           columns[i] = (int) NewLocalIds[columns[i] - X->invec_leng];
        }
     }
     if (NewLocalIds    != NULL) ML_free(NewLocalIds);
     if (WhichDOF       != NULL) ML_free(WhichDOF);
     if (GlobalBlockIds != NULL) ML_free(GlobalBlockIds);
  }


  /* Count the number of nonzeros in the new local matrix */

  NumNnz = 0;
  for (i = 0; i < X_or_XT->outvec_leng; i++) {
    CurBlock = -1;
    for (j = rowptr[i]; j < rowptr[i+1]; j++) {
      if (columns[j]/BlkSize != CurBlock) {
         NumNnz += BlkSize;
         CurBlock = columns[j]/BlkSize;
      }
    }
  }

  /* Make sure that all columns within a block are present (i.e. pad with  *
   * zeros if some are missing).                                           */


  if (NumNnz != rowptr[X_or_XT->outvec_leng]) {
     new_columns= (int    *) ML_allocate(sizeof(int)*(NumNnz+1));
     new_values = (double *) ML_allocate(sizeof(double)*(NumNnz+1));
     if ( (new_columns == NULL) || (new_values == NULL)) {
        printf("ML_AGG_DinvP: Not enough memory\n"); exit(1);
     }

     first = rowptr[0];
     nz_ptr = 0;
     for (i = 0; i < X_or_XT->getrow->Nrows; i++) {
       j = first;
       last = rowptr[i+1];

       while ( j < last ) {
         block = columns[j]/BlkSize;
         column = block*BlkSize;
         while (column < block*BlkSize+BlkSize) {
           new_columns[nz_ptr] = column;
           if (  (j < rowptr[i+1]) && (columns[j] == column))
             new_values[nz_ptr++] = values[j++];
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

  /* Fix the communication object as well so it corresponds to the padded *
   * matrix. This will also be used if D needs to be communicated.        */

  if (((ScalingType == COL_SCALE_WITH_D)||(ScalingType == COL_SCALE_WITH_DT))
       && (X_or_XT->comm->ML_nprocs > 1)) {

       ML_CommInfoOP_PopulateBlks(X->getrow->pre_comm, &Blkd_comm,
                                  X->invec_leng, BlkSize,X->comm);

       if (X->getrow->pre_comm != NULL)
            ML_CommInfoOP_Destroy( &(X->getrow->pre_comm) );

       X->getrow->pre_comm = Blkd_comm;
  }

  /* Communicate D's blocks if needed */

  if ( (X_or_XT->comm->ML_nprocs > 1) &&
      ((ScalingType == COL_SCALE_WITH_DT) || (ScalingType == COL_SCALE_WITH_D)))
  {
     /* Since we scale the columns, we must communicate the block diagonals. */

     Nghost = 0;
     if (X->getrow->pre_comm != NULL)
        Nghost = ML_CommInfoOP_Compute_TotalRcvLength(X->getrow->pre_comm);

     NBlocks = (X->invec_leng + Nghost)/BlkSize;

     if (  NBlocks*BlkSize != (X->invec_leng + Nghost) ) {
        printf("invec plus ghost not divisible by blksize %d %d %d\n",
               X->invec_leng,Nghost,BlkSize);
        exit(1);
     }

     /* Allocate space for ghost blocks. First copy the nonghost block  */
     /* pointers to a new pointer-pointer array.                        */

     DblPtr = (double **)  ML_allocate(sizeof(double *)*(NBlocks+1));
     for (i = 0; i < mls_widget->block_scaling->Nblocks; i++)
        DblPtr[i] = blockdata[i];
     if (blockdata != NULL) ML_free(blockdata);
     blockdata = DblPtr;
     mls_widget->block_scaling->blockfacts = blockdata;

     IntPtr = (int **)  ML_allocate(sizeof(int *)*(NBlocks+1));
     for (i = 0; i < mls_widget->block_scaling->Nblocks; i++)
        IntPtr[i] = perms[i];
     if (perms != NULL) ML_free(perms);
     perms = IntPtr;
     mls_widget->block_scaling->perms = perms;

     /* now allocate space for ghosts */

     for (i = mls_widget->block_scaling->Nblocks; i < NBlocks; i++) {
        blockdata[i] = (double *) ML_allocate(sizeof(double)*BlkSize*BlkSize);
        if (blockdata[i] == NULL) {
           printf("Not enough space for block factors\n"); exit(1);
        }
        perms[i] = (int *) ML_allocate(sizeof(int)*BlkSize);
        if (perms[i] == NULL) {
           printf("Not enough space for perms\n"); exit(1);
        }
     }
     mls_widget->block_scaling->Nblocks = NBlocks;

    /* Extract information corresponding to the ith column within each */
    /* diagonal block. Stuff this into a vector and do an ML_exch_bdry */
    /* to communicate the column. Then, pull the ghost information out */
    /* and put it back into blockdata.                                 */


    buffer = (double *) ML_allocate(sizeof(double)*(X->invec_leng+Nghost+1));
    if (buffer == NULL) {printf("Not enough space for buffer\n"); exit(1); }

    for (i = 0; i < BlkSize; i++) {
       for (j = 0; j < X->invec_leng; j++) {
          block = j/BlkSize;
          buffer[j] = blockdata[block][i*BlkSize + j - block*BlkSize];
       }

       ML_exchange_bdry(buffer,X->getrow->pre_comm, X->getrow->Nrows,
                        X->comm,ML_OVERWRITE,NULL);

       for (j = X->invec_leng; j < X->invec_leng+Nghost; j++){
          block = j/BlkSize;
          blockdata[block][i*BlkSize+j - block*BlkSize] = buffer[j];
       }
    }
    for (j = 0; j < X->invec_leng; j++) {
       block = j/BlkSize;
       buffer[j] = (double ) perms[block][j - block*BlkSize];
    }
    ML_exchange_bdry(buffer,X->getrow->pre_comm, X->getrow->Nrows,
                     X->comm,ML_OVERWRITE,NULL);

    for (j = X->invec_leng; j < X->invec_leng+Nghost; j++) {
       block = j/BlkSize;
       perms[block][j - block*BlkSize] = (int) buffer[j];
    }
    ML_free(buffer);

  }


  if ((ScalingType == ROW_SCALE_WITH_D) || (ScalingType == COL_SCALE_WITH_DT)){
    if (mls_widget->block_scaling->optimized == 1) {
       /* optimized inv(A^T): ML_permute_for_dgetrs_special()     */
      /* must be called after the factorization is computed.      */

      for (i = 0; i < X_or_XT->getrow->Nrows; i++) {
        for (j = csr_data->rowptr[i];j < csr_data->rowptr[i+1]; j += BlkSize) {
           block  = (csr_data->columns[j])/BlkSize;
           ML_dgetrs_trans_special(BlkSize, blockdata[block], perms[block],
                                   &(values[j]));
        }
      }
    }
    else {  /* unoptimized version of inv(A^T)  */
      for (i = 0; i < X_or_XT->getrow->Nrows; i++) {
        for (j = csr_data->rowptr[i];j < csr_data->rowptr[i+1]; j += BlkSize) {
           block  = (csr_data->columns[j])/BlkSize;
           DGETRS_F77(N,&BlkSize,&one,blockdata[block],&BlkSize, perms[block],
                      &(values[j]), &BlkSize, &info);
           if ( info != 0 ) {
              printf("dgetrs returns with %d at block %d\n",info,i);
              exit(1);
           }
        }
      }
    }
  }
  else {
    if (mls_widget->block_scaling->optimized == 1) {
      /* optimized inv(A): ML_permute_for_dgetrs_special()      */
      /* must be called after the factorization is computed.    */

      for (i = 0; i < X_or_XT->getrow->Nrows; i++) {
        for (j = csr_data->rowptr[i];j < csr_data->rowptr[i+1]; j += BlkSize) {
            block  = (csr_data->columns[j])/BlkSize;
            ML_dgetrs_special(BlkSize, blockdata[block], perms[block],
                            &(values[j]));
        }
      }
    }
    else { /* unoptimized version of inv(A) */
      for (i = 0; i < X_or_XT->getrow->Nrows; i++) {
        for (j = csr_data->rowptr[i];j < csr_data->rowptr[i+1]; j += BlkSize) {
            block  = (csr_data->columns[j])/BlkSize;
            DGETRS_F77(N,&BlkSize,&one,blockdata[block],&BlkSize,
                       perms[block], &(values[j]), &BlkSize, &info);
            if ( info != 0 ) {
              printf("dgetrs returns with %d at block %d\n",info,i);
              exit(1);
            }
        }
      }
    }
  }

  /* Transpose back if needed. Must use X->getrow->pre_comm as  */
  /* ML_Operator_Transpose(ML_Operator_Transpose(X)) messes up  */
  /* the communication object. This is because the first        */
  /* transpose creates a post communication object which the    */
  /* second transpose does not know how to deal with.           */

  if ((ScalingType == ROW_SCALE_WITH_D)||(ScalingType == ROW_SCALE_WITH_DT)) {
    TmpMat = ML_Operator_Create(X->comm);
    i = X_or_XT->outvec_leng;  /* need to make ML_Operator_Transpose work */
    X_or_XT->outvec_leng = X_or_XT->getrow->Nrows;
    ML_Operator_Transpose(X_or_XT, TmpMat); /* Do not change to by_row! */
    X_or_XT->outvec_leng = i;
    ML_Operator_Destroy(&X_or_XT);

    /* swap data pointers */
    csr_data     = (struct ML_CSR_MSRdata *) TmpMat->data;
    TmpMat->data = X->data;
    X->data      = csr_data;

    ML_Operator_Destroy(&TmpMat);
  }
  return 0;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
int ML_AGG_SemiCoarseP(ML *ml,int level, int clevel, void *data)
{
  /* See comments for MakeSemiCoarsenP */

  int                 Ncoarse;
  struct SemiCoarsen_Struct  *widget;
  struct  ML_CSR_MSRdata *csr_data;
  ML_Operator *Pmatrix;
  double *Pvals;
  int    *Pptr, *Pcols, Nglobal, Ncglobal;

  ML_Aggregate * ag = (ML_Aggregate *) data;

  ML_Operator* Amat = &(ml->Amat[level]);
#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  Pmatrix = &(ml->Pmat[clevel]);
  Amat->num_PDEs    = ag->num_PDE_eqns;
  widget            = (struct SemiCoarsen_Struct *) ag->field_of_values;

  // int NVertLines = Amat->invec_leng/(ag->num_PDE_eqns*widget->nz);
  Ncoarse = MakeSemiCoarsenP(Amat->invec_leng/ag->num_PDE_eqns, widget->nz,
                   widget->CoarsenRate, widget->LayerId, widget->VertLineId,
                   ag->num_PDE_eqns, Amat, &Pptr, &Pcols, &Pvals);

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  csr_data->rowptr  = Pptr;
  csr_data->columns = Pcols;
  csr_data->values  = Pvals;
  ML_Operator_Set_ApplyFuncData(Pmatrix,Ncoarse,Amat->invec_leng,csr_data,Amat->invec_leng,NULL,0);
  Pmatrix->data_destroy = ML_CSR_MSRdata_Destroy;
  ML_Operator_Set_Getrow(Pmatrix, Amat->invec_leng, CSR_getrow);
  ML_Operator_Set_ApplyFunc(Pmatrix, CSR_matvec);

  Pmatrix->Zorientation= 1;                                  /* for the next level coarsening      */
  Pmatrix ->num_PDEs    = ag->num_PDE_eqns;

  ML_Operator_Set_1Levels(&(ml->Pmat[clevel]), &(ml->SingleLevel[clevel]),
                          &(ml->SingleLevel[level]));

  Nglobal = Amat->invec_leng;
  Ncglobal = Ncoarse;
  Nglobal = ML_Comm_GsumInt( ml->comm, Nglobal);
  Ncglobal= ML_Comm_GsumInt( ml->comm, Ncglobal);

  /* doing it this way to avoid overflow (and to handle empty procs) */
  Pmatrix->NumZDir = (int) ( (((double)Ncglobal)/((double)Nglobal))*((double) widget->nz) + .001);
  if (ml->comm->ML_mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) {
       printf("SemiCoarsening: Coarsening from %d to %d\n",Nglobal,Ncglobal);
  }
  ag->curr_threshold = ag->threshold;


#ifdef ML_TIMING
  ml->Pmat[clevel].build_time   =  GetClock() - t0;
  ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
  return(0);

}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int FindCpts(int PtsPerLine, int CoarsenRate, int Thin, int **LayerCpts)
{
/*
 * Given the number of points in the z direction (PtsPerLine) and a
 * coarsening rate (CoarsenRate), determine which z-points will serve
 * as Cpts and return the total number of Cpts.
 *
 * Input
 *    PtsPerLine:   Number of fine level points in the z direction
 *
 *    CoarsenRate:  Roughly, number of Cpts  = (PtsPerLine+1)/CoarsenRate - 1
 *
 *    Thin:         Must be either 0 or 1. Thin decides what to do when
 *                  (PtsPerLine+1)/CoarsenRate is not an integer.
 *
 *                    Thin == 0  ==>   ceil() the above fraction
 *                    Thin == 1  ==>   floor() the above fraction
 *
 * Output
 *    LayerCpts     Array where LayerCpts[i] indicates that the
 *                  LayerCpts[i]th fine level layer is a Cpt Layer.
 *                  Note: fine level layers are assumed to be numbered starting
 *                        a one.
 */
   double temp, RestStride, di;
   int    NCpts, i;
   int    NCLayers = -1;
   int    FirstStride;

   temp =  ((double) (PtsPerLine+1))/((double) (CoarsenRate)) - 1.0;
   if  (Thin == 1) NCpts = (int) ceil(temp);
   else            NCpts = (int) floor(temp);

   if (PtsPerLine == 1) { printf("cannot coarsen further\n"); return -1; }
   if (NCpts < 1) NCpts = 1;



   FirstStride= (int) ceil( ((double) PtsPerLine+1)/( (double) (NCpts+1)));
   RestStride = ((double) (PtsPerLine-FirstStride+1))/((double) NCpts);

   NCLayers   = (int) floor((((double) (PtsPerLine-FirstStride+1))/RestStride)+.00001);

   if ( NCLayers != NCpts) { printf("sizes do not match %d %d\n",NCpts,NCLayers); exit(1); }
   *LayerCpts = (int *) malloc((NCLayers+1)*sizeof(int));

   di  = (double) FirstStride;
   for (i = 1; i <= NCpts; i++) {
      (*LayerCpts)[i] = (int) floor(di);
      di += RestStride;
   }

   return(NCLayers);
}


#define MaxHorNeighborNodes 75
#define HOOKED_TO_ML
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


int MakeSemiCoarsenP(int Ntotal, int nz, int CoarsenRate, int LayerId[],
                     int VertLineId[], int DofsPerNode,
#ifdef HOOKED_TO_ML
                     ML_Operator *Amat,
#else
                     int *OrigARowPtr, int *OrigAcols, double *OrigAvals,
#endif
                     int **ParamPptr, int **ParamPcols, double **ParamPvals)
{
/*
 * Given a CSR matrix (OrigARowPtr, OrigAcols, OrigAvals), information
 * describing the z-layer and vertical line (LayerId and VertLineId)
 * of each matrix block row, a coarsening rate, and dofs/node information,
 * construct a prolongator that coarsening to semicoarsening in the z-direction
 * using something like an operator-dependent grid transfer. In particular,
 * matrix stencils are collapsed to vertical lines. Thus, each vertical line
 * gives rise to a block tridiagonal matrix. BlkRows corresponding to
 * Cpts are replaced by identity matrices. This tridiagonal is solved
 * to determine each interpolation basis functions. Each Blk Rhs corresponds
 * to all zeros except at the corresponding C-pt which has an identity
 *
 * On termination, return the number of local prolongator columns owned by
 * this processor.
 *
 * Note: This code was adapted from a matlab code where offsets/arrays
 *       start from 1. In most parts of the code, this 1 offset is kept
 *       (in some cases wasting the first element of the array). The
 *       input and output matrices of this function has been changed to
 *       have offsets/rows/columns which start from 0. LayerId[] and
 *       VertLineId[] currently start from 1.
 *
 * Input
 * =====
 *    Ntotal       Number of fine level Blk Rows owned by this processor
 *
 *    nz           Number of vertical layers. Note: partitioning must be done
 *                 so that each processor owns an entire vertical line. This
 *                 means that nz is the global number of layers, which should
 *                 be equal to the local number of layers.
 *    CoarsenRate  Rate of z semicoarsening. Smoothed aggregation-like coarsening
 *                 would correspond to CoarsenRate = 3.
 *    LayerId      Array from 0 to Ntotal-1 + Ghost. LayerId(BlkRow) gives the
 *                 layer number associated with the dofs within BlkRow.
 *    VertLineId   Array from 1 to Ntotal, VertLineId(BlkRow) gives a unique
 *                 vertical line id (from 0 to Ntotal/nz-1) of BlkRow. All
 *                 BlkRows associated with nodes along the same vertical line
 *                 in the mesh should have the same LineId.
 *    DofsPerNode  Number of degrees-of-freedom per mesh node.
 *
 *    OrigARowPtr, CSR arrays corresponding to the fine level matrix.
 *    OrigAcols,
 *    OrigAvals
 *
 * Output
 * =====
 *    ParamPptr,   CSR arrays corresponding to the final prolongation matrix.
 *    ParamPcols,
 *    ParamsPvals
 */
 int    NLayers, NVertLines, MaxNnz, NCLayers, MyLine, MyLayer;
 int    *InvLineLayer=NULL, *CptLayers=NULL, StartLayer, NStencilNodes;
 int    BlkRow, dof_i, dof_j, node_k, *Sub2FullMap=NULL, *Acols=NULL, RowLeng;
 int    i, j, iii, col, count, index, loc, PtRow, PtCol;
 double *Avals=NULL, *BandSol=NULL, *BandMat=NULL, TheSum;
 int    *IPIV=NULL, KL, KU, KLU, N, NRHS, LDAB,INFO;
 int    *Pcols, *Pptr;
 double *Pvals, *dtemp = NULL;
 int    MaxStencilSize, MaxNnzPerRow;
 int    *LayDiff=NULL;
 int    CurRow, LastGuy = -1, NewPtr;
 int    allocated = 0, Ndofs;
 int    Nghost;
 int    *Layerdofs = NULL, *Col2Dof = NULL;


  MaxNnzPerRow = MaxHorNeighborNodes*DofsPerNode*3;
  LayDiff = (int *) malloc( (1+MaxNnzPerRow)*sizeof(int));

  *ParamPptr = NULL;
  *ParamPcols= NULL;
  *ParamPvals= NULL;

   Nghost = 0;
#ifdef HOOKED_TO_ML
   Nghost = ML_CommInfoOP_Compute_TotalRcvLength(Amat->getrow->pre_comm);
#endif
   dtemp    = (double *) malloc(sizeof(double)*(Ntotal*DofsPerNode+Nghost+1));
   Layerdofs= (int *) malloc(sizeof(int)*(Ntotal*DofsPerNode+Nghost+1));
   Col2Dof  = (int *) malloc(sizeof(int)*(Ntotal*DofsPerNode+Nghost+1));
   for (i = 0; i < Ntotal*DofsPerNode; i++)
      dtemp[i]= (double)LayerId[i/DofsPerNode];
#ifdef HOOKED_TO_ML
   ML_exchange_bdry(dtemp,Amat->getrow->pre_comm,Amat->outvec_leng,
                    Amat->comm, ML_OVERWRITE,NULL);
#endif
   for (i = 0; i < Ntotal*DofsPerNode+Nghost; i++) Layerdofs[i]=dtemp[i];
   for (i = 0; i < Ntotal*DofsPerNode;        i++) dtemp[i]= i%DofsPerNode;
#ifdef HOOKED_TO_ML
   ML_exchange_bdry(dtemp,Amat->getrow->pre_comm,Amat->outvec_leng,
                    Amat->comm, ML_OVERWRITE,NULL);
#endif
   for (i = 0; i < Ntotal*DofsPerNode+Nghost; i++) Col2Dof[i]=dtemp[i];
   if (dtemp != NULL) free(dtemp);

   if (Ntotal != 0) {
     NLayers   = LayerId[0];
     NVertLines= VertLineId[0];
   }
   else { NLayers = -1; NVertLines = -1; }

  for (i = 1; i < Ntotal; i++) {
      if ( VertLineId[i] > NVertLines ) NVertLines = VertLineId[i];
      if ( LayerId[i]    >   NLayers  ) NLayers    = LayerId[i];
  }
  NLayers++;
  NVertLines++;

 /*
  * Make an inverse map so that we can quickly find the dof
  * associated with a particular vertical line and layer.
  */

  InvLineLayer = (int *) malloc( (1+NVertLines*NLayers)*sizeof(int));
  for (i=0; i < Ntotal; i++) {
     InvLineLayer[ VertLineId[i]+1+LayerId[i]*NVertLines ] = i;
  }

 /*
  * Determine coarse layers where injection will be applied.
  */

  NCLayers = FindCpts(nz,CoarsenRate,0, &CptLayers);

 /*
  * Compute the largest possible interpolation stencil width based
  * on the location of the Clayers. This stencil width is actually
  * nodal (i.e. assuming 1 dof/node). To get the true max stencil width
  * one needs to multiply this by DofsPerNode.
  */

  if  (NCLayers < 2) MaxStencilSize = nz;
  else MaxStencilSize = CptLayers[2];

  for (i = 3; i <= NCLayers; i++) {
     if (MaxStencilSize < CptLayers[i]- CptLayers[i-2])
         MaxStencilSize = CptLayers[i]- CptLayers[i-2];
  }
  if (NCLayers > 1) {
     if (MaxStencilSize < nz - CptLayers[NCLayers-1]+1)
       MaxStencilSize =  nz - CptLayers[NCLayers-1]+1;
  }

 /*
  * Allocate storage associated with solving a banded sub-matrix needed to
  * determine the interpolation stencil. Note: we compute interpolation stencils
  * for all dofs within a node at the same time, and so the banded solution
  * must be large enough to hold all DofsPerNode simultaneously.
  */

  Sub2FullMap= (int    *) malloc(sizeof(int   )*(MaxStencilSize+1)*DofsPerNode);
  BandSol    = (double *) malloc(sizeof(double)*(MaxStencilSize+1)*DofsPerNode*
                                                 DofsPerNode);
 /*
  * Lapack variables. See comments for dgbsv().
  */
  KL     = 2*DofsPerNode-1;
  KU     = 2*DofsPerNode-1;
  KLU    = KL+KU;
  LDAB   = 2*KL+KU+1;
  NRHS = DofsPerNode;
  BandMat=(double *) malloc(sizeof(double)*(LDAB*MaxStencilSize*DofsPerNode+1));
  IPIV   =(int    *) malloc(sizeof(int   )*(MaxStencilSize+1)*DofsPerNode);

 /*
  * Allocate storage for the final interpolation matrix. Note: each prolongator
  * row might have entries corresponding to at most two nodes.
  * Note: the total fine level dofs equals DofsPerNode*Ntotal and the max
  *       nnz per prolongator row is DofsPerNode*2.
  */

  Ndofs  = DofsPerNode*Ntotal;
  MaxNnz = 2*DofsPerNode*Ndofs;
  Pvals  = (double *) malloc( (1+MaxNnz)*sizeof(double));
  Pptr   = (int    *) malloc( DofsPerNode*(2+Ntotal)*sizeof(int   ));
  Pptr[0] = 0; Pptr[1] = 0;
  Pcols  = (int    *) malloc( (1+MaxNnz)*sizeof(int   ));

  if (Pcols == NULL) {
     printf("MakeSemiCoarsenP: Not enough space \n");
     if (allocated       != 0) free(Avals);
     if (allocated       != 0) free(Acols);
     if (InvLineLayer != NULL) free(InvLineLayer);
     if (CptLayers    != NULL) free(CptLayers);
     if (Sub2FullMap  != NULL) free(Sub2FullMap);
     if (BandSol      != NULL) free(BandSol);
     if (BandMat      != NULL) free(BandMat);
     if (IPIV         != NULL) free(IPIV);
     if (LayDiff      != NULL) free(LayDiff);
     if (Layerdofs    != NULL) free(Layerdofs);
     if (Col2Dof      != NULL) free(Col2Dof);
     if (Pvals        != NULL) free(Pvals);
     if (Pptr         != NULL) free(Pptr);
     return -1;
  }

 /*
  * Setup P's rowptr as if each row had its maximum of 2*DofsPerNode nonzeros.
  * This will be useful while filling up P, and then later we will squeeze out
  * the unused nonzeros locations.
  */

  for (i = 1; i <= MaxNnz; i++) Pcols[i] = -1;  /* mark all entries as unused */
  count = 1;
  for (i = 1; i <= DofsPerNode*Ntotal+1; i++) {
     Pptr[i]  = count;
     count   += (2*DofsPerNode);
  }

 /*
  * Build P column by column. The 1st block column corresponds to the 1st coarse
  * layer and the first line. The 2nd block column corresponds to the 2nd coarse
  * layer and the first line. The NCLayers+1 block column corresponds to the
  * 1st coarse layer and the 2nd line, etc.
  */

  col = 0;
  for (MyLine=1; MyLine <= NVertLines; MyLine += 1) {
    for (iii=1; iii <= NCLayers;  iii+= 1) {
      col = col+1;
      MyLayer = CptLayers[iii];

      /*
       * StartLayer gives the layer number of the lowest layer that
       * is nonzero in the interpolation stencil that is currently
       * being computed. Normally, if we are not near a boundary this
       * is simply CptsLayers[iii-1]+1.
       *
       * NStencilNodes indicates the number of nonzero nodes in the
       * interpolation stencil that is currently being computed. Normally,
       * if we are not near a boundary this is CptLayers[iii+1]-StartLayer.
       */

      if (iii !=    1    ) StartLayer = CptLayers[iii-1]+1;
      else                 StartLayer = 1;

      if (iii != NCLayers) NStencilNodes= CptLayers[iii+1]-StartLayer;
      else                 NStencilNodes= NLayers - StartLayer + 1;


      N = NStencilNodes*DofsPerNode;

      /*
       *  dgbsv() does not require that the first KL rows be initialized,
       *  so we could avoid zeroing out some entries?
       */

       for (i = 0; i < NStencilNodes*DofsPerNode*DofsPerNode; i++)
          BandSol[ i] = 0.0;
       for (i = 0; i < LDAB*N; i++) BandMat[ i] = 0.0;

      /*
       *  Fill BandMat and BandSol (which is initially the rhs) for each
       *  node in the interpolation stencil that is being computed.
       */

      for (node_k=1; node_k <= NStencilNodes ; node_k++) {

         /*  Map a Line and Layer number to a BlkRow in the fine level  matrix
          *  and record the mapping from the sub-system to the BlkRow of the
          *  fine level matrix.
          */
         BlkRow  = InvLineLayer[MyLine+(StartLayer+node_k-2)*NVertLines]+1;
         Sub2FullMap[node_k] = BlkRow;

         /* Two cases:
          *    1) the current layer is not a Cpoint layer. In this case we
          *       want to basically stick the matrix couplings to other
          *       nonzero stencil rows into the band matrix. One way to do
          *       this is to include couplings associated with only MyLine
          *       and ignore all the other couplings. However, what we do
          *       instead is to sum all the coupling at each layer participating
          *       in this interpolation stencil and stick this sum into BandMat.
          *    2) the current layer is a Cpoint layer and so we
          *       stick an identity block in BandMat and rhs.
          */
         if (StartLayer+node_k-1 != MyLayer) {
            for (dof_i=0; dof_i < DofsPerNode; dof_i++) {

#ifdef HOOKED_TO_ML
                j = (BlkRow-1)*DofsPerNode+dof_i;
                ML_get_matrix_row(Amat,1,&j,&allocated,&Acols,&Avals,&RowLeng,0);
#else
                Acols = &(OrigAcols[OrigARowPtr[(BlkRow-1)*DofsPerNode+dof_i]]);
                Avals = &(OrigAvals[OrigARowPtr[(BlkRow-1)*DofsPerNode+dof_i]]);
                RowLeng= OrigARowPtr[(BlkRow-1)*DofsPerNode+dof_i+1] -
                                    OrigARowPtr[(BlkRow-1)*DofsPerNode+dof_i];
#endif

                if (RowLeng >= MaxNnzPerRow) {
                   printf("MakeSemiCoarsenP: recompile with larger Max(HorNeighborNodes)\n");
                    if (allocated       != 0) free(Avals);
                    if (allocated       != 0) free(Acols);
                    if (InvLineLayer != NULL) free(InvLineLayer);
                    if (CptLayers    != NULL) free(CptLayers);
                    if (Sub2FullMap  != NULL) free(Sub2FullMap);
                    if (BandSol      != NULL) free(BandSol);
                    if (BandMat      != NULL) free(BandMat);
                    if (IPIV         != NULL) free(IPIV);
                    if (LayDiff      != NULL) free(LayDiff);
                    if (Layerdofs    != NULL) free(Layerdofs);
                    if (Col2Dof      != NULL) free(Col2Dof);
                    if (Pvals        != NULL) free(Pvals);
                    if (Pptr         != NULL) free(Pptr);
                    if (Pcols        != NULL) free(Pcols);
                    return -1;
                }

                for (i = 0; i < RowLeng; i++) {
                   LayDiff[i]  = Layerdofs[Acols[i]]-StartLayer-node_k+2;

                   /* This is the main spot where there might be off- */
                   /* processor communication. That is, when we       */
                   /* average the stencil in the horizontal direction,*/
                   /* we need to know the layer id of some            */
                   /* neighbors that might reside off-processor.      */
                }
                PtRow = (node_k-1)*DofsPerNode+dof_i+1;
                for (dof_j=0; dof_j < DofsPerNode; dof_j++) {
                   PtCol = (node_k-1)*DofsPerNode+dof_j + 1;
                   /* Stick in entry corresponding to Mat(PtRow,PtCol) */
                   /* see dgbsv() comments for matrix format.          */
                   TheSum = 0.0;
                   for (i = 0; i < RowLeng; i++) {
                     if ((LayDiff[i] == 0)  && (Col2Dof[Acols[i]] == dof_j))
                          TheSum += Avals[i];
                   }
                   index = LDAB*(PtCol-1)+KLU+PtRow-PtCol;
                   BandMat[index] = TheSum;

                   if (node_k != NStencilNodes) {
                      /* Stick Mat(PtRow,PtCol+DofsPerNode) entry  */
                      /* see dgbsv() comments for matrix format.  */
                      TheSum = 0.0;
                      for (i = 0; i < RowLeng; i++) {
                         if ((LayDiff[i] == 1) &&(Col2Dof[Acols[i]]== dof_j))
                             TheSum += Avals[i];
                      }
                      j = PtCol+DofsPerNode;
                      index=LDAB*(j-1)+KLU+PtRow-j;
                      BandMat[index] = TheSum;

                   }
                   if (node_k != 1) {
                      /* Stick Mat(PtRow,PtCol-DofsPerNode) entry  */
                      /* see dgbsv() comments for matrix format.  */
                      TheSum = 0.0;
                      for (i = 0; i < RowLeng; i++) {
                         if ((LayDiff[i]== -1) &&(Col2Dof[Acols[i]]== dof_j))
                             TheSum += Avals[i];
                      }
                      j = PtCol-DofsPerNode;
                      index=LDAB*(j-1)+KLU+PtRow-j;
                      BandMat[index] = TheSum;
                   }
                }
            }
         }
         else {
             for (dof_i = 0; dof_i < DofsPerNode; dof_i++) {
                /* Stick Mat(PtRow,PtRow) and Rhs(PtRow,dof_i+1) */
                /* see dgbsv() comments for matrix format.     */
                PtRow = (node_k-1)*DofsPerNode+dof_i+1;
                index = LDAB*(PtRow-1)+KLU;
                BandMat[index] = 1.0;
                BandSol[(dof_i)*DofsPerNode*NStencilNodes+PtRow-1] = 1.;
             }
         }
      }

      /* Solve banded system and then stick result in Pmatrix arrays */

      DGBSV_F77( &N, &KL, &KU, &NRHS, BandMat, &LDAB, IPIV, BandSol, &N, &INFO );

      for (dof_j=0; dof_j < DofsPerNode; dof_j++) {
         for (dof_i=0; dof_i < DofsPerNode; dof_i++) {
            for (i =1; i <= NStencilNodes ; i++) {
               index = (Sub2FullMap[i]-1)*DofsPerNode+dof_i+1;
               loc = Pptr[index];
               Pcols[loc] = (col-1)*DofsPerNode+dof_j+1;
               Pvals[loc] = BandSol[dof_j*DofsPerNode*NStencilNodes +
                                           (i-1)*DofsPerNode + dof_i];
               Pptr[index]= Pptr[index] + 1;
            }
         }
      }
    }
  }

  if (allocated       != 0) free(Avals);
  if (allocated       != 0) free(Acols);
  if (InvLineLayer != NULL) free(InvLineLayer);
  if (CptLayers    != NULL) free(CptLayers);
  if (Sub2FullMap  != NULL) free(Sub2FullMap);
  if (BandSol      != NULL) free(BandSol);
  if (BandMat      != NULL) free(BandMat);
  if (IPIV         != NULL) free(IPIV);
  if (LayDiff      != NULL) free(LayDiff);
  if (Layerdofs    != NULL) free(Layerdofs);
  if (Col2Dof      != NULL) free(Col2Dof);

 /*
  * Squeeze the -1's out of the columns. At the same time convert Pcols
  * so that now the first column is numbered '0' as opposed to '1'.
  * Also, the arrays Pcols and Pvals should now use the zeroth element
  * as opposed to just starting with the first element. Pptr will be
  * fixed in the for loop below so that Pptr[0] = 0, etc.
  */
  CurRow = 1;
  NewPtr = 1;
  for (i=1; i <= Pptr[Ntotal*DofsPerNode]-1; i++) {
     if (i == Pptr[CurRow]) {
        Pptr[CurRow] = LastGuy;
        CurRow++;
        while (i > Pptr[CurRow]) {
           Pptr[CurRow] = LastGuy;
           CurRow++;
        }
     }
     if (Pcols[i] != -1) {
        Pcols[NewPtr-1] = Pcols[i]-1;   /* these -1's fix the offset and */
        Pvals[NewPtr-1] = Pvals[i];     /* start using the zeroth element */
        LastGuy = NewPtr;
        NewPtr++;
     }
  }
  for (i = CurRow; i <= Ntotal*DofsPerNode; i++) Pptr[CurRow] = LastGuy;

 /* Now move the pointers so that they now point to the beginning of each
  * row as opposed to the end of each row
  */
  for (i=-Ntotal*DofsPerNode+1; i>= 2 ; i--) {
     Pptr[i-1] = Pptr[i-2];  /* extra -1 added to start from 0 */
  }
  Pptr[0] = 0;

  Pcols = (int    *) realloc(Pcols, (LastGuy+1)*sizeof(int));
  Pvals = (double *) realloc(Pvals, (LastGuy+1)*sizeof(double));

  *ParamPptr  = Pptr;
  *ParamPcols = Pcols;
  *ParamPvals = Pvals;

  return NCLayers*NVertLines*DofsPerNode;
}

int ML_compute_line_info(int LayerId[], int VertLineId[],
                         int Ndof, int DofsPerNode, char semicoarsen_coordinate,
                         int MeshNumbering, int NumNodesPerVertLine,
                         ML_Aggregate_Viz_Stats *grid_info, ML_Comm *comm)
{
   double *xvals= NULL, *yvals = NULL, *zvals = NULL;
   int    Nnodes, NVertLines, MyNode;
   int    NumCoords, NumBlocks, index, next, subindex, subnext;
   double xfirst, yfirst;
   double *xtemp, *ytemp, *ztemp;
   int    *OrigLoc;
   int    i,j,count;
   int    RetVal, gRetVal;
   int    mypid;

   mypid = comm->ML_mypid;
   RetVal = 0;
   if ((MeshNumbering != 1) && (MeshNumbering != 2)) {
      if (grid_info != NULL) xvals = grid_info->x;
      if (grid_info != NULL) yvals = grid_info->y;
      if (grid_info != NULL) zvals = grid_info->z;

      ztemp = zvals;
      if (semicoarsen_coordinate == 'x') 
         { zvals = xvals; if (ztemp == NULL) xvals=yvals; else xvals= ztemp;}
      if (semicoarsen_coordinate == 'y') 
         { zvals = yvals; if (ztemp == NULL) yvals=xvals; else yvals= ztemp;}
      if (semicoarsen_coordinate == 'z') {
        if ( (zvals == NULL) && (xvals != NULL) && (yvals != NULL) ) {
         printf("Cannot coarsen 2D problems in z direction. Must set semicoarsen_coordinate to x or y\n");
         exit(1);
        }
      }

      if ( (xvals == NULL) || (yvals == NULL) || (zvals == NULL)) RetVal = -1;
   }
   else {
      if  (NumNodesPerVertLine == -1)                     RetVal = -4;
      if ( ((Ndof/DofsPerNode)%NumNodesPerVertLine) != 0) RetVal = -3;
   }
   if ( (Ndof%DofsPerNode) != 0) RetVal = -2;

   gRetVal = ML_gmax_int(RetVal, comm);
   if ( gRetVal < 0)  {
      i = comm->ML_nprocs;
      if (RetVal < 0) i = mypid;
      j = ML_gmin_int(i, comm);

      if (mypid == j) {
         if (RetVal == -1) printf("Not semicoarsening as no mesh numbering information or coordinates are given\n");
         if (RetVal == -4) printf("Not semicoarsening as the number of z nodes is not given.\n");
         if (RetVal == -3) printf("Not semicoarsening as the total number of nodes is not evenly divisible by the number of z direction nodes .\n");
         if (RetVal == -2) printf("Not semicoarsening as something is off with the number of degrees-of-freedom per node.\n");
      }
      return gRetVal;
   }

   Nnodes = Ndof/DofsPerNode;

   for (MyNode = 0; MyNode < Nnodes;  MyNode++) VertLineId[MyNode]= -1;
   for (MyNode = 0; MyNode < Nnodes;  MyNode++) LayerId[MyNode]   = -1;


   if (MeshNumbering == 1) {
      for (MyNode = 0; MyNode < Nnodes; MyNode++) {
         LayerId[MyNode]= MyNode%NumNodesPerVertLine;
         VertLineId[MyNode]= (MyNode- LayerId[MyNode])/NumNodesPerVertLine;
      }
   }
   else if (MeshNumbering == 2) {
      NVertLines = Nnodes/NumNodesPerVertLine;
      for (MyNode = 0; MyNode < Nnodes; MyNode++) {
         VertLineId[MyNode]   = MyNode%NVertLines;
         LayerId[MyNode]   = (MyNode- VertLineId[MyNode])/NVertLines;
      }
   }
   else {


          NumCoords = Ndof/DofsPerNode;

          /* sort coordinates so that we can order things according to lines */

          OrigLoc = (int    *) ML_allocate(sizeof(int   )*(NumCoords+1));
          xtemp   = (double *) ML_allocate(sizeof(double)*(NumCoords+1));
          ytemp   = (double *) ML_allocate(sizeof(double)*(NumCoords+1));
          ztemp   = (double *) ML_allocate(sizeof(double)*(NumCoords+1));

          if (ztemp == NULL) {
             printf("Not enough memory for line smoothers\n");
             exit(EXIT_FAILURE);
          }
          for (i = 0; i < NumCoords; i++) ytemp[i]= yvals[i];
          for (i = 0; i < NumCoords; i++) OrigLoc[i]= i;

          ML_az_dsort2(ytemp,NumCoords,OrigLoc);
          for (i = 0; i < NumCoords; i++) xtemp[i]= xvals[OrigLoc[i]];

          index = 0;

          while ( index < NumCoords ) {
             yfirst = ytemp[index];
             next   = index+1;
             while ( (next != NumCoords) && (ytemp[next] == yfirst))
             next++;
             ML_az_dsort2(&(xtemp[index]),next-index,&(OrigLoc[index]));
             for (i = index; i < next; i++) ztemp[i]= zvals[OrigLoc[i]];
             /* One final sort so that the ztemps are in order */
             subindex = index;
             while (subindex != next) {
                xfirst = xtemp[subindex]; subnext = subindex+1;
                while ( (subnext != next) && (xtemp[subnext] == xfirst)) subnext++;
                ML_az_dsort2(&(ztemp[subindex]),subnext-subindex,&(OrigLoc[subindex]));
                subindex = subnext;
             }
             index = next;
          }

         /* go through each vertical line and populate blockIndices so all   */
         /* dofs within a PDE within a vertical line correspond to one block.*/

         NumBlocks = 0;
         index = 0;

         while ( index < NumCoords ) {
            xfirst = xtemp[index];  yfirst = ytemp[index];
            next = index+1;
            while ( (next != NumCoords) && (xtemp[next] == xfirst) &&
                    (ytemp[next] == yfirst))
               next++;
            if (NumBlocks == 0) NumNodesPerVertLine = next-index;
            if (next-index != NumNodesPerVertLine) {
               printf("%d: Error code only works for constant block size now!!! A size of %d found instead of %d\n",mypid,next-index,NumNodesPerVertLine);
               exit(EXIT_FAILURE);
            }
            count = 0;
            for (j= index; j < next; j++) {
               VertLineId[OrigLoc[j]] = NumBlocks;
               LayerId[OrigLoc[j]] = count++;
            }
            NumBlocks++;
            index = next;
         }
         ML_free(ztemp);
         ML_free(ytemp);
         ML_free(xtemp);
         ML_free(OrigLoc);
       }

       /* check that everyone was assigned */

       for (i = 0; i < Nnodes;  i++) {
          if (VertLineId[i] == -1) {
             printf("Warning: did not assign %d to a vertical line?????\n",i);
          }
          if (LayerId[i] == -1) {
             printf("Warning: did not assign %d to a Layer?????\n",i);
          }
       }
       i = ML_gmax_int(NumNodesPerVertLine, comm);
       if (NumNodesPerVertLine == -1)  NumNodesPerVertLine = i;

       if (NumNodesPerVertLine != i)  {
          printf("%d: Different processors have different z direction line lengths? %d vs. %d\n",mypid,i,NumNodesPerVertLine);
          exit(EXIT_FAILURE);
       }
       return NumNodesPerVertLine;
}
