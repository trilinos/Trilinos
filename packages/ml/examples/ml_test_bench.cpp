// Marian, this should be a possible example to start working with ML.
// I tried to put everything you may need in one file -- so that it is easier
// to compile and run.
//
// I believe that you can do the following to configure & run:
//
// 1- Download Trilinos
// 2- cd Trilinos
// 3- mkdir LINUX_SERIAL
// 4- cd LINUX_SERIAL
// 5- configure --enable-teuchos --enable-triutils
// 6- make
// 7- cd LINUX_SERIAL/packages/ml/examples
// 8- <modify this file as desired>
// 9- make
// 10- ./ml_test_bench.exe
//
// If you change and compile the ML library, make will recompile all the example.
// In this case, I suggest the following:
// 1- cd ../src
// 2- make
// 3- cd ../examples
// 4- touch *.exe; \rm ml_test_bench.exe
// 5- make
// This should compile this example only.
//
//
// Note: you can compile with mpi as well. In this case, you will probably need
//       the option --with-mpi-compilers only.
//
// Note2: now the code is creating a symmetric matrix (a 2D Laplacian). Other
//        matrices are available through the Trilinos gallery. You can consult
//        the Trilinos tutorial (Trilinos/doc/Tutorial/tutorial.pdf. If you need
//        other matrices, we can add them to the gallery, or code an other example.
//
// Note3: This is C++, but it is essentially a C code. I used C++ only to take
//        advantage of the Gallery (and Epetra).
//
// Note4: Most of the functions you may need are reported in this file. I erased
//        several lines, to make them more readable. This means that only a 
//        limited subset of the ML capabilities is available here, but it should
//        be ok for the time being.
//
// Note5: The code of the file is `TB_Gen_MultiLevelHierarchy_Adaptive()'
//

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_epetra_preconditioner.h"

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "ml_epetra_utils.h"
#include "ml_epetra_operator.h"

using namespace ML_Epetra;
using namespace Teuchos;
using namespace Trilinos_Util;

// function declarations

int TB_Gen_MultiLevelHierarchy_Adaptive(ML *ml, int start, 
					ML_Aggregate *ag);
int TB_AGG_Gen_Prolongator(ML *ml,int level, int clevel, void *data);
int TB_Aggregate_Coarsen( ML_Aggregate *ag, ML_Operator *Amatrix, 
                          ML_Operator **Pmatrix, ML_Comm *comm);
			  											
#include <iostream>

// MAIN DRIVER -- example of use of ML_Epetra::MultiLevelOperator
//
// from the command line, you may try something like that:
// $ mpirun -np 4 ./ml_example_epetra_operator.exe -problem_type=laplace_3d
//          -problem_size=1000
//
// For more options for Trilinos_Util::CrsMatrixGallery, consult the
// Trilinos 4.0 tutorial

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  Epetra_Time Time(Comm);
  
  CommandLineParser CLP(argc,argv);
  CrsMatrixGallery Gallery("", Comm);

  // default values for problem type and size
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "900" ); 

  // initialize MatrixGallery object with options specified in the shell
  Gallery.Set(CLP);
  
  // get pointer to the linear system matrix
  Epetra_CrsMatrix * A = Gallery.GetMatrix();

  // get a pointer to the map
  const Epetra_Map * Map = Gallery.GetMap();

  // get a pointer to the linear system problem
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();
  
  // Construct a solver object for this problem
  AztecOO solver(*Problem);
  
  // ================= MultiLevelOperator SECTION ========================

  int maxMgLevels = 2;       // max number of levels
  int nLevels = 0;           // will contain the # of actual levels  
  ML_Set_PrintLevel(10);     // print out all
  ML *ml_handle;             // handle for ML structures (not aggregates)
  ML_Aggregate *agg_object;  // handle for aggregates-related structured
  
  ML_Create(&ml_handle, maxMgLevels);

  // convert to ML matrix, put finest matrix into position 0 f the hierarchy
  EpetraMatrix2MLMatrix(ml_handle, 0, A);
  
  // create an Aggregate object; this will contain information
  // about the aggregation process for each level
  ML_Aggregate_Create(&agg_object);
  
  // up to here is as in a normal ML file.
  
  // =============== N E W   P A R T   F R O M   H E R E ================ //

  // create the hierarchy. Use increasing -- finest matrix is at level 0
    
  nLevels = TB_Gen_MultiLevelHierarchy_Adaptive(ml_handle, 0, 
				                agg_object);

  if( Comm.MyPID() == 0 ) cout << endl << "Number of levels = " << nLevels << endl << endl;
  
  // ============== E N D   O F   N E W   P A R T  ====-================= //
    
  // create an Epetra_Operator based on the previously created
  // hierarchy
  MultiLevelOperator MLPrec(ml_handle, Comm, *Map, *Map);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(&MLPrec);

  // for cg, use AZ_cg or AZ_cg_condnum
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 16);

  // solve with AztecOO
  solver.Iterate(500, 1e-5);

  // compute the real residual

  double residual, diff, res2;
  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
  (Gallery.GetExactSolution())->Norm2(&res2);
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff/res2 << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
  
}

// ========================================================================= //
// included functions. Some of those functions are copied from the ML        //
// library. If so, I used the prefix TB (test-bench) prefix instead of ML.   //
// ========================================================================= //

int TB_Gen_MultiLevelHierarchy_Adaptive(ML *ml, int start, 
					ML_Aggregate *ag)
{
  int    level;

  ML_Aggregate_Set_MaxLevels( ag, ml->ML_num_levels);
  ML_Aggregate_Set_StartLevel( ag, start );
   
  int size = (ml->Amat[0]).invec_leng;
  double * r = new double[size];
  double * z = new double[size];

  /* -------------------------------------------------------------------- */
  /* create multilevel hierarchy                                          */
  /* -------------------------------------------------------------------- */

  /* erased out part on different transpose/restriction for non-symmetric */
  /* problems, and timing. The idea is to get something clean (possibly)  */
  /* Also, I suppose ML_INCREASING only.                                  */
  /* Operator complexity is no longer computed.                           */
  /* I don't call Gen_MultiLevelHierarchy, and I put everything here for  */
  /* simplicity.                                                          */
  
  int next, flag, count=1;

  ml->ML_finest_level = start;
  level = start;
  next  = ML_AGG_Increment_Level(ml, level, (void *)ag);

  int zzz = 0; // trash for now only

  bool satisfied = false;
  
  do {

   // Now the only aggregation scheme here supported is Uncoupled.
   // (see comments in function `TB_Aggregate_Coarsen()'.)  
   // create the hierarchy
   while (next >= 0) 
   {
    if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel()) 
       printf("ML_Gen_MultiLevelHierarchy (level %d) : Gen Restriction and Prolongator \n",
              level );

    // this generates the smoothed restriction
    // I avoid the check after the construction...
    // The constructed operator is in Pmat = ml->Pmat+next;
    flag = TB_AGG_Gen_Prolongator(ml, level, next, (void *)ag);
    if (flag < 0) break;

    // now we have to generate the prolongator; it is simply the
    // transpose of the restriction in this case.
    ML_MultiLevel_Gen_Restriction(ml, level, next, (void *)ag);

    if ( ml->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel())
       printf("ML_Gen_MultiLevelHierarchy (level %d) : Gen RAP\n", level);
    ML_Gen_AmatrixRAP(ml, level, next);

    level = next;
    next  = ML_AGG_Increment_Level(ml, next, (void *)ag);

    count++;
   }

   int nLevels = count;
  
   // Here I set up the other components of the hierarchy. Usually, it is
   // done outside this function.
  
   // set up some smoothers. Here we suppose a symmetric problem
   int nits = 1;
   for(int level = 0 ; level < nLevels-1 ; ++level )
     ML_Gen_Smoother_GaussSeidel(ml, level, ML_BOTH, nits, ML_DEFAULT);

   // simple coarse solver. You may want to use Amesos to access
   // to a large variety of direct solvers, serial and parallel
   ML_Gen_Smoother_GaussSeidel(ml, nLevels-1, ML_BOTH, 
		              nits, ML_DEFAULT);
 
   // generate the solver. `0' is the finest level, and `nLevels' the coarsest.
   // Here we are using a V cycle.
   ML_Gen_Solver(ml, ML_MGV, 0, nLevels-1);

   // This applies the cycle to a vector. Vector elements are
   // still to be defined.
   ML_Solve_MGV(ml, r, z);

   // I suppose some checks will occur here on the quality of coarsening
   
   // Pick any idea now to define when we can stop
   if( ++zzz == 10 ) satisfied = true;
   
   // need to clean the hierarchy if bad choice were made
   if( satisfied == false ) {
     // RAY: how to do it ????
   }
   
  } while( satisfied == false );
   
  // free memory and return
  
  delete [] r;
  delete [] z;
           
  return(count);
}

/* ************************************************************************* */
/* generate smooth prolongator                                               */
/* ------------------------------------------------------------------------- */

int TB_AGG_Gen_Prolongator(ML *ml,int level, int clevel, void *data)
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
   prev_P_tentatives = (ML_Operator **) ag->P_tentative;

   /*
   widget.near_bdry = (char *) ML_allocate(sizeof(char)*Amat->outvec_leng);
   ML_AGG_Compute_Near_Bdry(Amat, widget.near_bdry);
   */

   Nfine    = Amat->outvec_leng;

   gNfine   = ML_Comm_GsumInt( ml->comm, Nfine);
   ML_Aggregate_Set_CurrentLevel( ag, level );

   // tried to simplify things here. Many of the ML options have been
   // erased -- I don't believe they are important at this point
   
   Pmatrix = ML_Operator_Create(ml->comm);
   Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&Pmatrix,ml->comm);

   // smoothing the prolongator. This should be the classical approach
   
   if ( ag->smoothP_damping_factor != 0.0 )
   {
     if (ml->symmetrize_matrix == ML_TRUE) {
       t2 = ML_Operator_Create(Amat->comm);
       ML_Operator_Transpose_byrow(Amat,t2);
       t3 = ML_Operator_Create(Amat->comm);
       ML_Operator_Add(Amat,t2,t3,ML_CSR_MATRIX,1.);
       max_eigen = t3->lambda_max;
     }

     if ((max_eigen < -666.) && (max_eigen > -667)) {

       switch( ag->spectral_radius_scheme ) {

       case 1:  /* compute it using CG */
	 
         kdata = ML_Krylov_Create( ml->comm );
         ML_Krylov_Set_PrintFreq( kdata, 0 );
         ML_Krylov_Set_ComputeEigenvalues( kdata );

	 if (ml->symmetrize_matrix ==ML_TRUE) ML_Krylov_Set_Amatrix(kdata, t3);
         else ML_Krylov_Set_Amatrix(kdata, Amat);

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

       case 2:
         // This would be Anasazi
	 
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
   }
   
   if ( ag->smoothP_damping_factor != 0.0 ) {
     widget.drop_tol = ag->drop_tol_for_smoothing;
     if (ml->symmetrize_matrix == ML_TRUE) widget.Amat   = t3;
     else widget.Amat   = &(ml->Amat[level]);
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
     
     if ((ml->symmetrize_matrix == ML_TRUE) ||
         (ag->use_transpose ==ML_TRUE) ) {
       if (t3 != NULL) ML_Operator_Destroy(&t3);
       if (t2 != NULL) ML_Operator_Destroy(&t2);
     }
     if (ag->keep_P_tentative == ML_NO)  ML_Operator_Destroy(&Pmatrix);
     else {
       if (prev_P_tentatives == NULL) {
	 ag->P_tentative = ML_Operator_ArrayCreate(ag->max_levels);
         prev_P_tentatives = (ML_Operator **) ag->P_tentative;
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

   return 0;
}

/* ************************************************************************* */
/* Coarsening routine                                                        */
/* ------------------------------------------------------------------------- */

int TB_Aggregate_Coarsen( ML_Aggregate *ag, ML_Operator *Amatrix, 
                          ML_Operator **Pmatrix, ML_Comm *comm)
{
   int Ncoarse;
   int mypid;

   mypid = comm->ML_mypid;

   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Coarsen : wrong object. \n");
      exit(-1);
   }

   if (mypid == 0 && ag->print_flag < ML_Get_PrintLevel()) 
      printf("ML_Aggregate_Coarsen (level %d) begins\n", ag->cur_level);

/* #### moved this somewhere else ?? */
   Amatrix->num_PDEs = ag->num_PDE_eqns;
   Amatrix->num_rigid = ag->nullspace_dim;

   // one can change the aggregation scheme as required. Simply copy
   // `ML_Aggregate_CoarsenUncoupled' here, change its name,
   // and change the following line accordingly.
   Ncoarse = ML_Aggregate_CoarsenUncoupled(ag,Amatrix,Pmatrix,comm);

   return Ncoarse;
}


#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */

