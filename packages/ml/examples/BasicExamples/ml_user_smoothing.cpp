
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// The goal of this example is to show how a user may use her own smoothing
// function within the ML_Epetra::MultiLevelPreconditioner class.
// This example builds a 2D Poisson linear system and solves it via AztecOO,
// with ML as a preconditioner. ML uses a user-defined damped Jacobi smoother
// (defined in this file).
//
#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example also
// requires --enable-galeri (for the definition of the linear systems)
// and --enable-aztecoo (to solve the linear system)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

// epetra objects
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
// required to build the example matrix
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

extern int userSmoother(ML_Smoother *data, int x_length, double x[],
                   int rhs_length, double rhs[]);

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Creates the linear problem using the Galeri package.
  // Several matrix examples are supported; please refer to the
  // Galeri documentation for more details.
  // Most of the examples using the ML_Epetra::MultiLevelPreconditioner
  // class are based on Epetra_CrsMatrix. Example
  // `ml_EpetraVbr.cpp' shows how to define a Epetra_VbrMatrix.
  // `Laplace2D' is a symmetric matrix; an example of non-symmetric
  // matrices is `Recirc2D' (advection-diffusion in a box, with
  // recirculating flow). The grid has nx x ny nodes, divided into
  // mx x my subdomains, each assigned to a different processor.

  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 8;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);

  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    
  // Build a linear system with trivial solution, using a random vector
  // as starting solution.
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // As we wish to use AztecOO, we need to construct a solver object 
  // for this problem
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // Sets default parameters for classic smoothed aggregation. After this
  // call, MLList contains the default values for the ML parameters,
  // as required by typical smoothed aggregation for symmetric systems.
  // Other sets of parameters are available for non-symmetric systems
  // ("DD" and "DD-ML"), and for the Maxwell equations ("maxwell").
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // output level, 0 being silent and 10 verbose
  MLList.set("ML output", 10);
  // maximum number of levels
  MLList.set("max levels",5);
  // set finest level to 0
  MLList.set("increasing or decreasing","increasing");

  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");

  // smoother is symmetric Gauss-Seidel. Example file 
  // `ml/examples/TwoLevelDD/ml_2level_DD.cpp' shows how to use
  // AZTEC's preconditioners as smoothers
  //MLList.set("smoother: type","symmetric Gauss-Seidel");

  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");

  MLList.set("smoother: type", "user-defined");
  //The function pointer to the user-defined smoother.
  MLList.set("smoother: user-defined function",userSmoother);
  //(optional) The label for the user-defined smoother.
  MLList.set("smoother: user-defined name","myJacobi");
  //Use this smoother on the finest level, 0.
  MLList.set("smoother: type (level 0)", "user-defined");
  MLList.set("smoother: sweeps (level 0)", 1);
  //Use symmetric Gauss-Seidel on all subsequent levels.
  MLList.set("smoother: type (level 1)", "symmetric Gauss-Seidel");

  //user-defined smoothing for a 1-level method only
  //to use this, uncomment the next 3 lines and comment out the direct-solver
  //settings further below
  //MLList.set("coarse: type", "user-defined");
  //MLList.set("coarse: user-defined function",userSmoother);
  //MLList.set("coarse: user-defined name","myJacobi");


#ifdef HAVE_ML_AMESOS
  // solve with serial direct solver KLU
  MLList.set("coarse: type","Amesos-KLU");
#else
  // this is for testing purposes only, you should have 
  // a direct solver for the coarse problem (either Amesos, or the SuperLU/
  // SuperLU_DIST interface of ML)
  MLList.set("aggregation: type", "MIS");
  MLList.set("smoother: type","Jacobi");
  MLList.set("coarse: type","Jacobi");
#endif

  MLList.set("read XML", false); // skip XML file in this directory

  // Creates the preconditioning object. We suggest to use `new' and
  // `delete' because the destructor contains some calls to MPI (as
  // required by ML and possibly Amesos). This is an issue only if the
  // destructor is called **after** MPI_Finalize().
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // verify unused parameters on process 0 (put -1 to print on all
  // processes)
  MLPrec->PrintUnused(0);

  // =========================== end of ML part =============================
  
  // tell AztecOO to use the ML preconditioner, specify the solver 
  // and the output, then solve with 500 maximum iterations and 1e-12 
  // of tolerance (see AztecOO's user guide for more details)
  
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  // destroy the preconditioner
  delete MLPrec;
  
  // compute the real residual

  double residual;
  LHS.Norm2(&residual);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  // for testing purposes
  if (residual > 1e-5)
    exit(EXIT_FAILURE);

  delete A;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} //main

int userSmoother(ML_Smoother *data, int x_length, double x[],
                 int rhs_length, double rhs[])
{
   int i;
   double *ap, omega = .5; /* temp vector and damping factor */
   double *diag;
   ML_Operator *Amat;
   ML_Smoother *smoo;

   smoo    = (ML_Smoother *) data;
   Amat = (ML_Operator *) ML_Get_MySmootherData(smoo);
   ap = (double *) ML_allocate(Amat->outvec_leng * sizeof(double));
   ML_Operator_Apply(Amat, x_length, x, rhs_length, ap);
   ML_Operator_Get_Diag(Amat, x_length, &diag);
   
   for (i = 0; i < x_length; i++) x[i] = x[i] + omega*(rhs[i] - ap[i])/diag[i];

   ML_free(ap);

   return 0;
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
#if !defined(HAVE_ML_EPETRA)
  puts("--enable-epetra");
#endif
#if !defined(HAVE_ML_TEUCHOS)
  puts("--enable-teuchos");
#endif
#if !defined(HAVE_ML_AZTECOO)
  puts("--enable-aztecoo");
#endif
#if !defined(HAVE_ML_GALERI)
  puts("--enable-galeri");
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}

#endif
