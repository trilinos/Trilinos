// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Amesos_ConfigDefs.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <string>
//  #include "Trilinos_Util_ReadTriples2Epetra.h"
#include "Trilinos_Util.h"
//  #include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Amesos_Time.h"

#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SLUS
#include "Epetra_SLU.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef TEST_AZTEC
#include "AztecOO.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#if defined(HAVE_AMESOS_PARDISO) || defined(HAVE_AMESOS_PARDISO_MKL)
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_CSS_MKL
#include "Amesos_CssMKL.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif
#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)
#include "Amesos_Mumps.h"
#endif

#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h" 

//
//  Amesos_TestSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers and computes the error 
//  and residual.  
//
//  It reads the matrix in on a single processor and can pass that
//  matrix to the solver or it can convert that matrix to a 
//  distributed matrix and pass the distributed matrix to the solver.
//
//  Amesos_TestSolver can test either A x = b or A^T x = b.
//  This can be a bit confusing because sparse direct solvers 
//  use compressed column storage - the transpose of Trilinos'
//  sparse row storage.
//
//  Matrices:
//    readA - Serial.  As read from the file.
//    transposeA - Serial.  The transpose of readA.
//    serialA - if (transpose) then transposeA else readA 
//    distributedA - readA distributed to all processes
//    passA - if ( distributed ) then distributedA else serialA
//
//

int Amesos_TestSolver( Epetra_Comm &Comm, char *matrix_file,
                       SparseSolverType SparseSolver,
                       bool transpose,
                       int special, AMESOS_MatrixType matrix_type ) {


  Epetra_Map * readMap;
  Epetra_Vector * readx;
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
  Epetra_CrsMatrix *serialA;
   
  std::string FileName = matrix_file ;
  int FN_Size = FileName.size() ; 
  std::string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  std::string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );
  bool NonContiguousMap = false; 
  {
    Teuchos::RCP< Teuchos::Time > ioTimer = Teuchos::TimeMonitor::getNewCounter ("TestSolver::Matrix Read");
    Teuchos::TimeMonitor LocalTimer (*ioTimer);
    if ( LastFiveBytes == ".triU" ) {
      // Call routine to read in unsymmetric Triplet matrix
      NonContiguousMap = true; 
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, false, Comm, readMap, serialA, readx,
                                                        readb, readxexact, NonContiguousMap ) );
    } else {
      if ( LastFiveBytes == ".triS" ) { 
        NonContiguousMap = true; 
        // Call routine to read in symmetric Triplet matrix
        EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, true, Comm, readMap, serialA, readx,
                                                          readb, readxexact, NonContiguousMap ) );
      } else {
        if (  LastFourBytes == ".mtx" ) { 
          EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( matrix_file, Comm, readMap,
                                                                 serialA, readx, readb, readxexact) );
        } else {
          // Call routine to read in HB problem
          Trilinos_Util_ReadHb2Epetra( matrix_file, Comm, readMap, serialA, readx,
                                                     readb, readxexact);
        }
      }
    }
  }

  Epetra_CrsMatrix transposeA(Copy, *readMap, 0);

  if ( transpose ) {
    assert( CrsMatrixTranspose( serialA, &transposeA ) == 0 );
    serialA = &transposeA;
  }

  Epetra_RowMatrix * passA = 0; 
  Epetra_Vector * passx = 0; 
  Epetra_Vector * passb = 0;
  Epetra_Vector * passxexact = 0;
  Epetra_Vector * passresid = 0;
  Epetra_Vector * passtmp = 0;

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);
  Epetra_Map* map_;

  if( NonContiguousMap ) {
    //
    //  map gives us NumMyElements and MyFirstElement;
    //
    int NumGlobalElements =  readMap->NumGlobalElements();
    int NumMyElements = map.NumMyElements();
    int MyFirstElement = map.MinMyGID();
    std::vector<int> MapMap_( NumGlobalElements );
    readMap->MyGlobalElements( &MapMap_[0] ) ;
    Comm.Broadcast( &MapMap_[0], NumGlobalElements, 0 ) ; 
    map_ = new Epetra_Map( NumGlobalElements, NumMyElements, &MapMap_[MyFirstElement], 0, Comm);
  } else {
    map_ = new Epetra_Map( map ) ; 
  }


  Epetra_CrsMatrix A(Copy, *map_, 0);


  const Epetra_Map &OriginalMap = serialA->RowMatrixRowMap() ; 
  assert( OriginalMap.SameAs(*readMap) ); 
  Epetra_Export exporter(OriginalMap, *map_);
  Epetra_Export exporter2(OriginalMap, *map_);
  Epetra_Export MatrixExporter(OriginalMap, *map_);
  Epetra_CrsMatrix AwithDiag(Copy, *map_, 0);

  Epetra_Vector x(*map_);
  Epetra_Vector b(*map_);
  Epetra_Vector xexact(*map_);
  Epetra_Vector resid(*map_);
  Epetra_Vector readresid(*readMap);
  Epetra_Vector tmp(*map_);
  Epetra_Vector readtmp(*readMap);

  //  Epetra_Vector xcomp(*map_);      // X as computed by the solver
  bool distribute_matrix = ( matrix_type == AMESOS_Distributed ) ; 
  if ( distribute_matrix ) { 
    Teuchos::RCP< Teuchos::Time > distTimer = Teuchos::TimeMonitor::getNewCounter ("TestSolver::Matrix Distribute");
    Teuchos::TimeMonitor LocalTimer (*distTimer);
    // Create Exporter to distribute read-in matrix and vectors
    //
    //  Initialize x, b and xexact to the values read in from the file
    //
    x.Export(*readx, exporter, Add);
    b.Export(*readb, exporter, Add);
    xexact.Export(*readxexact, exporter, Add);
    Comm.Barrier();
    
    A.Export(*serialA, exporter, Add);
    //assert(A.FillComplete()==0);
    if(A.FillComplete()!=0) {
      exit(0);
    }
    delete serialA;
    
    Comm.Barrier();

    passA = &A; 

    passx = &x; 
    passb = &b;
    passxexact = &xexact;
    passresid = &resid;
    passtmp = &tmp;
  } else { 
    passA = serialA; 
    passx = readx; 
    passb = readb;
    passxexact = readxexact;
    passresid = &readresid;
    passtmp = &readtmp;
  }

  Epetra_MultiVector CopyB( *passb ) ;


  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA,
                                 (Epetra_MultiVector *) passx,
                                 (Epetra_MultiVector *) passb );
  

  for ( int i = 0; i < 1+special ; i++ ) { 
    Epetra_Time TotalTime( Comm ) ; 
    
    if ( false ) { 
      //  TEST_UMFPACK is never set by configure
#ifdef HAVE_AMESOS_SUPERLUDIST
    } else if ( SparseSolver == SUPERLUDIST ) {
        Teuchos::ParameterList ParamList ;
        ParamList.set( "MaxProcs", -3 );
        Amesos_Superludist A_Superludist( Problem );

  //ParamList.set( "Redistribute", true );
  //ParamList.set( "AddZeroToDiag", true );
  Teuchos::ParameterList& SuperludistParams = ParamList.sublist("Superludist") ;
  ParamList.set( "MaxProcs", -3 );

        EPETRA_CHK_ERR( A_Superludist.SetParameters( ParamList ) );
        EPETRA_CHK_ERR( A_Superludist.SetUseTranspose( transpose ) );
        EPETRA_CHK_ERR( A_Superludist.SymbolicFactorization(  ) );
        EPETRA_CHK_ERR( A_Superludist.NumericFactorization(  ) );
        EPETRA_CHK_ERR( A_Superludist.Solve(  ) );
#endif
#ifdef HAVE_AMESOS_DSCPACK
    } else if ( SparseSolver == DSCPACK ) {
      
      Teuchos::ParameterList ParamList ;
      ParamList.set( "MaxProcs", -3 );

      Amesos_Dscpack A_dscpack( Problem ) ; 
      EPETRA_CHK_ERR( A_dscpack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_dscpack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_dscpack.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_SCALAPACK
    } else if ( SparseSolver == SCALAPACK ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Scalapack A_scalapack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_scalapack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_scalapack.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_scalapack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_scalapack.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_scalapack.Solve(  ) ); 

#endif
#ifdef HAVE_AMESOS_TAUCS
    } else if ( SparseSolver == TAUCS ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Taucs A_taucs( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_taucs.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_taucs.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_taucs.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_taucs.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_taucs.Solve(  ) ); 

#endif
#if defined(HAVE_AMESOS_PARDISO) || defined(HAVE_AMESOS_PARDISO_MKL)
    } else if ( SparseSolver == PARDISO ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Pardiso A_pardiso( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_pardiso.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_pardiso.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_pardiso.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_pardiso.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_pardiso.Solve(  ) ); 

#endif
#ifdef HAVE_AMESOS_CSS_MKL
    } else if ( SparseSolver == CSS ) {

      Teuchos::ParameterList ParamList;
      Amesos_CssMKL A_css( Problem ); 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_css.SetParameters( ParamList ) );
      EPETRA_CHK_ERR( A_css.SetUseTranspose( transpose ) );
      {
        Teuchos::RCP< Teuchos::Time > symbolTimer = Teuchos::TimeMonitor::getNewCounter ("TestSolver::Symbolic(CSS)");
        Teuchos::TimeMonitor LocalTimer (*symbolTimer);
        EPETRA_CHK_ERR( A_css.SymbolicFactorization(  ) );
      }
      {
        Teuchos::RCP< Teuchos::Time > numericTimer = Teuchos::TimeMonitor::getNewCounter("TestSolver::Numeric (CSS)");
        Teuchos::TimeMonitor LocalTimer (*numericTimer);
        EPETRA_CHK_ERR( A_css.NumericFactorization(  ) );
      }
      {
        Teuchos::RCP< Teuchos::Time > solveTimer = Teuchos::TimeMonitor::getNewCounter  ("TestSolver::Solve   (CSS)");
        Teuchos::TimeMonitor LocalTimer (*solveTimer);
        EPETRA_CHK_ERR( A_css.Solve(  ) );
      }
      bool printTime = true;
      if (printTime) {
        A_css.PrintTiming();
        TimeMonitor::summarize();
      }
#endif
#ifdef HAVE_AMESOS_PARAKLETE
    } else if ( SparseSolver == PARAKLETE ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Paraklete A_paraklete( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_paraklete.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_paraklete.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_paraklete.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_paraklete.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_paraklete.Solve(  ) ); 

#endif
#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)
    } else if ( SparseSolver == MUMPS ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Mumps A_mumps( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_mumps.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_mumps.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_mumps.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_mumps.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_mumps.Solve(  ) ); 

#endif
#ifdef HAVE_AMESOS_SUPERLU
    } else if ( SparseSolver == SUPERLU ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Superlu A_superlu( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_superlu.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_superlu.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_superlu.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_superlu.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_superlu.Solve(  ) ); 

#endif
#ifdef HAVE_AMESOS_LAPACK
    } else if ( SparseSolver == LAPACK ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Lapack A_lapack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_lapack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_lapack.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_lapack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_lapack.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_lapack.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_UMFPACK
    } else if ( SparseSolver == UMFPACK ) {

      Teuchos::ParameterList ParamList ;
      Amesos_Umfpack A_umfpack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( A_umfpack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( A_umfpack.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( A_umfpack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( A_umfpack.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( A_umfpack.Solve(  ) ); 
#endif
#ifdef HAVE_AMESOS_KLU
    } else if ( SparseSolver == KLU ) {


      using namespace Teuchos;

      Amesos_Time AT; 
      int setupTimePtr = -1, symTimePtr = -1, numTimePtr = -1, refacTimePtr = -1, solveTimePtr = -1;
      AT.CreateTimer(Comm, 2);
      AT.ResetTimer(0);

      Teuchos::ParameterList ParamList ;
      // ParamList.set("OutputLevel",2);
      Amesos_Klu A_klu( Problem ); 
      ParamList.set( "MaxProcs", -3 );
      ParamList.set( "TrustMe", false );
      // ParamList.set( "Refactorize", true );
      EPETRA_CHK_ERR( A_klu.SetParameters( ParamList ) ) ; 
      EPETRA_CHK_ERR( A_klu.SetUseTranspose( transpose ) ); 
      setupTimePtr = AT.AddTime("Setup", setupTimePtr, 0);
      EPETRA_CHK_ERR( A_klu.SymbolicFactorization(  ) ); 
      symTimePtr = AT.AddTime("Symbolic", symTimePtr, 0);
      EPETRA_CHK_ERR( A_klu.NumericFactorization(  ) ); 
      numTimePtr = AT.AddTime("Numeric", numTimePtr, 0);
      EPETRA_CHK_ERR( A_klu.NumericFactorization(  ) ); 
      refacTimePtr = AT.AddTime("Refactor", refacTimePtr, 0);
      // for ( int i=0; i<100000 ; i++ ) 
      EPETRA_CHK_ERR( A_klu.Solve(  ) ); 
      solveTimePtr = AT.AddTime("Solve", solveTimePtr, 0);

      double SetupTime = AT.GetTime(setupTimePtr);
      double SymbolicTime = AT.GetTime(symTimePtr);
      double NumericTime = AT.GetTime(numTimePtr);
      double RefactorTime = AT.GetTime(refacTimePtr);
      double SolveTime = AT.GetTime(solveTimePtr);

      std::cout << __FILE__ << "::"  << __LINE__ << " SetupTime = " << SetupTime << std::endl ; 
      std::cout << __FILE__ << "::"  << __LINE__ << " SymbolicTime = " << SymbolicTime - SetupTime << std::endl ; 
      std::cout << __FILE__ << "::"  << __LINE__ << " NumericTime = " << NumericTime - SymbolicTime<< std::endl ; 
      std::cout << __FILE__ << "::"  << __LINE__ << " RefactorTime = " << RefactorTime - NumericTime << std::endl ; 
      std::cout << __FILE__ << "::"  << __LINE__ << " SolveTime = " << SolveTime - RefactorTime << std::endl ; 

#endif
    } else { 
      SparseDirectTimingVars::log_file << "Solver not implemented yet" << std::endl ;
      std::cerr << "\n\n####################  Requested solver not available on this platform ##################### ATS\n" << std::endl ;
      std::cout << " SparseSolver = " << SparseSolver << std::endl ; 
      std::cerr << " SparseSolver = " << SparseSolver << std::endl ; 
    }
    
    SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 
  }  // end for (int i=0; i<special; i++ ) 

  //
  //  Compute the error = norm(xcomp - xexact )
  //
  double error;
  passresid->Update(1.0, *passx, -1.0, *passxexact, 0.0);

  passresid->Norm2(&error);
  SparseDirectTimingVars::SS_Result.Set_Error(error) ;

  //  passxexact->Norm2(&error ) ; 
  //  passx->Norm2(&error ) ; 

  //
  //  Compute the residual = norm(Ax - b)
  //
  double residual ; 

  passA->Multiply( transpose, *passx, *passtmp);
  passresid->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
  //  passresid->Update(1.0, *passtmp, -1.0, CopyB, 0.0); 
  passresid->Norm2(&residual);

  SparseDirectTimingVars::SS_Result.Set_Residual(residual) ;
    
  double bnorm; 
  passb->Norm2( &bnorm ) ; 
  SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm) ;

  double xnorm; 
  passx->Norm2( &xnorm ) ; 
  SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm) ;

  delete readx;
  delete readb;
  if ( !distribute_matrix ) {
    delete serialA;
  }
  delete readxexact;
  delete readMap;
  delete map_;
  
  Comm.Barrier();

  return 0;
}
