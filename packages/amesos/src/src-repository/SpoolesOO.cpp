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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "SpoolesOO.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include <sstream>
#include <vector>

//
//  SPOOLES include file:
//
extern "C" {
#include "BridgeMPI.h"
}

//=============================================================================
SpoolesOO::SpoolesOO(Epetra_RowMatrix * A, 
		 Epetra_MultiVector * X,
		 Epetra_MultiVector * B) {
  //  AllocAzArrays();
  SetSpoolesDefaults();

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  SetUserMatrix(A);
  
  SetLHS(X);
  SetRHS(B);
  inConstructor_ = false;
}

//=============================================================================
SpoolesOO::SpoolesOO() {
  //  AllocAzArrays();
  SetSpoolesDefaults();
}

//=============================================================================
SpoolesOO::~SpoolesOO(void) {
  //  DeleteMemory();
  //  DeleteAzArrays();
}

//=============================================================================
int SpoolesOO::SetUserMatrix(Epetra_RowMatrix * UserMatrix) {

  if (UserMatrix == 0 && inConstructor_ == true) return(0);
  if (UserMatrix == 0) EPETRA_CHK_ERR(-1);

  UserMatrix_ = UserMatrix;

  return(0);
}

//=============================================================================
int SpoolesOO::SetLHS(Epetra_MultiVector * X) {

  if (X == 0 && inConstructor_ == true) return(0);
  if (X == 0) EPETRA_CHK_ERR(-1);
  X_ = X;
  X_->ExtractView(&x_, &x_LDA_);
  return(0);
}
//=============================================================================
int SpoolesOO::SetRHS(Epetra_MultiVector * B) {

  if (B == 0 && inConstructor_ == true) return(0);
  if (B == 0) EPETRA_CHK_ERR(-1);
  B_ = B;
  B_->ExtractView(&b_, &b_LDA_);

  return(0);
}
int SpoolesOO::SetSpoolesDefaults() {

 UserOperator_ = 0;
 UserMatrix_ = 0;
 // PrecOperator_ = 0;
 // PrecMatrix_ = 0;
 X_ = 0;
 B_ = 0;
 
 x_LDA_ = 0;
 x_ = 0;
 b_LDA_ = 0;
 b_ = 0;

 return(0);

}

//=============================================================================

int SpoolesOO::Solve() { 
  FILE *msgFile = fopen("msgFile", "w" ) ; 
  FILE *matFile = 0 ;

  Epetra_RowMatrix   *RowMatrixA = (GetUserMatrix()) ; 
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  int iam =  Comm.MyPID() ;
  bool verbose = false ; // Other option is (iam == 0 ) ; 
  verbose = ( iam == 0 ) ; 
  if ( verbose ) matFile = fopen("SPO.m", "w" ) ; 

  Epetra_CrsMatrix   *matA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  assert( matA != NULL ) ; // SuperludistOO shows how to convert a 
  // non-CrsMatrix into a CrsMatrix

  Epetra_MultiVector   *vecX = GetLHS() ; 
  Epetra_MultiVector   *vecB = GetRHS() ; 

  const Epetra_Map &matAmap = matA->RowMap() ; 
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPI_Comm MPIC = comm1.Comm() ;

  int IsLocal = ( matAmap.NumMyElements() == 
		  matAmap.NumGlobalElements() )?1:0;
  Comm.Broadcast( &IsLocal, 1, 0 ) ; 

  EPETRA_CHK_ERR( 1 - IsLocal  ); // SuperludistOO shows hows to 
  // deal with distributed matrices.  

  int nArows = matA->NumGlobalRows() ; 
  int nAcols = matA->NumGlobalCols() ; 
  int nproc =  Comm.NumProc() ;
  
  assert( vecX->NumVectors() == 1 ) ; 
  assert( vecB->NumVectors() == 1 ) ; 

  //  assert( nArows == vecX->MyLength() ) ; 
  //  assert( nAcols == vecB->MyLength() ) ;

  //  assert( matAmap.DistributedGlobal() == false ) ; 
  if ( iam == 0 ) { 
    assert( matAmap.NumMyElements() == matAmap.NumGlobalElements() ) ;
  } else { 
    assert( matAmap.NumMyElements() == 0 ) ;
  } 

  int numentries = matA->NumGlobalNonzeros();
  vector <int>rowindices( numentries ) ; 
  vector <int>colindices( numentries ) ; 
  vector <double>values( numentries ) ; 
  int NumEntries;
  double *RowValues;

  int *ColIndices;
  int numrows = matA->NumGlobalRows();
  int entrynum = 0 ; 

  //
  //  The following lines allow us time to attach the debugger
  //
  char hatever;
  //  if ( iam == 0 )  cin >> hatever ; 
  Comm.Barrier();


  InpMtx *mtxA; 
  if ( iam==0 ) {
    //
    // Create mtxA on proc 0 
    //
    for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
      assert( matA->ExtractMyRowView( MyRow, NumEntries, RowValues, ColIndices ) == 0 ) ;
      for ( int j = 0; j < NumEntries; j++ ) { 
	rowindices[entrynum] = MyRow ; 
	colindices[entrynum] = ColIndices[j] ; 
	values[entrynum] = RowValues[j] ; 
	entrynum++; 
      }
    }
    mtxA = InpMtx_new() ; 
    
    InpMtx_init( mtxA, 1, SPOOLES_REAL, 0, 0 ) ;   // I can't figure out 
    //  what the right bound for the number of vectors is - Ken 
    
    if ( GetTrans() ) { 
      InpMtx_inputRealTriples( mtxA, numentries, &colindices[0], 
			       &rowindices[0], &values[0] ) ; 
    } else { 
      InpMtx_inputRealTriples( mtxA, numentries, &rowindices[0], 
			       &colindices[0], &values[0] ) ; 
    }
  } else {
    mtxA = 0 ; 
  }

  //  return 1;    OK to here
  
  DenseMtx *mtxX = 0 ; 
  DenseMtx *mtxY = 0 ;
  double *bValues ;
  double *xValues ;
  if ( iam == 0 ) { 
    //
    //  Convert Epetra_Vector x and Epetra_Vector b arrays
    //
    int bLda, xLda ; 
    
    assert( vecB->ExtractView( &bValues, &bLda ) == 0 )  ; 
    assert( vecX->ExtractView( &xValues, &xLda ) == 0 ) ; 
  
  
    //
    //  SPOOLES matrices
    //
    mtxX = DenseMtx_new() ; 
    mtxY = DenseMtx_new() ; 
#ifdef OLDWAY
    DenseMtx_init( mtxX, SPOOLES_REAL, -1, -1, numrows, 1, 1, numrows );
    DenseMtx_init( mtxY, SPOOLES_REAL, -1, -1, numrows, 1, 1, numrows );
#else
    DenseMtx_init( mtxX, SPOOLES_REAL, 0, 0, numrows, 1, 1, numrows );
    DenseMtx_init( mtxY, SPOOLES_REAL, 0, 0, numrows, 1, 1, numrows );
#endif
    int nrow, nnrhs ; 
    DenseMtx_dimensions( mtxY, &nrow, &nnrhs ) ;
    assert( nrow = numrows ) ; 
    assert( nnrhs == 1 ) ; 
    
  
   if (verbose) InpMtx_writeForMatlab( mtxA, "mtxA", matFile ) ; 
    //
    //  This is a maximally inefficient way to create the right hand side
    //  matrix, but I could not find anything better in the documentation.
    //
    for ( int i = 0 ; i < numrows; i ++ ) 
      {
	DenseMtx_setRealEntry( mtxY, i, 0, bValues[i] );
      }
    if ( verbose ) DenseMtx_writeForMatlab( mtxY, "mtxY", matFile ) ; 
    DenseMtx_zero(mtxX) ;
  }
  
  
  int type = SPOOLES_REAL ;
  int symmetryflag = SPOOLES_NONSYMMETRIC ;
  // SPOOLES requires a message level and a message File 
  int msglvl = 0;        
  int rc;                // return code 
  
  /*
    --------------------------------
    create and setup a Bridge object
    --------------------------------
  */
  BridgeMPI     *bridgeMPI = BridgeMPI_new() ;
  BridgeMPI_setMPIparams(bridgeMPI, nproc, iam, MPIC ) ;
  BridgeMPI_setMatrixParams(bridgeMPI, numrows, type, symmetryflag) ;
  double tau = 100.0 ; 
  double droptol = 0.0 ; 
  int lookahead = 0 ; 
  PatchAndGoInfo *PatchAndGo = 0 ; 
  BridgeMPI_setFactorParams(bridgeMPI, FRONTMTX_DENSE_FRONTS, 
			    SPOOLES_PIVOTING, tau, droptol, lookahead, 
			    PatchAndGo ) ; 
  BridgeMPI_setMessageInfo(bridgeMPI, msglvl, msgFile) ;
  assert( type == SPOOLES_REAL ) ; 
  assert( msglvl >= 0 && msglvl <= 2 ) ; 
    
  rc = BridgeMPI_setup(bridgeMPI, mtxA) ;
  if ( rc != 1 ) {
    fprintf(stderr, "\n error return %d from BridgeMPI_setup()", rc) ;
    exit(-1) ;
  }
  Comm.Barrier();
  
  int nfront, nfind, nfent, nsolveops;
  double nfactorops;
  rc = BridgeMPI_factorStats(bridgeMPI, type, symmetryflag, &nfront,
			     &nfind, &nfent, &nsolveops, &nfactorops) ;
  if ( rc != 1 ) {
    fprintf(stderr,
	    "\n error return %d from BridgeMPI_factorStats()", rc) ;
    exit(-1) ;
  }
  
  /*
    --------------------------------
    setup the parallel factorization
    --------------------------------
  */
  rc = BridgeMPI_factorSetup(bridgeMPI, 0, 0.0) ;
  if (rc != 1 ) {
    std::stringstream Message ; 
    
    Message <<  " SPOOLES factorsetup failed with return code " << 
      rc ;
    string mess = Message.str() ; 
    throw( mess ) ; 
  }
  
  /*
    -----------------
    factor the matrix
    -----------------
  */
  int permuteflag = 1 ;
  int errorcode ;
  rc = BridgeMPI_factor(bridgeMPI, mtxA, permuteflag, &errorcode) ;
  assert( permuteflag == 1 ) ; 
  
  if (rc != 1 ) {
    std::stringstream Message ; 
    
    Message <<  " SPOOLES factorization failed with return code " << 
      rc << " and error code " << errorcode  ; 
    string mess = Message.str() ; 
    throw( mess ) ; 
  }
  /*
    ----------------
    solve the system
    ----------------
  */
  
  rc = BridgeMPI_solveSetup(bridgeMPI) ;
  if (rc != 1 ) {
    std::stringstream Message ; 
    
    Message <<  " SPOOLES BridgeMPI_solveSetup failed with return code" << 
      rc << endl ; 
    string mess = Message.str() ; 
    throw( mess ) ; 
  }

  assert( permuteflag == 1 ) ; 
  rc = BridgeMPI_solve(bridgeMPI, permuteflag, mtxX, mtxY) ;
  assert( permuteflag == 1 ) ; 
  if (rc != 1 ) {
    std::stringstream Message ; 
    
    Message <<  " SPOOLES BridgeMPI_solve failed with return code" << 
      rc << endl ; 
    string mess = Message.str() ; 
    throw( mess ) ; 
  }
  
  if ( verbose ) DenseMtx_writeForMatlab( mtxX, "mtxXX", matFile ) ; 
  if ( verbose ) fclose(matFile);

  //  Result->SolveTime().Time_First( ) ;
  //
  //  Here we copy the values of B back from mtxX into xValues
  //
  if ( iam == 0 ) { 
    for ( int i = 0 ; i < numrows; i ++ ) 
      {
	DenseMtx_realEntry( mtxX, i, 0, &xValues[i] );
      }
    //  DenseMtx_writeForMatlab( mtxX, "mtxX", matFile ) ; 
  }
  

  if ( iam == 0 ) {
    InpMtx_free(mtxA) ;
    DenseMtx_free(mtxX) ;
    DenseMtx_free(mtxY) ;
  }
  BridgeMPI_free(bridgeMPI) ;

  Comm.Barrier();

  return(0) ; 
}

