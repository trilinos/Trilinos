
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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

// Trilinos Tutorial
// -----------------
// Use of Epetra_Operator.
// This code must be run with one process
//
// (output reported at the end of the file)
//
// Marzio Sala, SNL, 9214, 19-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

// ==================== //
// TriDiagonal Operator //
// -------------------- //

class TriDiagonalOperator : public Epetra_Operator 
{

public:
  
  // constructor
  TriDiagonalOperator( double diag_minus_one,
		       double diag,
		       double diag_plus_one,
		       Epetra_Map & Map) :
    diag_minus_one_(diag_minus_one),
    diag_(diag),
    diag_plus_one_(diag_plus_one),
    Map_( Map )
  {}

  // application of the tridiagonal operator
  int Apply( const Epetra_MultiVector & X,
	     Epetra_MultiVector & Y ) const
  {
    int Length = X.MyLength();
    
    // maybe some error checks on MultiVector Lenghts
    // for the future...
    
    for( int vec=0 ; vec<X.NumVectors() ; ++vec ) {
      
      // one-dimensional problems here
      if( Length == 1 ) {
	Y[vec][0] = diag_ * X[vec][0];
	break;
      }
      
      // more general case (Lenght >= 2)

      // first row
      Y[vec][0] = diag_ * X[vec][0] + diag_plus_one_ * X[vec][1];
      
      // intermediate rows
      for( int i=1 ; i<Length-1 ; ++i ) {
	Y[vec][i] = diag_ * X[vec][i] + diag_plus_one_ * X[vec][i+1]
	  + diag_minus_one_ * X[vec][i-1];
      }
      // final row
      Y[vec][Length-1] = diag_ * X[vec][Length-1]
	+ diag_minus_one_ * X[vec][Length-2];
    }
    
    return true;
  }

  // other function
  int SetUseTranspose( bool UseTranspose) 
  {}

  int ApplyInverse( const Epetra_MultiVector & X,
		    Epetra_MultiVector & Y ) const
  {
    return 0;
  }

  double NormInf() const
  {
    return( abs(diag_) + abs(diag_minus_one_) + abs(diag_plus_one_) );
  }

  char * Label () const
  {
    return "TriDiagonalOperator";
  }

  bool UseTranspose() const
  {
    return false;
  }

  bool HasNormInf () const
  {
    return true;
  }
  
  
  const Epetra_Comm & Comm() const
  {
    return( Map_.Comm() );
  }

  const Epetra_Map & OperatorDomainMap() const
  {
    return( Map_ );
  }
  
  const Epetra_Map & OperatorRangeMap() const
  {
    return( Map_ );
  }

  
private:

  Epetra_Map Map_;
  double diag_minus_one_;   // value in the sub-diagonal
  double diag_;             // value in the diagonal
  double diag_plus_one_;    // value in the super-diagonal
  
};

// =========== //
// main driver //
// ----------- //

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if( Comm.NumProc() != 1 ) {
    if( Comm.MyPID() == 0 ) {
      cerr << "This is mono-process example\n"
	   << "Please run with one processo only\n";
    }
    exit( EXIT_FAILURE );
  }
  
  int ierr;

  // global dimension of the problem, could be any positive number
  int NumGlobalElements( 5 );

  // linear decomposition (for simplicity, could be general)
  Epetra_Map Map(NumGlobalElements,0,Comm );

  // define two vectors based on Map
  Epetra_Vector x(Map);
  Epetra_Vector y(Map);
  x.PutScalar(1.0);
  
  // define a linear operator, as previously defined in class
  // TriDiagonalOperator

  TriDiagonalOperator TriDiagOp(-1.0,2.0,-1.0,Map);

  TriDiagOp.Apply(x,y);

  cout << x;
  cout << y;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

/*
output of this code

[s850675:epetra]> ex21.exe
Epetra::Vector
     MyPID           GID               Value
         0             0                       1
         0             1                       1
         0             2                       1
         0             3                       1
         0             4                       1
Epetra::Vector
     MyPID           GID               Value
         0             0                       1
         0             1                       0
         0             2                       0
         0             3                       0
         0             4                       1

*/	 
