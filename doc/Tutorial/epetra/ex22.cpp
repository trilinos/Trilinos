
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
// Use of Epetra_Operator
// This code should be run with at least two processes
//
// (output reported at the end of the file)

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
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"

// auxiliary function to local an index in an integer vector

int find( int key, int vector[], int Length ) 
{
  for( int i=0 ; i<Length ; ++i ) {
    if( vector[i] == key ) return i;
  }
  return -1;
}

// ==================== //
// TriDiagonal Operator //
// -------------------- //

// NOTE: Will not work with IndexBase != 0
class TriDiagonalOperator : public Epetra_Operator 
{

public:
  
  // constructor
  TriDiagonalOperator( double diag_minus_one,
		       double diag,
		       double diag_plus_one,
		       const Epetra_Map & Map) :
    diag_minus_one_(diag_minus_one),
    diag_(diag),
    diag_plus_one_(diag_plus_one),
    Map_( Map )
  {
    // build the importer
    // Each local node will need the node+1 and node-1
    // (except for global node 0 and global node NumGlobalElemenets-1
    NumMyElements_ = Map_.NumMyElements();
    NumGlobalElements_ = Map_.NumGlobalElements();
    int MyGlobalElements[NumMyElements_];
    Map_.MyGlobalElements(MyGlobalElements);

    // count the nodes required from other processors
    // (this will be an upper bound of the # of required nodes
    // because I may count twice an external node)
    int count=0;
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int globalIndex = MyGlobalElements[i];
      // no -1 node for the first node of the grid
      if( globalIndex>0 )
	if( Map.LID(globalIndex-1) == -1 ) ++count;
      // now +1 node for the last node of the grid
      if( globalIndex<NumGlobalElements_-1 )
	if( Map.LID(globalIndex+1) == -1 ) ++count;
      ++count;
    }

    // now allocate space for local nodes and external nodes
    // (an external node is a node required for the matrix-vector
    // product, but owned by another process)
    int Length = count;
    int ListOfNodes[Length];

    count=0;
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int globalIndex = MyGlobalElements[i];
      // no -1 node for the first node of the grid
      if( globalIndex>0 ) {
	if( Map.LID(globalIndex-1) == -1 )
	  if( find( globalIndex-1, ListOfNodes, Length) == -1 ) {
	    ListOfNodes[count] = globalIndex-1;
	    ++count;
	  }
      }
      // now +1 node for the last node of the grid
      if( globalIndex<NumGlobalElements_-1 ) {
	if( Map.LID(globalIndex+1) == -1 ) {
	  if( find( globalIndex+1, ListOfNodes, Length) == -1 ) {
	    ListOfNodes[count] = globalIndex+1;
	    ++count;
	  }
	}
      }
      ListOfNodes[count] = globalIndex;
      ++count;
    }
    /*
    cout << "count = " << count << endl;
    for( int i=0 ; i<count ; i++ ) {
      cout << "ListOfNodes[" << i << "] = " << ListOfNodes[i] << endl;
    }
    */
    // create a Map defined using ListOfNodes

    ImportMap_ = new Epetra_Map(-1,count,ListOfNodes,0,Map_.Comm());

    Importer_ = new  Epetra_Import(*ImportMap_,Map_);

    return;
    
  }

  // application of the tridiagonal operator
  int Apply( const Epetra_MultiVector & X,
	     Epetra_MultiVector & Y ) const
  {

    cout << X;
    
    // maybe some error checks on MultiVector Lenghts
    // for the future...
    
    Epetra_MultiVector Xext((*ImportMap_),X.NumVectors());

    // this will contain local nodes and the required extenal nodes
    Xext.Import(X,*Importer_,Insert);
    
    for( int i=0 ; i<X.MyLength() ; ++i ) {
      
      int globalRow = Map_.GID(i);
      int iMinusOne = (*ImportMap_).LID(globalRow-1);
      int iPlusOne = (*ImportMap_).LID(globalRow+1);
      
      printf("%d %d %d\n", globalRow, iMinusOne, iPlusOne);

      for( int vec=0 ; vec<X.NumVectors() ; ++vec ) {
      
	Y[vec][i] = diag_ * X[vec][i];
	
	if( iMinusOne != -1 )
	  Y[vec][i] += diag_minus_one_ * Xext[vec][iMinusOne];
	
	if( iPlusOne != -1 )
	  Y[vec][i] += diag_plus_one_ * Xext[vec][iPlusOne];
	
      }
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
  int NumMyElements_;
  int NumGlobalElements_;
  double diag_minus_one_;   // value in the sub-diagonal
  double diag_;             // value in the diagonal
  double diag_plus_one_;    // value in the super-diagonal
  Epetra_Import *Importer_;
  Epetra_Map *ImportMap_;

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

  int ierr;

  // global dimension of the problem, could be any positive number
  int NumGlobalElements( 5 );

  // linear decomposition (for simplicity, could be general)
  Epetra_Map Map(NumGlobalElements,0,Comm );

  // define two vectors based on Map
  Epetra_Vector x(Map);
  Epetra_Vector y(Map);
  int NumMyElements = Map.NumMyElements();
  Epetra_IntSerialDenseVector MyGlobalElements(NumMyElements);
  Map.MyGlobalElements( MyGlobalElements.Values() );

  // x is a linear function, Laplace applied to it
  // should be zero except for the boundary nodes
  for( int i=0 ; i<NumMyElements ; ++i ) 
    x[i] = 1.0*MyGlobalElements[i];
  
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

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 2 ./ex22.exe
Epetra::Vector
     MyPID           GID               Value
         0             0                       0
         0             1                       1
         0             2                       2
Epetra::Vector
         1             3                       3
         1             4                       4
0 -1 1
1 0 3
2 1 2
Epetra::Vector
     MyPID           GID               Value
         0             0                       0
         0             1                       1
         0             2                       2
3 0 2
4 1 -1
Epetra::Vector
         1             3                       3
         1             4                       4
Epetra::Vector
     MyPID           GID               Value
         0             0                      -1
         0             1                       0
         0             2                       0
Epetra::Vector
         1             3                       0
         1             4                       5
*/
