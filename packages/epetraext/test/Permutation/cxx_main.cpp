//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
//@HEADER

//Permutation Test routine
#include <Epetra_ConfigDefs.h>
#include "EpetraExt_Version.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "EpetraExt_Permutation.h"
#include "../epetra_test_err.h"

int check_rowpermute_crsmatrix_local_diagonal(Epetra_Comm& Comm, bool verbose);
int check_rowpermute_crsmatrix_global_diagonal(Epetra_Comm& Comm, bool verbose);
int check_rowpermute_crsgraph_local_diagonal(Epetra_Comm& Comm, bool verbose);
int check_colpermute_crsgraph(Epetra_Comm& Comm, bool verbose);
int check_colpermute_crsmatrix(Epetra_Comm& Comm, bool verbose);
int check_rowpermute_multivector_local(Epetra_Comm& Comm, bool verbose);

int main(int argc, char *argv[]) {

  int returnierr=0;

  bool verbose = false;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else
  Epetra_SerialComm Comm;
#endif

  // Check if we should print results to standard out
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  }

  //Make sure the value of verbose is consistent across processors.
  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  if (!verbose) {
    Comm.SetTracebackMode(0); // This should shut down error traceback reporting
  }

  if (verbose && Comm.MyPID()==0)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  EPETRA_CHK_ERR( check_rowpermute_crsmatrix_local_diagonal( Comm, verbose ) );

  EPETRA_CHK_ERR( check_rowpermute_crsmatrix_global_diagonal( Comm, verbose) );

  EPETRA_CHK_ERR( check_rowpermute_crsgraph_local_diagonal( Comm, verbose) );

  EPETRA_CHK_ERR( check_colpermute_crsgraph( Comm, verbose) );

  EPETRA_CHK_ERR( check_colpermute_crsmatrix( Comm, verbose) );

  EPETRA_CHK_ERR( check_rowpermute_multivector_local( Comm, verbose) );


#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

//------------------------------------------------------------------------------
int check_rowpermute_crsmatrix_local_diagonal(Epetra_Comm& Comm,
					   bool verbose)
{
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) {
    cerr << "================check_rowpermute_crsmatrix_local_diagonal=========="
	 <<endl;
  }

  int NumMyElements = 5;
  int NumGlobalElements = NumMyElements*NumProc;
 
  Epetra_Map Map(NumGlobalElements, NumMyElements, 0, Comm);

  int* p = new int[NumMyElements];
  int firstGlobalRow = MyPID*NumMyElements;

  //Set up a permutation that will reverse the order of all LOCAL rows. (i.e.,
  //this test won't cause any inter-processor data movement.)

  if (verbose) {
    cout << "Permutation P:"<<endl;
  }

  int i;

  for(i=0; i<NumMyElements; ++i) {
    p[i] = firstGlobalRow+NumMyElements-1-i;
    if (verbose1) {
      cout << "p["<<firstGlobalRow+i<<"]: "<<p[i]<<endl;
    }
  }

  Epetra_CrsMatrix A(Copy, Map, 1);

  int col;
  double val;

  //set up a diagonal matrix A. It's diagonal because that's the easiest
  //to fill and to examine output before and after permutation...

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    val = 1.0*row;
    col = row;

    A.InsertGlobalValues(row, 1, &val, &col);
  }
 
  A.FillComplete();

  if (verbose1) {
    cout << "********** matrix A: **************"<<endl;
    cout << A << endl;
  }

  EpetraExt::Permutation<Epetra_CrsMatrix> P(Copy, Map, p);

  Epetra_CrsMatrix& B = P(A);

  if (verbose1) {
    cout <<"************ permuted matrix B: ***************"<<endl;
    cout << B << endl;
  }

  return(0);
}

//------------------------------------------------------------------------------
int check_rowpermute_crsgraph_local_diagonal(Epetra_Comm& Comm,
					  bool verbose)
{
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) {
    cerr << "================check_rowpermute_crsgraph_local_diagonal=========="
	 <<endl;
  }

  int NumMyElements = 5;
  int NumGlobalElements = NumMyElements*NumProc;
 
  Epetra_Map Map(NumGlobalElements, NumMyElements, 0, Comm);

  int* p = new int[NumMyElements];
  int firstGlobalRow = MyPID*NumMyElements;

  //Set up a permutation that will reverse the order of all LOCAL rows. (i.e.,
  //this test won't cause any inter-processor data movement.)

  if (verbose) {
    cout << "Permutation P:"<<endl;
  }

  int i;

  for(i=0; i<NumMyElements; ++i) {
    p[i] = firstGlobalRow+NumMyElements-1-i;
    if (verbose1) {
      cout << "p["<<firstGlobalRow+i<<"]: "<<p[i]<<endl;
    }
  }

  Epetra_CrsGraph Agrph(Copy, Map, 1);

  int col;

  //set up a diagonal graph. It's diagonal because that's the easiest
  //to fill and to examine output before and after permutation...

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    col = row;

    Agrph.InsertGlobalIndices(row, 1, &col);
  }
 
  Agrph.FillComplete();

  if (verbose1) {
    cout << "*************** graph Agrph: ********************"<<endl;
    cout << Agrph << endl;
  }

  EpetraExt::Permutation<Epetra_CrsGraph> P(Copy, Map, p);

  Epetra_CrsGraph& Bgrph = P(Agrph);

  if (verbose1) {
    cout <<"************* permuted graph Bgrph: ****************"<<endl;
    cout << Bgrph << endl;
  }

  return(0);
}

//------------------------------------------------------------------------------
int check_colpermute_crsgraph(Epetra_Comm& Comm,
			      bool verbose)
{
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) {
    cerr << "================check_colpermute_crsgraph=========="
	 <<endl;
  }

  int NumMyElements = 5;
  int NumGlobalElements = NumMyElements*NumProc;
 
  Epetra_Map Map(NumGlobalElements, NumMyElements, 0, Comm);

  int* p = new int[NumMyElements];
  int firstGlobalRow = MyPID*NumMyElements;

  if (verbose) {
    cout << "Permutation P:"<<endl;
  }

  int i;

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    p[i] = NumGlobalElements - row - 1;
    if (verbose1) {
      cout << "p["<<firstGlobalRow+i<<"]: "<<p[i]<<endl;
    }
  }

  Epetra_CrsGraph Agrph(Copy, Map, 1);

  int col;

  //set up a tri-diagonal graph.

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    col = NumGlobalElements - row - 1;

    Agrph.InsertGlobalIndices(row, 1, &col);

    if (col > 0) {
      int colm1 = col-1;
      Agrph.InsertGlobalIndices(row, 1, &colm1);
    }

    if (col < NumGlobalElements-1) {
      int colp1 = col+1;
      Agrph.InsertGlobalIndices(row, 1, &colp1);
    }
  }
 
  Agrph.FillComplete();

  if (verbose1) {
    cout << "*************** graph Agrph: ********************"<<endl;
    cout << Agrph << endl;
  }

  EpetraExt::Permutation<Epetra_CrsGraph> P(Copy, Map, p);

  bool column_permutation = true;
  Epetra_CrsGraph& Bgrph = P(Agrph, column_permutation);

  if (verbose1) {
    cout <<"************* column-permuted graph Bgrph: ****************"<<endl;
    cout << Bgrph << endl;
  }

  delete [] p;

  return(0);
}

//-------------------------------------------------------------------------------
int check_rowpermute_crsmatrix_global_diagonal(Epetra_Comm& Comm,
			bool verbose)
{
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) {
    cerr << "================check_rowpermute_crsmatrix_global_diagonal=========="
	 <<endl;
  }

  int NumMyElements = 5;
  int NumGlobalElements = NumMyElements*NumProc;
 
  Epetra_Map Map(NumGlobalElements, NumMyElements, 0, Comm);

  int* p = new int[NumMyElements];
  int firstGlobalRow = MyPID*NumMyElements;

  //Now set up a permutation that will GLOBALLY reverse the order of all rows.
  //(i.e., if there are multiple processors, there will be inter-processor
  //data movement as rows are migrated.)

  int i;

  Epetra_CrsMatrix A(Copy, Map, 1);

  int col;
  double val;

  //set up a diagonal matrix A. It's diagonal because that's the easiest
  //to fill and to examine output before and after permutation...

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    val = 1.0*row;
    col = row;

    A.InsertGlobalValues(row, 1, &val, &col);
  }
 
  A.FillComplete();

  if (verbose1) {
    cout << "******************* matrix A: ****************************"<<endl;
    cout << A << endl;
  }

  if (verbose) {
    cout << "Permutation P:"<<endl;
  }

  for(i=0; i<NumMyElements; ++i) {
    int globalrow = NumGlobalElements-(firstGlobalRow+i)-1;
    p[i] = globalrow;
    if (verbose1) {
      cout << "p["<<firstGlobalRow+i<<"]: "<<p[i]<<endl;
    }
  }

  EpetraExt::Permutation<Epetra_CrsMatrix> Pglobal(Copy, Map, p);

  Epetra_CrsMatrix& Bglobal = Pglobal(A);

  if (verbose1) {
    cout << "******************* permuted matrix Bglobal: *******************" <<endl;
    cout << Bglobal << endl;
  }

  return(0);
}

//-------------------------------------------------------------------------------
int check_colpermute_crsmatrix(Epetra_Comm& Comm,
			       bool verbose)
{
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) {
    cerr << "================check_colpermute_crsmatrix=========="
	 <<endl;
  }

  int NumMyElements = 5;
  int NumGlobalElements = NumMyElements*NumProc;
 
  Epetra_Map Map(NumGlobalElements, NumMyElements, 0, Comm);

  int* p = new int[NumMyElements];
  int firstGlobalRow = MyPID*NumMyElements;

  if (verbose) {
    cout << "Permutation P:"<<endl;
  }

  int i;

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    p[i] = NumGlobalElements - row - 1;
    if (verbose1) {
      cout << "p["<<firstGlobalRow+i<<"]: "<<p[i]<<endl;
    }
  }

  Epetra_CrsMatrix A(Copy, Map, 1);

  int col;
  double val;

  //set up a tri-diagonal graph.

  for(i=0; i<NumMyElements; ++i) {
    int row = firstGlobalRow+i;
    col = NumGlobalElements - row - 1;
    val = 1.0*col;

    A.InsertGlobalValues(row, 1, &val, &col);

    if (col > 0) {
      int colm1 = col-1;
      val = 1.0*colm1;
      A.InsertGlobalValues(row, 1, &val, &colm1);
    }

    if (col < NumGlobalElements-1) {
      int colp1 = col+1;
      val = 1.0*colp1;
      A.InsertGlobalValues(row, 1, &val, &colp1);
    }
  }
 
  A.FillComplete();

  if (verbose1) {
    cout << "*************** matrix A: ********************"<<endl;
    cout << A << endl;
  }

  EpetraExt::Permutation<Epetra_CrsMatrix> P(Copy, Map, p);

  bool column_permutation = true;
  Epetra_CrsMatrix& B = P(A, column_permutation);

  if (verbose1) {
    cout <<"************* column-permuted matrix B: ****************"<<endl;
    cout << B << endl;
  }

  delete [] p;

  return(0);
}

//------------------------------------------------------------------------------
int check_rowpermute_multivector_local(Epetra_Comm& Comm,
				       bool verbose)
{
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) {
    cerr << "================check_rowpermute_multivector_local=========="
	 <<endl;
  }

  int NumMyElements = 5;
  int NumGlobalElements = NumMyElements*NumProc;
 
  Epetra_Map Map(NumGlobalElements, NumMyElements, 0, Comm);

  int* p = new int[NumMyElements];
  int firstGlobalRow = MyPID*NumMyElements;

  //Set up a permutation that will reverse the order of all LOCAL rows. (i.e.,
  //this test won't cause any inter-processor data movement.)

  if (verbose) {
    cout << "Permutation P:"<<endl;
  }

  int i;

  for(i=0; i<NumMyElements; ++i) {
    p[i] = firstGlobalRow+NumMyElements-1-i;
    if (verbose1) {
      cout << "p["<<firstGlobalRow+i<<"]: "<<p[i]<<endl;
    }
  }

  Epetra_MultiVector v(Map, 3);

  double* v0 = v[0];
  double* v1 = v[1];
  double* v2 = v[2];

  for(i=0; i<NumMyElements; ++i) {
    v0[i] = 1.0*(firstGlobalRow+i) + 0.1;
    v1[i] = 1.0*(firstGlobalRow+i) + 0.2;
    v2[i] = 1.0*(firstGlobalRow+i) + 0.3;
  }
 
  if (verbose1) {
    cout << "*************** MultiVector v: ********************"<<endl;
    cout << v << endl;
  }

  EpetraExt::Permutation<Epetra_MultiVector> P(Copy, Map, p);

  Epetra_MultiVector& Pv = P(v);

  if (verbose1) {
    cout <<"************* permuted MultiVector Pv: ****************"<<endl;
    cout << Pv << endl;
  }

  return(0);
}

