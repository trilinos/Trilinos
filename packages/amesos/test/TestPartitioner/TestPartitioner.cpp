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

// Test class Amesos_Partitioner
//
#include "Amesos_config.h"
#include "Amesos.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos_Partitioner.h"
#include <vector>
#include "Amesos_Utils.h"

using namespace Trilinos_Util;

//==============================================================================
// General test utility. It can be run with any number of processors.
//
int Test(Epetra_Comm& Comm, int NumPoints,
	 char* PartitionerType, int OverlapLevel, 
	 int NumParts, bool UseGraph) 
{

  if (Comm.MyPID() == 0) {
    cout << "Test         = " << PartitionerType << endl;
    cout << "NumPoints    = " << NumPoints << endl;
    cout << "NumParts     = " << NumParts << endl;
    cout << "UseGraph     = " << UseGraph << endl;
    cout << "OverlapLevel = " << OverlapLevel << endl;
  }

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type","linear");

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // build the Amesos_Partitioner object

  Teuchos::ParameterList List;
  List.set("overlap level", OverlapLevel);
  List.set("partitioner", PartitionerType);
  List.set("local parts", NumParts);
  
  Amesos_Partitioner Partitioner(&(CrsA->Graph()),List);
  // FIXME

  AMESOS_CHK_ERR(Partitioner.Compute());

  // check the decomposition

  int NumMyOverlappingRows = Partitioner.NumMyOverlappingRows();
  vector<int> tmp;
  tmp.resize(NumMyOverlappingRows);
  for (int i = 0 ; i < NumMyOverlappingRows ; ++i)
    tmp[i] = -1;
  
  int count = 0;
  for (int i = 0 ; i < Partitioner.NumLocalParts() ; ++i) {
    int RowsInPart = Partitioner.NumRowsInPart(i);
    for (int j = 0 ; j < RowsInPart ; ++j) {
      int k = Partitioner(i,j,true);
      tmp[k] = 1;
    }
    count += RowsInPart;
  }

  // verify that all rows in overlapped graph 
  // are included in at least one part
  for (int i = 0 ; i < NumMyOverlappingRows ; ++i) {
    if (tmp[i] == -1) {
      AMESOS_CHK_ERR(-10);
    }
  }

  // check that we do not exceed maximum rrows
  if (count < Partitioner.NumMyOverlappingRows()) {
    AMESOS_CHK_ERR(-10);
  }

  // check that all rows are counted at least one
  if (OverlapLevel == 0) {
    if (count != NumMyOverlappingRows) {
      AMESOS_CHK_ERR(-10);
    }
  }
    
  // FIXME: I don't work with singletons
  return(0);

}
//==============================================================================
int Test1D(Epetra_Comm& Comm) 
{

  if (Comm.MyPID() == 0)
    cout << "Test1D" << endl;

  // size of the global matrix. 
  const int NumPoints = 256;
  
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type","linear");

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // build the Amesos_Partitioner object

  Teuchos::ParameterList List;
  List.set("overlap level", 1);
  List.set("partitioner", "1D");
  List.set("local parts", 4);
  
  Amesos_Partitioner Partitioner(&(CrsA->Graph()),List);

  AMESOS_CHK_ERR(Partitioner.Compute());

  // check the decomposition
  // verify that all rows in overlapped graph 
  // are included in at least one part

  int NumMyOverlappingRows = Partitioner.NumMyOverlappingRows();
  vector<int> tmp;
  tmp.resize(NumMyOverlappingRows);
  for (int i = 0 ; i < NumMyOverlappingRows ; ++i)
    tmp[i] = -1;
  
  bool ok = true;
  int count = 0;
  int NumParts = Partitioner.NumLocalParts();
  for (int i = 0 ; i < NumParts ; ++i) {
    int RowsInPart = Partitioner.NumRowsInPart(i);
    for (int j = 0 ; j < RowsInPart ; ++j) {
      int k = Partitioner(i,j,true);
      tmp[k] = 1;
    }
    if ((RowsInPart != 32) && (RowsInPart != 48)) 
      ok = false;
  }

  for (int i = 0 ; i < NumMyOverlappingRows ; ++i) {
    if (tmp[i] == -1) {
      AMESOS_CHK_ERR(-10);
    }
  }

  if (count > Partitioner.NumMyOverlappingRows()) {
    AMESOS_CHK_ERR(-10);
  }

  if (!ok) {
    AMESOS_CHK_ERR(-10);
  }
  else
    return(0);

}

//==============================================================================
int Test2D(Epetra_Comm& Comm)
{

  if (Comm.MyPID() == 0)
    cout << "Test2D" << endl;

  // size of the global matrix. 
  const int NumPoints = 256;
  
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type","box");

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // build the Amesos_Partitioner object

  Teuchos::ParameterList List;
  List.set("overlap level", 0);
  List.set("partitioner", "2D");
  List.set("local parts", 4);
  
  Amesos_Partitioner Partitioner(&(CrsA->Graph()),List);

  AMESOS_CHK_ERR(Partitioner.Compute());

  // check the decomposition
  bool ok = true;
  int count = 0;
  int NumParts = Partitioner.NumLocalParts();
  for (int i = 0 ; i < NumParts ; ++i) {
    int RowsInPart = Partitioner.NumRowsInPart(i);
    count += RowsInPart;
    if (RowsInPart != 16) 
      ok = false;
  }

  if (count != A->NumMyRows()) {
    AMESOS_CHK_ERR(-10);
  }

  if (!ok) {
    AMESOS_CHK_ERR(-10);
  }
  else
    return(0);

}

//==============================================================================
int TestMETIS(Epetra_Comm& Comm) 
{

  if (Comm.MyPID() == 0)
    cout << "TestMETIS" << endl;
 
  // size of the global matrix. 
  const int NumPoints = 4096;
  
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type","box");

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // build the Amesos_Partitioner object

  Teuchos::ParameterList List;
  List.set("overlap level", 0);
  List.set("partitioner", "METIS");
  List.set("local parts", 4);
  
  Amesos_Partitioner Partitioner(&(CrsA->Graph()),List);

  AMESOS_CHK_ERR(Partitioner.Compute());

  // check the decomposition
  bool ok = true;
  int count = 0;
  int NumParts = Partitioner.NumLocalParts();
  for (int i = 0 ; i < NumParts ; ++i) {
    int RowsInPart = Partitioner.NumRowsInPart(i);
    count += RowsInPart;
    // should be always 256, but maybe some random numbers can change it.
    if ((RowsInPart < 230) || (RowsInPart > 280)) 
      ok = false;
  }

  if (count != A->NumMyRows()) {
    AMESOS_CHK_ERR(-10);
  }

  if (!ok) {
    AMESOS_CHK_ERR(-10);
  }
  else
    return(0);

  return(0);


}

//==============================================================================
int TestGreedy(Epetra_Comm& Comm) 
{

  if (Comm.MyPID() == 0)
    cout << "TestGreddy" << endl;

  return(0);
}

//==============================================================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // general test, that can be run with any number of processes

  // 1D partition
  for (int overlap = 0 ; overlap < 5 ; ++overlap) {
    for (int NumParts = 1 ; NumParts < 32 ; NumParts += 4) {
      AMESOS_CHK_ERR(Test(Comm,4096, "1D", overlap, NumParts, true));
    }
  }

  // METIS partition
  for (int overlap = 0 ; overlap < 5 ; ++overlap) {
    for (int NumParts = 1 ; NumParts < 32 ; NumParts += 4) {
      AMESOS_CHK_ERR(Test(Comm,4096, "greedy", overlap, NumParts, true));
    }
  }

#define HAVE_AMESOS_METIS
#ifdef HAVE_AMESOS_METIS
  // METIS partition
  for (int overlap = 0 ; overlap < 5 ; ++overlap) {
    for (int NumParts = 1 ; NumParts < 32 ; NumParts += 4) {
      Test(Comm,4096, "METIS", overlap, NumParts, true);
    }
  }
#endif
  
  // more specialized tests, that has to be run with 4 processes
  if (Comm.NumProc() == 4) {

    // test 1D decomposition. 
    AMESOS_CHK_ERR(Test1D(Comm));

    // test 2D decomposition. This requires
    // four processors
    AMESOS_CHK_ERR(Test2D(Comm));

#ifdef HAVE_AMESOS_METIS
    // test METIS decomposition (if --enable-amesos_metis has
    // been enabled)
    AMESOS_CHK_ERR(TestMETIS(Comm));
#endif

    // test greedy decomposition
    AMESOS_CHK_ERR(TestGreedy(Comm));

  }
#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  // return
  exit(EXIT_SUCCESS);

}
