// @HEADER
// ***********************************************************************
// 
//                Amesos: Interface to Direct Solver Libraries
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

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Amesos_MC64.h"

class Ifpack_ReorderOperator 
{
 public:
  Ifpack_ReorderOperator(const int size, int* perm, double* scale) :
    size_(size),
    perm_(perm),
    scale_(scale)
  { }

  int Apply(double* x, double* y)
  {
    for (int i = 0 ; i < size_ ; ++i)
      y[i] = scale_[i] * x[perm_[i]];
  }
 private:
  const int size_;
  int* perm_;
  double* scale_;
};

// Creates the following matrix:
//
// | 3.00  5.00           |
// | 2.00  3.00      0.00 |
// | 3.00       4.00 6.00 |
// |            1.00 2.00 |
//
// as reported in the HSL MC64 version 1.3.0 user's guide. This simple driver
// calls MC64 using JOB=5, then prints on screen the CPERM and DW vectors. The
// reordered and scaled matrix looks like:
//
// | 1.00 1.00           |
// | 0.90 1.00      0.00 |
// |      1.00 1.00 1.00 |
// |           0.75 1.00 |
//
// \warning The interface between Amesos and MC64 works for serial
// computations only.
//
// \author Marzio Sala, ETHZ.
//
// \date Last updated on 05-Feb-06.

int main(int argc, char *argv[])
{
  // initialize MPI and Epetra communicator
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  if (Comm.NumProc() != 1)
  {
    std::cerr << "This example can be run with one processor only!" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  int NumMyRows = 4;
  Epetra_Map Map(NumMyRows, 0, Comm);
  Epetra_CrsMatrix A(Copy, Map, 0);

  int col;
  double val;

  col = 0; val = 3.0;
  A.InsertGlobalValues(0, 1, &val, &col);

  col = 1; val = 5.0;
  A.InsertGlobalValues(0, 1, &val, &col);

  col = 0; val = 2.0;
  A.InsertGlobalValues(1, 1, &val, &col);

  col = 1; val = 3.0;
  A.InsertGlobalValues(1, 1, &val, &col);

  col = 3; val = 0.0;
  A.InsertGlobalValues(1, 1, &val, &col);

  col = 0; val = 3.0;
  A.InsertGlobalValues(2, 1, &val, &col);

  col = 2; val = 4.0;
  A.InsertGlobalValues(2, 1, &val, &col);

  col = 3; val = 6.0;
  A.InsertGlobalValues(2, 1, &val, &col);

  col = 2; val = 1.0;
  A.InsertGlobalValues(3, 1, &val, &col);

  col = 3; val = 2.0;
  A.InsertGlobalValues(3, 1, &val, &col);

  A.FillComplete();

  // Creates an instance of the MC64 reordering and scaling interface
  // and computes the (column) reordering and scaling, using JOB=5
  Amesos_MC64 MC64(A, 5);

  // checks the return value
  std::cout << "INFO(1) = " << MC64.GetINFO(1) << std::endl;

  // Gets the pointer to reordering (CPERM) and scaling (DW). Both
  // vectors are allocated and free'd within Amesos_MC64.
  int*    CPERM = MC64.GetCPERM();
  double* DW    = MC64.GetDW();

  for (int i = 0 ; i < A.NumMyRows() ; ++i)
    std::cout << "CPERM[" << i << "] = " << CPERM[i] << std::endl;

  for (int i = 0 ; i < A.NumMyRows() * 2 ; ++i)
    std::cout << "DW[" << i << "] = " << DW[i] << std::endl;


  Ifpack_ReorderOperator RowPerm(4, MC64.GetRowPerm(), MC64.GetRowScaling());
  Ifpack_ReorderOperator ColPerm(4, MC64.GetColPerm(), MC64.GetColScaling());

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}
