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

//
//  CreateTridi populates an empty EpetraCrsMatrix with a tridiagonal with 
//  -1 on the off-diagonals and 2 on the diagonal.  
//
//  CreateTridiPlus creates the same matrix as CreateTridi except that it adds
//  -1 in the two off diagonal corners.
//
//  This code was plaguerized from epetra/example/petra_power_method/cxx_main.cpp
//  presumably written by Mike Heroux.
//
//  Adapted by Ken Stanley - Aug 2003 
//
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

int CreateTridi(Epetra_CrsMatrix& A)
{

  Epetra_Map Map = A.RowMap();
  int NumMyElements = Map.NumMyElements();
  int NumGlobalElements = Map.NumGlobalElements();

  int * MyGlobalElements = new int[NumMyElements];
    Map.MyGlobalElements(MyGlobalElements);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[3];
  int *Indices = new int[3];
#ifndef NDEBUG
  int NumEntries;
#endif
  for (int i=0; i<NumMyElements; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 0;
	Indices[1] = 1;
	Values[0] = 2.0;
	Values[1] = -1.0;
#ifndef NDEBUG
	NumEntries = 2;
#endif
      }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
      {
	Indices[0] = NumGlobalElements-1;
	Indices[1] = NumGlobalElements-2;
	Values[0] = 2.0;
	Values[1] = -1.0;
#ifndef NDEBUG
	NumEntries = 2;
#endif
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i];
	Indices[2] = MyGlobalElements[i]+1;
	Values[0] = -1.0; 
	Values[1] = 2.0;
	Values[2] = -1.0;
#ifndef NDEBUG
	NumEntries = 3;
#endif
      }
    
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
     // Put in the diagonal entry
     //     assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, &MyGlobalElements[i])==0);
    }
  
  // Finish up
  assert(A.FillComplete()==0);


  delete[] MyGlobalElements;
  delete[] Values;
  delete[] Indices;
  return 0;
}

int CreateTridiPlus(Epetra_CrsMatrix& A)
{

  Epetra_Map Map = A.RowMap();
  int NumMyElements = Map.NumMyElements();
  int NumGlobalElements = Map.NumGlobalElements();

  int * MyGlobalElements = new int[NumMyElements];
    Map.MyGlobalElements(MyGlobalElements);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[3];
  int *Indices = new int[3];
#ifndef NDEBUG
  int NumEntries;
#endif
  for (int i=0; i<NumMyElements; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 0;
	Indices[1] = 1;
	Indices[2] = NumGlobalElements-1;
	Values[0] = 2.0;
	Values[1] = -1.0;
	Values[2] = -0.5;
#ifndef NDEBUG
	NumEntries = 3;
#endif
      }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
      {
	Indices[0] = NumGlobalElements-1;
	Indices[1] = NumGlobalElements-2;
	Indices[2] = 0;
	Values[0] = 2.0;
	Values[1] = -1.0;
	Values[2] = -0.5;
#ifndef NDEBUG
	NumEntries = 3;
#endif
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i];
	Indices[2] = MyGlobalElements[i]+1;
	Values[0] = -1.0; 
	Values[1] = 2.0;
	Values[2] = -1.0;
#ifndef NDEBUG
	NumEntries = 3;
#endif
      }
    
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
     // Put in the diagonal entry
     //     assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, &MyGlobalElements[i])==0);
    }
  
  // Finish up
  assert(A.FillComplete()==0);


  delete[] MyGlobalElements;
  delete[] Values;
  delete[] Indices;
  return 0;
}

