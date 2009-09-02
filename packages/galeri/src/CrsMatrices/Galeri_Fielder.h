// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_FIEDLER_H
#define GALERI_FIEDLER_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* Fielder(const Epetra_Map* Map)
{
  // this is actually a dense matrix, stored into Crs format
  int NumGlobalElements = Map->NumGlobalElements();
  int NumMyElements     = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, NumGlobalElements);

  vector<double> Values(NumGlobalElements);
  vector<int>    Indices(NumGlobalElements);

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = NumGlobalElements;
    int iGlobal = MyGlobalElements[i];
    
    for (int j = 0 ; j < NumGlobalElements ; ++j) 
    {
      Indices[j] = j;
      Values[j]  = (double)abs(iGlobal - j);
    }

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
  }

  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
