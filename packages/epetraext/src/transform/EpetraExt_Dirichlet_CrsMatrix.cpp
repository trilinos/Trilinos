#if 0 // dead source
// @HEADER
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
// @HEADER

#include <EpetraExt_Dirichlet_CrsMatrix.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_IntVector.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>

namespace EpetraExt {

bool
CrsMatrix_Dirichlet::
fwd()
{
  Epetra_CrsMatrix & Matrix = *origObj_;

  const Epetra_Map & RowMap = Matrix.RowMap();
  const Epetra_Map & ColMap = Matrix.ColMap();

  int NumMyElements = RowMap.NumMyElements();
  int NumMyColElements = ColMap.NumMyElements();

  if( symmetric_ && colSet_.empty() ) // get non-local column info
  {
    if( Matrix.IndicesAreGlobal() )
    {
      Epetra_Import Importer( ColMap, RowMap );
      Epetra_IntVector colLocations( ColMap );
      colLocations.Import( locations_, Importer, Insert );
      for( int i = 0; i < NumMyColElements; ++ i )
        if( colLocations[i] ) colSet_.insert(i);
    }
    else
    {
      for( int i = 0; i < NumMyElements; ++i )
        if( locations_[i] ) colSet_.insert(i);
    }
  }

  for( int i = 0; i < NumMyElements; ++i ) 
  {
    int * Indices;
    double * Vals;
    int NumIndices;
    if( locations_[i] ) //this is a Dirichlet BC location
    {
      Matrix.ExtractMyRowView( i, NumIndices, Vals, Indices );
      for( int j = 0; j < NumIndices; ++j )
      {
        if( Indices[j] == i ) Vals[i] = 1.0;
        else                  Vals[i] = 0.0;
      }
    }
    else if( symmetric_ )
    {
      Matrix.ExtractMyRowView( i, NumIndices, Vals, Indices );
      for( int j = 0; j < NumIndices; ++j )
        if( colSet_.count( Indices[j] ) ) Vals[j] = 0.0;
    }
  }

  return true;
}

bool
CrsMatrix_Dirichlet::
rvs()
{
  return true;
}

} //namespace EpetraExt

#endif // if 0
