//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef HAVE_NOX_EPETRAEXT

// Use local version of needed coloring file
#include <nox_EDT_CrsGraph_MapColoringIndex.H>

#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>

#include <map>

namespace Epetra_Transform {

vector<Epetra_IntVector>* CrsGraph_MapColoringIndex::operator()( const Epetra_CrsGraph & original )
{
  const Epetra_BlockMap & RowMap = original.RowMap();
  int nRows = RowMap.NumMyElements();

  int NumColors = ColorMap_.NumColors();
  int * ListOfColors = ColorMap_.ListOfColors();

  map<int,int> MapOfColors;
  for( int i = 0; i < NumColors; ++i ) MapOfColors[ ListOfColors[i] ] = i;

  //initial setup of stl vector of IntVectors for indexing
  vector<int> dummy( nRows, -1 );
  vector<Epetra_IntVector>* IndexVec(
	new vector<Epetra_IntVector>( NumColors, Epetra_IntVector( Copy, RowMap, &dummy[0] ) ) );

  int MaxNumIndices = original.MaxNumIndices();
  int NumIndices;
  vector<int> Indices( MaxNumIndices );

  for( int i = 0; i < nRows; ++i )
  {
    original.ExtractGlobalRowCopy( original.GRID(i), MaxNumIndices, NumIndices, &Indices[0] );

    for( int j = 0; j < NumIndices; ++j )
     (*IndexVec)[ MapOfColors[ColorMap_(Indices[j])] ][i] = Indices[j];
  }

  return IndexVec;
}

}
#endif
