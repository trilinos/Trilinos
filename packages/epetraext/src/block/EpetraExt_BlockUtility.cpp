//@HEADER
// ************************************************************************
// 
//               EpetraExt: Extended Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

namespace EpetraExt {

using std::vector;

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector<int> & RowStencil,
        int RowIndex,
        const Epetra_Comm & GlobalComm ) 
{
  const Epetra_BlockMap & BaseMap = BaseGraph.RowMap();
  int MaxGID = BaseMap.MaxAllGID();
  int BaseIndex = BaseMap.IndexBase();

  //Setup block offset
  int Offset = 1;
  while( Offset < MaxGID ) Offset *= 10;

  //Get Base Global IDs
  int Size = BaseMap.NumMyElements();
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  for( int i = 0; i < Size; ++i )
    GIDs[i] += RowIndex * Offset;

  int GlobalSize;
  GlobalComm.SumAll( &Size, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, Size, &GIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int> Indices(MaxIndices);
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  int StencilSize = RowStencil.size();
  for( int i = 0; i < Size; ++i )
  {
    int BaseRow = BaseMap.GID(i);
    int GlobalRow = GlobalMap.GID(i);

    BaseGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
    for( int j = 0; j < StencilSize; ++j )
    {
      int ColOffset = (RowIndex+RowStencil[j]) * Offset;
      if( j > 0 ) ColOffset -= (RowIndex+RowStencil[j-1]) * Offset;

      for( int k = 0; k < NumIndices; ++k )
        Indices[k] += ColOffset;

      GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices[0] );
    }
  }
  GlobalGraph->TransformToLocal();
  
  return GlobalGraph;
}

} //namespace EpetraExt
