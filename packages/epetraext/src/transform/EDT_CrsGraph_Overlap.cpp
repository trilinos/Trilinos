//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include <EDT_CrsGraph_Overlap.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

namespace EpetraExt {

CrsGraph_Overlap::
~CrsGraph_Overlap()
{
  if( newObj_ ) delete newObj_;

  if( OverlapMap_ ) delete OverlapMap_;
}

CrsGraph_Overlap::NewTypeRef
CrsGraph_Overlap::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int err;

  //check that this is a distributed graph and overlap level is not zero
  if( orig.DomainMap().DistributedGlobal() && levelOverlap_ )
  {
    Epetra_CrsGraph * OverlapGraph = new Epetra_CrsGraph( orig );
    OverlapMap_ = new Epetra_BlockMap( orig.RowMap() );

    Epetra_BlockMap * DomainMap = &(const_cast<Epetra_BlockMap&>(orig.DomainMap()));
    Epetra_BlockMap * RangeMap = &(const_cast<Epetra_BlockMap&>(orig.RangeMap()));

    for( int level = 0; level < levelOverlap_; ++level )
    {
      Epetra_BlockMap * OldRowMap = OverlapMap_;
      Epetra_CrsGraph * OldGraph = OverlapGraph;

      Epetra_Import & OverlapImporter = *(const_cast<Epetra_Import *>( OldGraph->Importer() ));
      OverlapMap_ = new Epetra_BlockMap( OverlapImporter.TargetMap() );

      //filter to local square block on last level if required
      if( squareLocalBlock_ && level==(levelOverlap_-1) )
        OverlapGraph = new Epetra_CrsGraph( Copy, *OverlapMap_, *OverlapMap_, 0 );
      else
        OverlapGraph = new Epetra_CrsGraph( Copy, *OverlapMap_, 0 );

      OverlapGraph->Import( *OldGraph, OverlapImporter, Insert );
      OverlapGraph->TransformToLocal( DomainMap, RangeMap );

      delete OldGraph;
      delete OldRowMap;
    }

    newObj_ = OverlapGraph;
  }
  else //just create a copy since this is not a InPlaceTransform
    newObj_ = new Epetra_CrsGraph( orig );

  return *newObj_;
}

} // namespace EpetraExt
