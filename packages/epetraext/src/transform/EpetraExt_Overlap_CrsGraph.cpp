
#include <EpetraExt_Overlap_CrsGraph.h>

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
