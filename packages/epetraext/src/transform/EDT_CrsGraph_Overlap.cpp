
#include <EDT_CrsGraph_Overlap.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

namespace EpetraExt {
namespace Transform {

NewTypePtr CrsGraph_Overlap::operator()( OriginalTypeRef original )
{
  int err;

  Epetra_CrsGraph * OverlapGraph = 0;

  if( original.DomainMap().DistributedGlobal() || !levelOverlap_ )
  {
    Epetra_Import * OverlapImporter =
                    const_cast<Epetra_Import *>( original.Importer() );
    Epetra_BlockMap * OverlapMap =
                    new Epetra_BlockMap( OverlapImporter->TargetMap() );
    Epetra_BlockMap * DomainMap =
                    const_cast<Epetra_BlockMap *>( &original.DomainMap() );

    OverlapGraph = new Epetra_CrsGraph( Copy, *OverlapMap, 0 );
    OverlapGraph->Import( original, *OverlapImporter, Insert );
    OverlapGraph->TransformToLocal( DomainMap, OverlapMap );

    for( int level = 1; level < levelOverlap_; ++level )
    {
      Epetra_BlockMap * OldMap = OverlapMap;
      Epetra_CrsGraph * OldGraph = OverlapGraph;

      OverlapImporter = const_cast<Epetra_Import *>( OldGraph->Importer() );
      OverlapMap = new Epetra_BlockMap( OverlapImporter->TargetMap() );

      OverlapGraph = new Epetra_CrsGraph( Copy, *OverlapMap, 0 );
      OverlapGraph->Import( *OldGraph, *OverlapImporter, Insert );
      OverlapGraph->TransformToLocal( DomainMap, OverlapMap );

      delete OldGraph;
      delete OldMap;
    }

  }
  
  return OverlapGraph;
}

} //namespace Transform
} //namespace EpetraExt
