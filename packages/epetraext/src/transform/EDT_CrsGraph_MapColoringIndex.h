
#ifndef EDT_CRSGRAPH_MAPCOLORINGINDEX_H
#define EDT_CRSGRAPH_MAPCOLORINGINDEX_H

#include <Epetra_Transform.h>

#include <vector>

using std::vector;

class Epetra_CrsGraph;
class Epetra_MapColoring;
class Epetra_IntVector;

namespace EpetraExt {

class CrsGraph_MapColoringIndex : public StructuralTransform< Epetra_CrsGraph,vector<Epetra_IntVector> > {

 const Epetra_MapColoring & ColorMap_;

 public:

  ~CrsGraph_MapColoringIndex() {}

  CrsGraph_MapColoringIndex( const Epetra_MapColoring & ColorMap )
  : ColorMap_( ColorMap )
  {}

  NewTypePtr operator()( OriginalTypeRef original );
};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_MAPCOLORINGINDEX_H
