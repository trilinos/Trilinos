
#ifndef EDT_CRSGRAPH_MAPCOLORINGINDEX_H
#define EDT_CRSGRAPH_MAPCOLORINGINDEX_H

#include <Epetra_DistTransform.h>

#include <vector>

class Epetra_CrsGraph;
class Epetra_MapColoring;
class Epetra_IntVector;

namespace Epetra_Transform {

class CrsGraph_MapColoringIndex : public DistTransform< Epetra_CrsGraph,vector<Epetra_IntVector> > {

 const Epetra_MapColoring & ColorMap_;

 public:

  ~CrsGraph_MapColoringIndex() {}

  CrsGraph_MapColoringIndex( const Epetra_MapColoring & ColorMap )
  : ColorMap_( ColorMap )
  {}

  std::auto_ptr< vector<Epetra_IntVector> > operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_MAPCOLORINGINDEX_H
