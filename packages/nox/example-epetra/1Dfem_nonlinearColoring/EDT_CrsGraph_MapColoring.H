
#ifndef EDT_CRSGRAPH_MAPCOLORING_H
#define EDT_CRSGRAPH_MAPCOLORING_H

#include <Epetra_DistTransform.H>

class Epetra_CrsGraph;
class Epetra_MapColoring;

namespace Epetra_Transform {

class CrsGraph_MapColoring : public DistTransform<Epetra_CrsGraph,Epetra_MapColoring> {

  bool verbose_;

 public:

  ~CrsGraph_MapColoring() {}

  CrsGraph_MapColoring( bool verbose = false )
  : verbose_(verbose)
  {}

  std::auto_ptr<Epetra_MapColoring> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_MAPCOLORING_H
