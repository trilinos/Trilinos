
#ifndef NOX_EDT_CRSGRAPH_MAPCOLORING_H
#define NOX_EDT_CRSGRAPH_MAPCOLORING_H

#include <nox_Epetra_DistTransform.H>

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

  Epetra_MapColoring* operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //NOX_EDT_CRSGRAPH_MAPCOLORING_H
