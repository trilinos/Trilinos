
#ifndef EDT_CRSGRAPH_OVERLAP_H
#define EDT_CRSGRAPH_OVERLAP_H

#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_Overlap : public SameTypeDistTransform<Epetra_CrsGraph> {

  const int levelOverlap_;

 public:

  ~CrsGraph_Overlap() {}

  CrsGraph_Overlap( const int overlap )
  : levelOverlap_(overlap)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_OVERLAP_H
