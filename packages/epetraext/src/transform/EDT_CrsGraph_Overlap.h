
#ifndef EDT_CRSGRAPH_OVERLAP_H
#define EDT_CRSGRAPH_OVERLAP_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;

namespace EpetraExt {

class CrsGraph_Overlap : public SameTypeTransform<Epetra_CrsGraph> {

  const int levelOverlap_;

 public:

  ~CrsGraph_Overlap() {}

  CrsGraph_Overlap( const int overlap )
  : levelOverlap_(overlap)
  {}

  NewTypePtr operator()( OriginalTypeRef original );

  bool fwd() { return true; }
  bool rvs() { return false; }

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_OVERLAP_H
