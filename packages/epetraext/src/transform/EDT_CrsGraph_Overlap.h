
#ifndef EDT_CRSGRAPH_OVERLAP_H
#define EDT_CRSGRAPH_OVERLAP_H

#include <Epetra_Transform.h>

class Epetra_BlockMap;
class Epetra_CrsGraph;

namespace EpetraExt {

class CrsGraph_Overlap : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  const int levelOverlap_;

  const bool squareLocalBlock_;

  Epetra_BlockMap * OverlapMap_;

 public:

  ~CrsGraph_Overlap();

  CrsGraph_Overlap( int overlap, bool squareLocalBlock = false )
  : levelOverlap_(overlap),
    squareLocalBlock_(squareLocalBlock),
    OverlapMap_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_OVERLAP_H
