
#ifndef EpetraExt_CRSGRAPH_TRANSPOSE_H
#define EpetraExt_CRSGRAPH_TRANSPOSE_H

#include <EpetraExt_Transform.h>

class Epetra_CrsGraph;

namespace EpetraExt {

class CrsGraph_Transpose : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 bool ignoreNonLocalCols_;

 public:

  ~CrsGraph_Transpose();

  CrsGraph_Transpose( bool IgnoreNonLocalCols = false )
  : ignoreNonLocalCols_(IgnoreNonLocalCols)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );
};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_TRANSPOSE_H
