
#ifndef EDT_CRSGRAPH_TRANSPOSE_H
#define EDT_CRSGRAPH_TRANSPOSE_H

#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_Transpose : public SameTypeDistTransform<Epetra_CrsGraph> {

 bool ignoreNonLocalCols_;

 public:

  ~CrsGraph_Transpose() {}

  CrsGraph_Transpose( bool IgnoreNonLocalCols = false )
  : ignoreNonLocalCols_(IgnoreNonLocalCols)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_TRANSPOSE_H
