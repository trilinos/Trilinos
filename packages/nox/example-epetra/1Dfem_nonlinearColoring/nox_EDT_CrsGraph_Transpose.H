
#ifndef NOX_EDT_CRSGRAPH_TRANSPOSE_H
#define NOX_EDT_CRSGRAPH_TRANSPOSE_H

#include <nox_Epetra_DistTransform.H>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_Transpose : public SameTypeDistTransform<Epetra_CrsGraph> {

 bool ignoreNonLocalCols_;

 public:

  ~CrsGraph_Transpose() {}

  CrsGraph_Transpose( bool IgnoreNonLocalCols = false )
  : ignoreNonLocalCols_(IgnoreNonLocalCols)
  {}

  Epetra_CrsGraph* operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //NOX_EDT_CRSGRAPH_TRANSPOSE_H
