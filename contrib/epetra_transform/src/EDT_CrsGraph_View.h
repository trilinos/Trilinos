
#ifndef EDT_CRSGRAPH_VIEW_H
#define EDT_CRSGRAPH_VIEW_H

#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;
class Epetra_BlockMap;

namespace Epetra_Transform {

class CrsGraph_View : public SameTypeDistTransform<Epetra_CrsGraph> {

  const Epetra_BlockMap * NewRowMap_;
  const Epetra_BlockMap * NewColMap_;

 public:

  ~CrsGraph_View() {}

  CrsGraph_View( const Epetra_BlockMap * new_row_map,
                 const Epetra_BlockMap * new_col_map = 0 )
  : NewRowMap_(new_row_map),
    NewColMap_(new_col_map)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_VIEW_H
