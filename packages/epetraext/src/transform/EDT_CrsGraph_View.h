
#ifndef EDT_CRSGRAPH_VIEW_H
#define EDT_CRSGRAPH_VIEW_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;
class Epetra_BlockMap;

namespace EpetraExt {
namespace Transform {

class CrsGraph_View : public ViewTransform<Epetra_CrsGraph> {

  const Epetra_BlockMap * NewRowMap_;
  const Epetra_BlockMap * NewColMap_;

 public:

  ~CrsGraph_View() {}

  CrsGraph_View( const Epetra_BlockMap * new_row_map,
                 const Epetra_BlockMap * new_col_map = 0 )
  : NewRowMap_(new_row_map),
    NewColMap_(new_col_map)
  {}

  NewTypePtr operator()( OriginalTypeRef original );

};

} //namespace Transform
} //namespace EpetraExt

#endif //EDT_CRSGRAPH_VIEW_H
