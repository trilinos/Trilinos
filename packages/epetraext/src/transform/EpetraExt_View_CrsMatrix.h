
#ifndef EpetraExt_CRSMATRIX_VIEW_H
#define EpetraExt_CRSMATRIX_VIEW_H

#include <EpetraExt_Transform.h>

class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace EpetraExt {

class CrsMatrix_View : public ViewTransform<Epetra_CrsMatrix> {

  const Epetra_CrsGraph & OrigGraph_;
  const Epetra_CrsGraph & NewGraph_;

 public:

  ~CrsMatrix_View();

  CrsMatrix_View( const Epetra_CrsGraph & orig_graph,
                  const Epetra_CrsGraph & new_graph )
  : OrigGraph_(orig_graph),
    NewGraph_(new_graph)
  { /*Should test graphs for requirements*/ }

  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_CRSMATRIX_VIEW_H
