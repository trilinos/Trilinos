
#ifndef EDT_CRSMATRIX_VIEW_H
#define EDT_CRSMATRIX_VIEW_H

#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace Epetra_Transform {

class CrsMatrix_View : public SameTypeDistTransform<Epetra_CrsMatrix> {

  const Epetra_CrsGraph & OrigGraph_;
  const Epetra_CrsGraph & NewGraph_;

 public:

  ~CrsMatrix_View() {}

  CrsMatrix_View( const Epetra_CrsGraph & orig_graph,
                  const Epetra_CrsGraph & new_graph )
  : OrigGraph_(orig_graph),
    NewGraph_(new_graph)
  { /*Should test graphs for requirements*/ }

  std::auto_ptr<Epetra_CrsMatrix> operator()( const Epetra_CrsMatrix & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSMATRIX_VIEW_H
