
#ifndef EDT_MULTIVECTOR_VIEW_H
#define EDT_MULTIVECTOR_VIEW_H

#include <Epetra_Transform.h>

class Epetra_MultiVector;
class Epetra_BlockMap;

namespace EpetraExt {

class MultiVector_View : public ViewTransform<Epetra_MultiVector> {

  const Epetra_BlockMap & OrigMap_;
  const Epetra_BlockMap & NewMap_;

  const int NumVec_;

 public:

  ~MultiVector_View() {}

  MultiVector_View( const Epetra_BlockMap & orig_map,
                    const Epetra_BlockMap & new_map,
                    const int num_vec = -1 )
  : OrigMap_(orig_map),
    NewMap_(new_map),
    NumVec_(num_vec)
  {}

  NewTypePtr operator()( OriginalTypeRef original );

};

} //namespace EpetraExt

#endif //EDT_MULTIVECTOR_VIEW_H
