
#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_BTF : public SameTypeDistTransform<Epetra_CrsGraph> {

 public:

  CrsGraph_BTF() {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

