
#include <Epetra_Transform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_BTF : public SameTypeTransform<Epetra_CrsGraph> {

 public:

  CrsGraph_BTF() {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

