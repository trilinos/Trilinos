
#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_Transpose : public SameTypeDistTransform<Epetra_CrsGraph> {

 public:

  ~CrsGraph_Transpose() {}

  CrsGraph_Transpose() {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

