
#ifndef ET_CRSGRAPH_BTF_H
#define ET_CRSGRAPH_BTF_H

#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_BTF : public SameTypeDistTransform<Epetra_CrsGraph> {

 public:

  ~CrsGraph_BTF() {}

  CrsGraph_BTF() {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //ET_CRSGRAPH_BTF_H
