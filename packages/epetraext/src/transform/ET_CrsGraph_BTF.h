
#ifndef ET_CRSGRAPH_BTF_H
#define ET_CRSGRAPH_BTF_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;

namespace EpetraExt {
namespace Transform {

class CrsGraph_BTF : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ~CrsGraph_BTF() {}

  CrsGraph_BTF() {}

  NewTypePtr operator()( OriginalTypeRef original );

};

} //namespace Transform
} //namespace EpetraExt

#endif //ET_CRSGRAPH_BTF_H
