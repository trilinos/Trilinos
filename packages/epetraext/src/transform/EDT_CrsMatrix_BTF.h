
#ifndef ET_CRSGRAPH_BTF_H
#define ET_CRSGRAPH_BTF_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

class CrsGraph_BTF : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ~CrsGraph_BTF();

  CrsGraph_BTF()
  : NewRowMap_(0),
    NewDomainMap_(0),
    NewGraph_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

 private:

  Epetra_Map * NewRowMap_;
  Epetra_Map * NewDomainMap_;
  
  Epetra_CrsGraph * NewGraph_;

};

} //namespace EpetraExt

#endif //ET_CRSGRAPH_BTF_H
