
#ifndef EDT_CRSGRAPH_ZOLTANORDER_H
#define EDT_CRSGRAPH_ZOLTANORDER_H

#ifdef ZOLTAN_ORDER

#include <Epetra_Transform.h>

class Zoltan_LoadBalance;

class Epetra_Map;
class Epetra_CrsGraph;

namespace EpetraExt {

class CrsGraph_ZoltanOrder : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  Epetra_Map * NewRowMap_;

 public:

  ~CrsGraph_ZoltanOrder();

  CrsGraph_ZoltanOrder()
  : NewRowMap_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //ZOLTAN_ORDER

#endif //EDT_CRSGRAPH_ZOLTANORDER_H
