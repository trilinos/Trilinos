
#ifndef EDT_CRSGRAPH_ZOLTANORDER_H
#define EDT_CRSGRAPH_ZOLTANORDER_H

#include <Epetra_Transform.h>

class Zoltan_LoadBalance;

class Epetra_CrsGraph;

namespace EpetraExt {
namespace Transform {

class CrsGraph_ZoltanOrder : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  Zoltan_LoadBalance * lb_;

 public:

  ~CrsGraph_ZoltanOrder() {}

  CrsGraph_ZoltanOrder( Zoltan_LoadBalance * lb = 0 )
  : lb_(lb)
  {}

  NewTypePtr operator()( OriginalTypeRef original );

};

} //namespace Transform
} //namespace EpetraExt

#endif //EDT_CRSGRAPH_ZOLTANORDER_H
