
#ifndef EDT_CRSGRAPH_ZOLTANORDER_H
#define EDT_CRSGRAPH_ZOLTANORDER_H

#include <Epetra_DistTransform.h>

class Zoltan_LoadBalance;

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_ZoltanOrder : public SameTypeDistTransform<Epetra_CrsGraph> {

  Zoltan_LoadBalance * lb_;

 public:

  ~CrsGraph_ZoltanOrder() {}

  CrsGraph_ZoltanOrder( Zoltan_LoadBalance * lb = 0 )
  : lb_(lb)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_ZOLTANORDER_H
