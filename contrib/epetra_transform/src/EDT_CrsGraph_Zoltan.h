
#ifndef EDT_CRSGRAPH_ZOLTAN_H
#define EDT_CRSGRAPH_ZOLTAN_H

#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;
class Zoltan_LoadBalance;

namespace Epetra_Transform {

class CrsGraph_Zoltan : public SameTypeDistTransform<Epetra_CrsGraph> {

  const std::string partitionMethod_;
  const bool reorder_;
  Zoltan_LoadBalance * lb_;

 public:

  ~CrsGraph_Zoltan() {}

  CrsGraph_Zoltan( const std::string & part_method = std::string("PartKway"),
                   const bool reorder = true,
                   Zoltan_LoadBalance * lb = 0 )
  : partitionMethod_(part_method),
    reorder_(reorder),
    lb_(lb)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

#endif //EDT_CRSGRAPH_ZOLTAN_H
