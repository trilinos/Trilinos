
#ifndef EDT_CRSGRAPH_ZOLTAN_H
#define EDT_CRSGRAPH_ZOLTAN_H

#include <Epetra_Transform.h>

#include <string>

class Epetra_CrsGraph;
class Zoltan_LoadBalance;

namespace EpetraExt {

class CrsGraph_Zoltan : public StructuralSameTypeTransform<Epetra_CrsGraph> {

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

  NewTypePtr operator()( OriginalTypeRef original );

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_ZOLTAN_H
