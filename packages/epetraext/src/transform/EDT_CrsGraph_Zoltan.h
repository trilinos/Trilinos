
#ifndef EDT_CRSGRAPH_ZOLTAN_H
#define EDT_CRSGRAPH_ZOLTAN_H

#include <Epetra_Transform.h>

#include <string>

class Epetra_CrsGraph;
class Epetra_Map;

class Zoltan_LoadBalance;

namespace EpetraExt {

class CrsGraph_Zoltan : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  const std::string partitionMethod_;

  Epetra_Map * NewRowMap_;

 public:

  ~CrsGraph_Zoltan();

  CrsGraph_Zoltan( const std::string & part_method = std::string("PartKway") )
  : partitionMethod_(part_method),
    NewRowMap_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_ZOLTAN_H
