
#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;

namespace Epetra_Transform {

class CrsGraph_Zoltan : public SameTypeDistTransform<Epetra_CrsGraph> {

  const std::string partitionMethod_;

 public:

  CrsGraph_Zoltan( const std::string & part_method = std::string("PartKway") )
  : partitionMethod_(part_method)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

