
#include <Epetra_DistTransform.h>

class Epetra_CrsGraph;
class Epetra_BlockMap;

namespace Epetra_Transform {

class CrsGraph_View : public SameTypeDistTransform<Epetra_CrsGraph> {

  const Epetra_BlockMap & OrigMap_;
  const Epetra_BlockMap & NewMap_;

 public:

  ~CrsGraph_View() {}

  CrsGraph_View( const Epetra_BlockMap & orig_map,
                 const Epetra_BlockMap & new_map )
  : OrigMap_(orig_map),
    NewMap_(new_map)
  {}

  std::auto_ptr<Epetra_CrsGraph> operator()( const Epetra_CrsGraph & original );

};

} //namespace Epetra_Transform

