#include "Pike_Solver_BlockJacobi.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_DataTransfer.hpp"

namespace pike {

  BlockJacobi::BlockJacobi()
  {
    this->getNonconstValidParameters()->set("Type","Block Jacobi");
  }

  void BlockJacobi::stepImplementation()
  {
    for (TransferIterator t = transfers_.begin(); t != transfers_.end(); ++t)
      (*t)->doTransfer(*this);
    
    for (ModelIterator m = models_.begin(); m != models_.end(); ++m)
      (*m)->solve();
  }

}
