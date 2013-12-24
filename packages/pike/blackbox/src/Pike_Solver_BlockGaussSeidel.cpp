#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"

namespace pike {

  void BlockGaussSeidel::stepImplementation()
  {
    for (ModelIterator m = models_.begin(); m != models_.end(); ++m) {
     
      // for the model about to be solved, transfer all data to the model

      (*m)->solve();
    }
  }

}
