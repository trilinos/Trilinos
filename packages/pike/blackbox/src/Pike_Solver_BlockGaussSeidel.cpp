#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"

namespace pike {

  void BlockGaussSeidel::stepImplementation()
  {
    typedef std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::iterator ModelIterator;
    for (ModelIterator m = models_.begin(); m != models_.end(); ++m)
      (*m)->solve();
  }

}
