#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_DataTransfer.hpp"

namespace pike {

  void BlockGaussSeidel::completeRegistration()
  {
    this->pike::SolverDefaultImplementation::completeRegistration();

    modelAndTransfers_.resize(models_.size());
    for (std::size_t  i = 0; i < modelAndTransfers_.size(); ++i) {
      modelAndTransfers_[i].first = models_[i];
      modelNameToIndex_[models_[i]->name()] = i;
    }

    for (TransferIterator t = transfers_.begin(); t != transfers_.end(); ++t) {
      const std::vector<std::string> targetModels = (*t)->getTargetModelNames();
      for (std::vector<std::string>::const_iterator n = targetModels.begin(); 
	   n != targetModels.end(); ++n) {
	modelAndTransfers_[modelNameToIndex_[*n]].second.push_back(*t);
      }
    }
  }

  void BlockGaussSeidel::stepImplementation()
  {
    typedef std::vector<std::pair<Teuchos::RCP<pike::BlackBoxModelEvaluator>,std::vector<Teuchos::RCP<pike::DataTransfer> > > >::iterator GSIterator;

    for (GSIterator m = modelAndTransfers_.begin(); m != modelAndTransfers_.end(); ++m) {
 
      // for the model about to be solved, transfer all data to the model
      for (std::vector<Teuchos::RCP<pike::DataTransfer> >::iterator t = m->second.begin();
	   t != m->second.end(); ++t)
	(*t)->doTransfer(*this);

      m->first->solve();
    }
  }

}
