#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace pike {

  BlockGaussSeidel::BlockGaussSeidel() 
  {
    validParameters_ = Teuchos::parameterList("pike:BlockGaussSeidel Valid Parameters");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }

  void BlockGaussSeidel::step()
  {}
  
  void BlockGaussSeidel::solve()
  {}

  int BlockGaussSeidel::getNumberOfIterations() const
  { return numberOfIterations_; }
  
  void BlockGaussSeidel::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setParameterList(paramList);
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> BlockGaussSeidel::getValidParameters() const
  { return validParameters_; }

}
