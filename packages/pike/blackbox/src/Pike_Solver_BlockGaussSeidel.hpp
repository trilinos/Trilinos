#ifndef PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
#define PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP

#include "Pike_Solver_DefaultImpl.hpp"
#include <vector>
#include <map>
#include <utility>

namespace pike {

  class BlockGaussSeidel : public pike::SolverDefaultImpl {
    
  public:
    
    void completeRegistration();

    void stepImplementation();

  private:

    //! Maps the name of a model to the corresponding index in the models vector.
    std::map<std::string,std::size_t> modelNameToIndex_;

    //! Binds each model to a vector of tranfers where the target of the data transfer is the corresponding model.
    std::vector<std::pair<Teuchos::RCP<pike::BlackBoxModelEvaluator>,std::vector<Teuchos::RCP<pike::DataTransfer> > > > modelAndTransfers_;
    
  };

}

#endif
