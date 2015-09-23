#ifndef PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
#define PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP

#include "Pike_Solver_DefaultBase.hpp"
#include <vector>
#include <map>
#include <utility>

namespace pike {

  class BlockGaussSeidel : public pike::SolverDefaultBase {
    
  public:

    BlockGaussSeidel();

    void completeRegistration();

    void stepImplementation();

    void registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  private:

    //! Maps the name of a model to the corresponding index in the models vector.
    std::map<std::string,std::size_t> modelNameToIndex_;

    //! Binds each model to a vector of tranfers where the target of the data transfer is the corresponding model.
    std::vector<std::pair<Teuchos::RCP<pike::BlackBoxModelEvaluator>,std::vector<Teuchos::RCP<pike::DataTransfer> > > > modelAndTransfers_;
    
    bool barrierTransfers_;
    bool barrierSolves_;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  };

}

#endif
