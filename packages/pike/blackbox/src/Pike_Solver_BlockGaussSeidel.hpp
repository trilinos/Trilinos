#ifndef PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
#define PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP

#include "Pike_Solver.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace pike {

  class BlockGaussSeidel : public pike::Solver {

  public:

    BlockGaussSeidel();

    void step();

    void solve();

    int getNumberOfIterations() const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:

    int numberOfIterations_;

    Teuchos::RCP<Teuchos::ParameterList> validParameters_;

  };

}

#endif
