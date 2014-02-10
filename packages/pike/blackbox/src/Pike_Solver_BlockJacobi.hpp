#ifndef PIKE_SOLVER_BLOCK_JACOBI_HPP
#define PIKE_SOLVER_BLOCK_JACOBI_HPP

#include "Pike_Solver_DefaultBase.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos { template<typename> class Comm; }

namespace pike {

  class BlockJacobi : public pike::SolverDefaultBase {
    
  public:

    BlockJacobi();
    void completeRegistration();
    void stepImplementation();

    void registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  private:
    
    bool barrierTransfers_;
    bool barrierSolves_;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

  };

}

#endif
