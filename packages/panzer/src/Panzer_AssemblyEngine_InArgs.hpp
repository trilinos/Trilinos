#ifndef PANZER_ASSEMBLY_ENGINE_INARGS_HPP
#define PANZER_ASSEMBLY_ENGINE_INARGS_HPP

#include "Teuchos_RCP.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;

namespace panzer {

  struct AssemblyEngineInArgs {

    //! Input arguments for a phalanx fill evaluation.
    /*!

      @param x - solution vector
      @param dxdt - derivative of solution vector wrt time
      @param f - redisual vector
      @param j - jacobian matrix

    */
    AssemblyEngineInArgs(
       const Teuchos::RCP<Epetra_Vector>& in_x = Teuchos::null,
       const Teuchos::RCP<Epetra_Vector>& in_dxdt = Teuchos::null,
       const Teuchos::RCP<Epetra_Vector>& in_f = Teuchos::null,
       const Teuchos::RCP<Epetra_CrsMatrix>& in_j = Teuchos::null,
       double in_alpha = 1.0, double in_beta = 0.0)
      :
      x(in_x),
      dxdt(in_dxdt),
      f(in_f),
      j(in_j),
      alpha(in_alpha),
      beta(in_beta)
    {

    }

    Teuchos::RCP<Epetra_Vector> x;
    Teuchos::RCP<Epetra_Vector> dxdt;
    Teuchos::RCP<Epetra_Vector> f;
    Teuchos::RCP<Epetra_CrsMatrix> j;
    double alpha;
    double beta;

  };

}

#endif
