#ifndef PANZER_ASSEMBLY_ENGINE_INARGS_HPP
#define PANZER_ASSEMBLY_ENGINE_INARGS_HPP

#include "Teuchos_RCP.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;
class Epetra_Map;

namespace Thyra {
template <typename ScalarT> class MultiVectorBase;
template <typename ScalarT> class LinearOpBase;
}

namespace panzer {

  class AssemblyEngineInArgs {
    public:

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

    //! Input arguments for a Thyra phalanx fill evaluation.
    /*!

      @param x - solution vector
      @param dxdt - derivative of solution vector wrt time
      @param f - redisual vector
      @param j - jacobian matrix

    */
    AssemblyEngineInArgs(
       const Teuchos::RCP<Thyra::MultiVectorBase<double> >& in_x,
       const Teuchos::RCP<Thyra::MultiVectorBase<double> >& in_dxdt,
       const Teuchos::RCP<Thyra::MultiVectorBase<double> >& in_f,
       const Teuchos::RCP<Thyra::LinearOpBase<double> >& in_j,
       double in_alpha = 1.0, double in_beta = 0.0)
      :
      th_x(in_x),
      th_dxdt(in_dxdt),
      th_f(in_f),
      th_j(in_j),
      alpha(in_alpha),
      beta(in_beta)
    {
       // thyraToEpetra();
    }

    Teuchos::RCP<Epetra_Vector> x;
    Teuchos::RCP<Epetra_Vector> dxdt;
    Teuchos::RCP<Epetra_Vector> f;
    Teuchos::RCP<Epetra_CrsMatrix> j;

    Teuchos::RCP<Thyra::MultiVectorBase<double> > th_x;
    Teuchos::RCP<Thyra::MultiVectorBase<double> > th_dxdt;
    Teuchos::RCP<Thyra::MultiVectorBase<double> > th_f;
    Teuchos::RCP<Thyra::LinearOpBase<double> > th_j;

    double alpha;
    double beta;

    void epetraToThyra(const Teuchos::RCP<const Epetra_Map> & range,
                       const Teuchos::RCP<const Epetra_Map> & domain);

  private:
    void thyraToEpetra();
  };

}

#endif
