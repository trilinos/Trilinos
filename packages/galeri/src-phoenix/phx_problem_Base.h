#ifndef PHX_PROBLEM_BASE_H
#define PHX_PROBLEM_BASE_H

namespace phx {
namespace problem {
  
class Base : public phx::core::Object 
{
  public:
    virtual ~Base() {}

    virtual void integrate(phx::grid::Loadable& domain,
                         Epetra_FECrsMatrix& A,
                         Epetra_FEVector& RHS) = 0;

    virtual void computeNorms(phx::grid::Loadable& domain,
                              const Epetra_MultiVector& solution) = 0;

    virtual void computeNorms(phx::grid::Loadable& domain,
                              const Epetra_MultiVector& solution,
                              double& solNormL2, double& solSemiNormH1,
                              double& exaNormL2, double& exaSemiNormH1,
                              double& errNormL2, double& errSemiNormH1) = 0;

    virtual void integrateOverElement(phx::quadrature::Element& QE,
                                      Epetra_SerialDenseMatrix& ElementLHS, 
                                      Epetra_SerialDenseMatrix& ElementRHS) = 0;

    virtual void computeNormOverElement(phx::quadrature::Element& QE,
                                        Epetra_SerialDenseMatrix& elementSol,
                                        Epetra_SerialDenseMatrix& elementNorm) = 0;

    virtual void computeErrorOverElement(phx::quadrature::Element& QE,
                                         Epetra_SerialDenseMatrix& elementSol,
                                         Epetra_SerialDenseMatrix& elementNorm) = 0;

    virtual void computeNormOverElement(phx::quadrature::Element& QE,
                                        Epetra_SerialDenseMatrix& elementNorm) = 0;
};

} // namespace problem
} // namespace phx
#endif
