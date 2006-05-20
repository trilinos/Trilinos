#ifndef PHX_PROBLEM_BASE_H
#define PHX_PROBLEM_BASE_H

namespace phx {
namespace problem {
  
class Base : public phx::core::Object 
{
  public:
    virtual ~Base() {}

    virtual void integrate(RefCountPtr<phx::grid::Loadable> domain,
                           const int numDimensions) = 0;

    virtual void computeNorms(RefCountPtr<phx::grid::Loadable> domain,
                              const int numDimensions,
                              const Epetra_MultiVector& solution,
                              const bool print = true) = 0;

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
