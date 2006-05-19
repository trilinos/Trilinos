#ifndef PHX_PROBLEM_BASE_H
#define PHX_PROBLEM_BASE_H

namespace phx {
namespace problem {
  
class Base : public phx::core::Object 
{
  public:
    virtual ~Base() {}

    virtual void integrate(phx::quadrature::Element& QE,
                           Epetra_SerialDenseMatrix& ElementLHS, 
                           Epetra_SerialDenseMatrix& ElementRHS) = 0;

    virtual void computeNorm(phx::quadrature::Element& QE,
                             Epetra_SerialDenseMatrix& elementSol,
                             Epetra_SerialDenseMatrix& elementNorm) = 0;

    virtual void computeNorm(phx::quadrature::Element& QE,
                             Epetra_SerialDenseMatrix& elementSol,
                             double (*exactSolution)(const char& what, const double& x, 
                                                     const double& y, const double& z),
                             Epetra_SerialDenseMatrix& elementNorm) = 0;

    virtual void computeNorm(phx::quadrature::Element& QE,
                             double (*exactSolution)(const char& what, const double& x, 
                                                     const double& y, const double& z),
                             Epetra_SerialDenseMatrix& elementNorm) = 0;
};

} // namespace problem
} // namespace phx
#endif
