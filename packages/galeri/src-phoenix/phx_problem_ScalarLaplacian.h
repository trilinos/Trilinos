#ifndef PHX_PROBLEM_SCALAR_LAPLACIAN_H
#define PHX_PROBLEM_SCALAR_LAPLACIAN_H

#include "phx_problem_Base.h"

namespace phx {
namespace problem {
  
class ScalarLaplacian : public Base
{
  public:
  ScalarLaplacian() {}

  ~ScalarLaplacian() {}

  virtual void integrate(phx::quadrature::Element& QE,
                         Epetra_SerialDenseMatrix& ElementLHS, 
                         Epetra_SerialDenseMatrix& ElementRHS)
  {
    for (int i = 0; i < ElementLHS.M(); ++i)
      for (int j = 0; j < ElementLHS.N(); ++j)
        ElementLHS(i, j) = 0.0;

    for (int i = 0; i < ElementRHS.M(); ++i)
      for (int j = 0; j < ElementRHS.N(); ++j)
        ElementRHS(i, j) = 0.0;

    // cycle over all quadrature nodes

    for (int ii = 0 ; ii < QE.getNumQuadrNodes() ; ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      for (int i = 0 ; i < QE.getNumBasisFunctions() ; ++i) 
      {
        const double& phi_i = QE.getPhi(i);
        const double& phi_x_i = QE.getPhiX(i);
        const double& phi_y_i = QE.getPhiY(i);
        const double& phi_z_i = QE.getPhiZ(i);

        for (int j = 0 ; j < QE.getNumBasisFunctions() ; ++j) 
        {
          const double& phi_j = QE.getPhi(j);

          const double& phi_x_j = QE.getPhiX(j);
          const double& phi_y_j = QE.getPhiY(j);
          const double& phi_z_j = QE.getPhiZ(j);

          double contrib = 0.0;
          contrib += phi_x_i * phi_x_j + phi_y_i * phi_y_j + phi_z_i * phi_z_j;
          contrib *= weight * det;

          ElementLHS(i, j) += contrib;
        }

        ElementRHS(i, 0) += weight * det * phi_i;
      }
    }
  }

  virtual void computeNorm(phx::quadrature::Element& QE,
                           Epetra_SerialDenseMatrix& elementSol,
                           Epetra_SerialDenseMatrix& elementNorm)
  {
    elementNorm(0, 0) = 0.0;
    elementNorm(1, 0) = 0.0;

    for (int ii = 0 ; ii < QE.getNumQuadrNodes() ; ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      double sol      = 0.0, sol_derx = 0.0;
      double sol_dery = 0.0, sol_derz = 0.0;

      for (int k = 0 ; k < QE.getNumBasisFunctions() ; ++k)
      {
        sol      += QE.getPhi(k)  * elementSol(k, 0);
        sol_derx += QE.getPhiX(k) * elementSol(k, 0);
        sol_dery += QE.getPhiY(k) * elementSol(k, 0);
        sol_derz += QE.getPhiZ(k) * elementSol(k, 0);
      }

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

  virtual void computeNorm(phx::quadrature::Element& QE,
                           Epetra_SerialDenseMatrix& elementSol,
                           double (*exactSolution)(const char& what, const double& x, 
                                                   const double& y, const double& z),
                           Epetra_SerialDenseMatrix& elementNorm)
  {
    elementNorm(0, 0) = 0.0;
    elementNorm(1, 0) = 0.0;

    for (int ii = 0 ; ii < QE.getNumQuadrNodes() ; ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      double sol      = 0.0, sol_derx = 0.0;
      double sol_dery = 0.0, sol_derz = 0.0;

      for (int k = 0 ; k < QE.getNumBasisFunctions() ; ++k)
      {
        sol      += QE.getPhi(k)  * elementSol(k, 0);
        sol_derx += QE.getPhiX(k) * elementSol(k, 0);
        sol_dery += QE.getPhiY(k) * elementSol(k, 0);
        sol_derz += QE.getPhiZ(k) * elementSol(k, 0);
      }

      sol      -= exactSolution('f', xq, yq, zq);
      sol_derx -= exactSolution('x', xq, yq, zq);
      sol_dery -= exactSolution('y', xq, yq, zq);
      sol_derz -= exactSolution('z', xq, yq, zq);

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

  virtual void computeNorm(phx::quadrature::Element& QE,
                           double (*exactSolution)(const char& what, const double& x, 
                                                   const double& y, const double& z),
                           Epetra_SerialDenseMatrix& elementNorm)
  {
    elementNorm(0, 0) = 0.0;
    elementNorm(1, 0) = 0.0;

    for (int ii = 0 ; ii < QE.getNumQuadrNodes() ; ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      double sol      = exactSolution('f', xq, yq, zq);
      double sol_derx = exactSolution('x', xq, yq, zq);
      double sol_dery = exactSolution('y', xq, yq, zq);
      double sol_derz = exactSolution('z', xq, yq, zq);

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

}; // class Base

} // namespace problem
} // namespace phx
#endif

