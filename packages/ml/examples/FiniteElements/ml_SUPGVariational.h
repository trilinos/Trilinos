#ifndef ML_SUPGVARIATIONAL_H
#define ML_SUPGVARIATIONAL_H

/*!
 * \file ml_SUPGVariationa.h
 */

#include "ml_Workspace.h"
#include "ml_AbstractVariational.h"

namespace ML_FiniteElements {

/*!
 * \class SUPGVariational
 *
 * \brief SUPG discretization of an advection-diffusion PDE.
 *
 * This class performs the finite element discretization of a scalar,
 * advection-diffusion PDE, using the SUPG stabilization and the coth
 * formula for the definition of tau. This class works only with triangles
 * and tetrahedra.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on Apr-05.
 */

template<class T>
class SUPGVariational : public AbstractVariational, public T
{
public:

  //! Constructor.
  SUPGVariational(const int NumQuadratureNodes,
                  double (*diff)(const double&, const double&, const double&),
                  double (*bx)(const double&, const double&, const double&),
                  double (*by)(const double&, const double&, const double&),
                  double (*bz)(const double&, const double&, const double&),
                  double (*source)(const double&, const double&, const double&),
                  double (*force)(const double&, const double&, const double&),
                  double (*bc)(const double&, const double&, const double&, const int&)) :
  T(NumQuadratureNodes),
  diff_(diff),
  source_(source),
  force_(force),
  conv_x_(bx),
  conv_y_(by),
  conv_z_(bz),
  bc_(bc)
  {}

  //! Destructor.
  ~SUPGVariational() {}

  //! Evaluates the diffusion coefficient at point (x, y, z).
  inline double diff(const double x, const double y, const double z) const
  {
    return (diff_(x, y, z));
  }

  //! Evaluates the source term at point (x, y, z).
  inline double source(const double x, const double y, const double z) const
  {
    return (source_(x, y, z));
  }

  //! Evaluates the force term at point (x, y, z).
  inline double force(const double x, const double y, const double z) const
  {
    return (force_(x, y, z));
  }

  //! Evaluates the x-component of the convective term at point (x, y, z).
  inline double conv_x(const double x, const double y, const double z) const
  {
    return (conv_x_(x, y, z));
  }

  //! Evaluates the y-component of the convective term at point (x, y, z).
  inline double conv_y(const double x, const double y, const double z) const
  {
    return (conv_y_(x, y, z));
  }

  //! Evaluates the z-component of the convective term at point (x, y, z).
  inline double conv_z(const double x, const double y, const double z) const
  {
    return (conv_z_(x, y, z));
  }

  virtual int IntegrateOverElement(const AbstractVariational& Variational,
				   const double* x, const double* y, const double* z,
                                   const double* data,
				   double* ElementMatrix, double* ElementRHS) const
  {
    double xq, yq, zq;
    int size = NumPhiFunctions();
    double h = data[0];
    
    // zero out local matrix and rhs
    
    for (int i = 0 ; i < size * size ; i++) ElementMatrix[i] = 0.0;
    for (int i = 0 ; i < size ; i++)        ElementRHS[i] = 0.0;

    // cycle over all quadrature nodes

    for (int ii = 0 ; ii < NumQuadrNodes() ; ii++) 
    {
      ComputeQuadrNodes(ii,x, y, z, xq, yq, zq);
      ComputeJacobian(ii,x, y, z);
      ComputeDerivatives(ii);

      // compute local Peclet number for this element
      // using the coth formula
      double Norm = ConvNorm(xq, yq, zq);
      double Pe = 0.5 * Norm * h / diff(xq, yq, zq);
      if (Norm == 0.0)
        tau_ = 0.0;
      else
        tau_ = 0.5 * (h / Norm) * (cosh(Pe) / sinh(Pe) - 1.0 / Pe);

      for (int i = 0 ; i < NumPhiFunctions() ; ++i) 
      {
        for (int j = 0 ; j < NumPsiFunctions() ; ++j) 
        {
          ElementMatrix[j + size * i] +=
            QuadrWeight(ii) * DetJacobian(ii) * 
            Variational.LHS(Phi(i), Psi(j), PhiX(i), PsiX(j),
                            PhiY(i), PsiY(j), PhiZ(i), PsiZ(j),
                            xq, yq, zq);
        }
        ElementRHS[i] += QuadrWeight(ii) * DetJacobian(ii) *
          Variational.RHS(Psi(i), PsiX(i), PsiY(i), PsiZ(i), 
                          xq, yq, zq);
      }
    }

    return 0;
  }

  virtual int ElementNorm(const double* LocalSol, const double* x, 
                          const double* y, const double* z, double* Norm) const
  {
    double xq, yq, zq;
    double exact[4];

    for (int ii = 0 ; ii < NumQuadrNodes() ; ii++) 
    {
      ComputeQuadrNodes(ii,x, y, z, xq, yq, zq );
      ComputeJacobian(ii,x, y, z);
      ComputeDerivatives(ii);

      double GlobalWeight = QuadrWeight(ii) * DetJacobian(ii);

      double sol      = 0.0, sol_derx = 0.0;
      double sol_dery = 0.0, sol_derz = 0.0;

      for (int k = 0 ; k < NumPhiFunctions() ; ++k)
      {
        sol      += Phi(k)  * LocalSol[k];
        sol_derx += PhiX(k) * LocalSol[k];
        sol_dery += PhiY(k) * LocalSol[k];
        sol_derz += PhiZ(k) * LocalSol[k];
      }

      Norm[0] += GlobalWeight*sol*sol;
      Norm[1] += GlobalWeight*(sol_derx*sol_derx +
                               sol_dery*sol_dery +
                               sol_derz*sol_derz);
    }

    return 0;
  }

  virtual int ElementNorm(int (*ExactSolution)(double, double, double, double *),
			  const double* x, const double* y, const double* z,
			  double* Norm) const
  {
    double xq, yq, zq;
    double exact[4];

    for (int ii = 0 ; ii < NumQuadrNodes() ; ii++) 
    {
      ComputeQuadrNodes(ii, x, y, z, xq, yq, zq );
      ComputeJacobian(ii, x, y, z);
      ComputeDerivatives(ii);

      double GlobalWeight = QuadrWeight(ii) * DetJacobian(ii);

      (*ExactSolution)(xq, yq, zq, exact);

      Norm[0] += GlobalWeight * exact[0] * exact[0];
      Norm[1] += GlobalWeight * (exact[1] * exact[1] +
                                 exact[2] * exact[2] +
                                 exact[3] * exact[3]);
    }

    return 0;
  }
  
  virtual int ElementNorm(const double* LocalSol,
			  int (*ExactSolution)(double, double, double, double *),
			  const double* x, const double* y, const double* z, double * Norm) const
  {
    double xq, yq, zq;
    double exact[4];

    for (int ii = 0 ; ii < NumQuadrNodes() ; ii++) 
    {
      ComputeQuadrNodes(ii, x, y, z, xq, yq, zq );
      ComputeJacobian(ii, x, y, z);
      ComputeDerivatives(ii);

      double GlobalWeight = QuadrWeight(ii) * DetJacobian(ii);

      double diff      = 0.0, diff_derx = 0.0;
      double diff_dery = 0.0, diff_derz = 0.0;

      for (int k = 0 ; k < NumPhiFunctions() ; ++k) 
      {
        diff      += Phi(k)  * LocalSol[k];
        diff_derx += PhiX(k) * LocalSol[k];
        diff_dery += PhiY(k) * LocalSol[k];
        diff_derz += PhiZ(k) * LocalSol[k];
      }

      (*ExactSolution)(xq, yq, zq,exact);

      diff      -= exact[0];
      diff_derx -= exact[1];
      diff_dery -= exact[2];
      diff_derz -= exact[3];

      Norm[0] += GlobalWeight * diff * diff;
      Norm[1] += GlobalWeight * (diff_derx * diff_derx +
                                 diff_dery * diff_dery +
                                 diff_derz * diff_derz);
    } 
    return(0);
  }

  inline double LHS(const double Phi, const double Psi,
                    const double PhiX, const double PsiX,
                    const double PhiY, const double PsiY,
                    const double PhiZ, const double PsiZ,
                    const double x, const double y, const double z) const
  {
    double res;
    // Galerkin contribution
    res = diff(x,y,z) * PhiX * PsiX +
          diff(x,y,z) * PhiY * PsiY +
          diff(x,y,z) * PhiZ * PsiZ +
          source(x,y,z) * Phi * Psi +
          conv_x(x,y,z) * PhiX * Psi + 
          conv_y(x,y,z) * PhiY * Psi +
          conv_z(x,y,z) * PhiZ * Psi;
    // SUPG stabilization
    res += tau_ * ((conv_x(z, y, z) * PsiX + 
                    conv_y(x, y, z) * PsiY + 
                    conv_z(x, y, z) * PsiZ) *
                   (conv_x(z, y, z) * PhiX + 
                    conv_y(x, y, z) * PhiY + 
                    conv_z(x, y, z) * PhiZ));
    return(res);
  }

  inline double RHS(const double Psi, const double PsiX, 
                    const double PsiY, const double PsiZ,
                    const double x, const double y, const double z) const
  {
    double res;
    // Galerkin contribution
    res = force(x,y,z)*Psi;
    // SUPG stabilization
    res += tau_ * (conv_x(x, y, z) * PsiX +
                   conv_y(x, y, z) * PsiY +
                   conv_z(x, y, z) * PsiZ);
    return(res);
  }

  int BC(const int PatchID) const
  {
    return(ML_DIRICHLET);
  }

  double BC(const double x, const double y, const double z, const int Patch) const
  {
    return(bc_(x, y, z, Patch));
  }

private:
  //! Computes the norm of the convective term a point (x, y, z).
  inline double ConvNorm(const double x, const double y, const double z) const
  {
    double cx = conv_x(z, y, z);
    double cy = conv_x(z, y, z);
    double cz = conv_x(z, y, z);
    return (sqrt(cx * cx + cy * cy + cz * cz));
  }

  double (*diff_)(const double& x, const double& y, const double& z);
  double (*source_)(const double& x, const double& y, const double& z);
  double (*conv_x_)(const double& x, const double& y, const double& z);
  double (*conv_y_)(const double& x, const double& y, const double& z);
  double (*conv_z_)(const double& x, const double& y, const double& z);
  double (*force_)(const double& x, const double& y, const double& z);
  double (*bc_)(const double& x, const double& y, const double& z, const int& Patch);
  mutable double tau_;
};

} // namespace ML_FiniteElements
#endif
