// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_GALERKINVARIATIONAL_H
#define GALERI_GALERKINVARIATIONAL_H

/*!
 * \file Galeri_GalerkinVariational.h
 */

#include "Galeri_Workspace.h"
#include "Galeri_AbstractVariational.h"

namespace Galeri {
namespace FiniteElements {

/*!
 * \class GalerkinVariational
 *
 * \brief Defines a pure Galerkin variational form of a scalar PDE.
 *
 * This class defines a pure Galerkin variational form of a second order,
 * symmetric scalar PDE, discretized using Lagrange finite elements. The class is
 * templated with an AbstractQuadrature class, which will be used to 
 * specify the quadrature formula, and the values of test and basis functions
 * at the quadrature node. The constructor requires function pointers, that
 * specify the values of the coefficients.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on Apr-05.
 */

template<class T>
class GalerkinVariational : public AbstractVariational, public T
{
public:

  //! Constructor.
  GalerkinVariational(const int NumQuadratureNodes,
                      double (*diff)(const double&, const double&, const double&),
                      double (*source)(const double&, const double&, const double&),
                      double (*force)(const double&, const double&, const double&),
                      double (*bc)(const double&, const double&, const double&, const int&),
                      int (*bc_type)(const int&)):
    T(NumQuadratureNodes),
    diff_(diff),
    source_(source),
    force_(force),
    bc_(bc),
    bc_type_(bc_type)
  {}

  //! Destructor.  
  ~GalerkinVariational() {}

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

  //! Integrates the variational form and the right-hand side.
  virtual int IntegrateOverElement(const AbstractVariational& Variational,
				   const double* x, const double* y, const double* z,
                                   const double* data,
				   double* ElementMatrix, double* ElementRHS) const
  {
    double xq, yq, zq;
    int size = T::NumPhiFunctions();
    //double h = data[0];
    
    // zero out local matrix and rhs
    
    for (int i = 0 ; i < size * size ; i++) ElementMatrix[i] = 0.0;
    for (int i = 0 ; i < size ; i++)        ElementRHS[i] = 0.0;

    // cycle over all quadrature nodes

    for (int ii = 0 ; ii < T::NumQuadrNodes() ; ii++) 
    {
      T::ComputeQuadrNodes(ii,x, y, z, xq, yq, zq);
      T::ComputeJacobian(ii,x, y, z);
      T::ComputeDerivatives(ii);

      for (int i = 0 ; i < T::NumPhiFunctions() ; ++i) 
      {
        for (int j = 0 ; j < T::NumPsiFunctions() ; ++j) 
        {
          ElementMatrix[j + size * i] +=
            T::QuadrWeight(ii) * T::DetJacobian(ii) * 
            Variational.LHS(T::Phi(i), T::Psi(j), T::PhiX(i), T::PsiX(j),
                            T::PhiY(i), T::PsiY(j), T::PhiZ(i), T::PsiZ(j),
                            xq, yq, zq);
        }
        ElementRHS[i] += T::QuadrWeight(ii) * T::DetJacobian(ii) *
          Variational.RHS(T::Psi(i), T::PsiX(i), T::PsiY(i), T::PsiZ(i), 
                          xq, yq, zq);
      }
    }

    return 0;
  }

  //! Computes the norm of the numerical solution over an element.
  virtual int ElementNorm(const double* LocalSol, const double* x, 
                          const double* y, const double* z, double* Norm) const
  {
    double xq, yq, zq;
    //double exact[4];

    for (int ii = 0 ; ii < T::NumQuadrNodes() ; ii++) 
    {
      T::ComputeQuadrNodes(ii,x, y, z, xq, yq, zq );
      T::ComputeJacobian(ii,x, y, z);
      T::ComputeDerivatives(ii);

      double GlobalWeight = T::QuadrWeight(ii) * T::DetJacobian(ii);

      double sol      = 0.0, sol_derx = 0.0;
      double sol_dery = 0.0, sol_derz = 0.0;

      for (int k = 0 ; k < T::NumPhiFunctions() ; ++k)
      {
        sol      += T::Phi(k)  * LocalSol[k];
        sol_derx += T::PhiX(k) * LocalSol[k];
        sol_dery += T::PhiY(k) * LocalSol[k];
        sol_derz += T::PhiZ(k) * LocalSol[k];
      }

      Norm[0] += GlobalWeight*sol*sol;
      Norm[1] += GlobalWeight*(sol_derx*sol_derx +
                               sol_dery*sol_dery +
                               sol_derz*sol_derz);
    }

    return 0;
  }

  //! Computes the norm of the exact solution over an element.
  virtual int ElementNorm(int (*ExactSolution)(double, double, double, double *),
			  const double* x, const double* y, const double* z,
			  double* Norm) const
  {
    double xq, yq, zq;
    double exact[4];

    for (int ii = 0 ; ii < T::NumQuadrNodes() ; ii++) 
    {
      T::ComputeQuadrNodes(ii, x, y, z, xq, yq, zq );
      T::ComputeJacobian(ii, x, y, z);
      T::ComputeDerivatives(ii);

      double GlobalWeight = T::QuadrWeight(ii) * T::DetJacobian(ii);

      (*ExactSolution)(xq, yq, zq, exact);

      Norm[0] += GlobalWeight * exact[0] * exact[0];
      Norm[1] += GlobalWeight * (exact[1] * exact[1] +
                                 exact[2] * exact[2] +
                                 exact[3] * exact[3]);
    }

    return 0;
  }
  
  //! Computes the norm of the error over an element.
  virtual int ElementNorm(const double* LocalSol,
			  int (*ExactSolution)(double, double, double, double *),
			  const double* x, const double* y, const double* z, double * Norm) const
  {
    double xq, yq, zq;
    double exact[4];

    for (int ii = 0 ; ii < T::NumQuadrNodes() ; ii++) 
    {
      T::ComputeQuadrNodes(ii, x, y, z, xq, yq, zq );
      T::ComputeJacobian(ii, x, y, z);
      T::ComputeDerivatives(ii);

      double GlobalWeight = T::QuadrWeight(ii) * T::DetJacobian(ii);

      double diff      = 0.0, diff_derx = 0.0;
      double diff_dery = 0.0, diff_derz = 0.0;

      for (int k = 0 ; k < T::NumPhiFunctions() ; ++k) 
      {
        diff      += T::Phi(k)  * LocalSol[k];
        diff_derx += T::PhiX(k) * LocalSol[k];
        diff_dery += T::PhiY(k) * LocalSol[k];
        diff_derz += T::PhiZ(k) * LocalSol[k];
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

  //! Evaluates the left-hand side at point (x, y, z).
  inline double LHS(const double Phi, const double Psi,
                    const double PhiX, const double PsiX,
                    const double PhiY, const double PsiY,
                    const double PhiZ, const double PsiZ,
                    const double x, const double y, const double z) const
  {
    return(diff(x,y,z) * PhiX * PsiX +
           diff(x,y,z) * PhiY * PsiY +
           diff(x,y,z) * PhiZ * PsiZ +
           source(x,y,z) * Phi * Psi);
  }

  //! Evaluates the right-hand side at point (x, y, z).
  inline double RHS(const double Psi, const double PsiX, 
                    const double PsiY, const double PsiZ,
                    const double x, const double y, const double z) const
  {
    return(force(x,y,z)*Psi);
  }

  //! Returns the boundary condition type of the specified patch.
  int BC(const int PatchID) const
  {
    return(bc_type_(PatchID));
  }

  //! Returns the value of the boundary condition at point (x, y, z).
  double BC(const double x, const double y, const double z, const int PatchID) const
  {
    return(bc_(x, y, z, PatchID));
  }

private:
  double (*diff_)(const double& x, const double& y, const double& z);
  double (*source_)(const double& x, const double& y, const double& z);
  double (*force_)(const double& x, const double& y, const double& z);
  double (*bc_)(const double& x, const double& y, const double& z, const int& Patch);
  int (*bc_type_)(const int& Patch);
};

} // namespace FiniteElements
} // namespace Galeri
#endif
