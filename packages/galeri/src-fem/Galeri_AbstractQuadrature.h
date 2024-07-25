// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_ABSTRACTQUADRATURE_H
#define GALERI_ABSTRACTQUADRATURE_H

/*!
 * \file Galeri_AbstractQuadrature.h
 */

#include "Galeri_AbstractVariational.h"

namespace Galeri {
namespace FiniteElements {

/*!
 * \class AbstractQuadrature
 *
 * \brief Interfaces for quadrature over elements.
 *
AbstractQuadrature is a pure virtual class that defines a set of
abstract interfaces to basis and test functions (and their
derivatives), and also furnishes all the tools required to
numerically integrate over an element.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on Apr-05.
 */

class AbstractQuadrature
{
  
public:

  // @{ \name Constructor and destructor
  
  //! Destructor.
  virtual ~AbstractQuadrature() {}
  
  // @}
  // @{ \name Query methods

  //! Returns the number of quadrature node per element.
  virtual int NumQuadrNodes() const = 0;

  //! Returns the number of basis function on the reference element.
  virtual int NumPhiFunctions() const = 0;

  //! Returns the number of test function on the reference element.
  virtual int NumPsiFunctions() const = 0;

  // @}
  // @{ \name Integration methods

  //! Computes the Jacobian at the specified quadrature node.
  virtual void ComputeJacobian(const int QuadrNode,
			      const double* x,
			      const double* y,
			      const double* z) const = 0;
  
  //! Maps the quadrature nodes from the reference element to the actual one.
  virtual void ComputeQuadrNodes(const int QuadrNode, const double* x,
				const double* y, const double* z,
				double& xq, double& yq, double& zq) const = 0;
    
  //! Computes the derivatives at the specified quadrature node.
  virtual void ComputeDerivatives(const int QuadrNode) const = 0;

  //! Computes the weight at the specified quadrature node.
  virtual double QuadrWeight(const int QuadrNode) const = 0;

  //! Computes the determinant of the Jacobian matrix at the quadrature node.
  virtual double DetJacobian(const int QuadrNode) const = 0;

  // @}
  // @{ \name Basis and test functions.
  
  //! Returns the value of the i-th basis function on the reference element.
  virtual double Phi(const int i) const = 0;

  //! Returns the value of the x-derivative i-th basis function on the reference element.
  virtual double PhiX(const int i) const = 0;

  //! Returns the value of the y-derivative i-th basis function on the reference element.
  virtual double PhiY(const int i) const = 0;

  //! Returns the value of the z-derivative i-th basis function on the reference element.
  virtual double PhiZ(const int i) const = 0;

  //! Returns the value of the i-th test function on the reference element.
  virtual double Psi(const int i) const = 0;

  //! Returns the value of the z-derivative i-th test function on the reference element.
  virtual double PsiX(const int i) const = 0;

  //! Returns the value of the y-derivative i-th test function on the reference element.
  virtual double PsiY(const int i) const = 0;

  //! Returns the value of the z-derivative i-th test function on the reference element.
  virtual double PsiZ(const int i) const = 0;
  
  // @} 
};
 
} // namespace FiniteElements
} // namespace Galeri
#endif
