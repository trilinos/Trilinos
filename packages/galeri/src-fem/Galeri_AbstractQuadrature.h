// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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
