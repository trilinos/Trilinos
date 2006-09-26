// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef GAUSSIANQUADRATURE2_HPP
#define GAUSSIANQUADRATURE2_HPP

#include "AbstractQuadrature.hpp"

/*!
 * \brief Two point Gaussian quadrature for integrating functions over the
 * interval \f$(-1,1)\f$.
 */
class GaussianQuadrature2 : public AbstractQuadrature {
public:

  //! Default constructor
  GaussianQuadrature2();

  //! Destructor
  virtual ~GaussianQuadrature2();

  //! Return the number of quadrature points
  virtual unsigned int numPoints() const;

  //! Return the quadrature points
  virtual const std::vector<double>& quadPoints() const;

  //! Return the weights
  virtual const std::vector<double>& weights() const;

private:

  //! Private to prohibit copying
  GaussianQuadrature2(const GaussianQuadrature2&);

  //! Private to prohibit copying
  GaussianQuadrature2& operator=(const GaussianQuadrature2&);

protected:

  //! Quad points
  std::vector<double> qp;

  //! Weights
  std::vector<double> w;

};

#endif // GAUSSIANQUADRATURE2_HPP
