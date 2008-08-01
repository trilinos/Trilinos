// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_HERMITEEBASIS2_HPP
#define STOKHOS_HERMITEEBASIS2_HPP

#include "Stokhos_OrthogPolyBasisBase.hpp"

namespace Stokhos {

  template <typename T>
  class HermiteEBasis2 : public OrthogPolyBasisBase<T> {
  public:
    
    //! Typename of values
    typedef typename OrthogPolyBasisBase<T>::value_type value_type;
    
    //! Constructor
    HermiteEBasis2(unsigned int p);
    
    //! Destructor
    ~HermiteEBasis2();
    
    //! Project a polynomial into this basis
    void projectPoly(const Polynomial<T>& poly, std::vector<T>& coeffs) const;

    //! Project derivative of basis polynomial into this basis
    void projectDerivative(unsigned int i, std::vector<T>& coeffs) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const std::vector<T>& point,
			    std::vector<T>& basis_pts) const;

  private:

    // Prohibit copying
    HermiteEBasis2(const HermiteEBasis2&);

    // Prohibit Assignment
    HermiteEBasis2& operator=(const HermiteEBasis2& b);
    
  }; // class HermiteEBasis2

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_HermiteEBasis2Imp.hpp"

#endif
