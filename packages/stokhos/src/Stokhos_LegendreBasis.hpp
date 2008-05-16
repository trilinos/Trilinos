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

#ifndef STOKHOS_LEGENDREBASIS_HPP
#define STOKHOS_LEGENDREBASIS_HPP

#include "Stokhos_OrthogPolyBasisBase.hpp"

namespace Stokhos {

  template <typename T>
  class LegendreBasis : public OrthogPolyBasisBase<T> {
  public:

    //! Typename of values
    typedef typename OrthogPolyBasisBase<T>::value_type value_type;

    //! Constructor
    LegendreBasis(unsigned int p);

    //! Destructor
    ~LegendreBasis();

    //! Project a polynomial into this basis
    void projectPoly(const Polynomial<T>& poly, std::vector<T>& coeffs) const;

    //! Project derivative of basis polynomial into this basis
    void projectDerivative(unsigned int i, std::vector<T>& coeffs) const;

  private:

    // Prohibit copying
    LegendreBasis(const LegendreBasis&);

    // Prohibit Assignment
    LegendreBasis& operator=(const LegendreBasis& b);

  protected:

    //! Derivative coefficients
    std::vector< std::vector<T> > deriv_coeffs;

  }; // class LegendreBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_LegendreBasisImp.hpp"

#endif
