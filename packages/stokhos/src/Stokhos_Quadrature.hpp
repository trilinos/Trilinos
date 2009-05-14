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

#ifndef STOKHOS_QUADRATURE
#define STOKHOS_QUADRATURE

#include <vector>

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class Quadrature {
  public:

    //! Constructor
    Quadrature() {}

    //! Destructor
    virtual ~Quadrature() {}

    //! Get quadrature points
    virtual const std::vector< std::vector<value_type> >& 
    getQuadPoints() const = 0;

    //! Get quadrature weights
    virtual const std::vector<value_type>& 
    getQuadWeights() const = 0;

    //! Get values of basis at quadrature points
    virtual const std::vector< std::vector<value_type> > & 
    getBasisAtQuadPoints() const = 0;

  private:

    // Prohibit copying
    Quadrature(const Quadrature&);

    // Prohibit Assignment
    Quadrature& operator=(const Quadrature& b);

  }; // class Quadrature

} // namespace Stokhos

#endif // STOKHOS_QUADRATURE
