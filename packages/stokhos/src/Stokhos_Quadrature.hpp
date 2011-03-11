// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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

#include <ostream>
#include "Teuchos_Array.hpp"

namespace Stokhos {

  //! Abstract base class for quadrature methods
  template <typename ordinal_type, typename value_type>
  class Quadrature {
  public:

    //! Constructor
    Quadrature() {}

    //! Destructor
    virtual ~Quadrature() {}

    //! Get number of quadrature points
    virtual ordinal_type size() const = 0;

    //! Get quadrature points
    /*!
     * Array is dimensioned Q-by-d where Q is the number of quadrature
     * points and d is the dimension of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> >& 
    getQuadPoints() const = 0;

    //! Get quadrature weights
    /*!
     * Array is of size Q where Q is the number of quadrature points.
     */
    virtual const Teuchos::Array<value_type>& 
    getQuadWeights() const = 0;

    //! Get values of basis at quadrature points
    /*!
     * Array is dimensioned Q-by-P where Q is the number of quadrature
     * points and P is the size of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> > & 
    getBasisAtQuadPoints() const = 0;

    //! Print quadrature data
    virtual std::ostream& print(std::ostream& os) const = 0;

  private:

    // Prohibit copying
    Quadrature(const Quadrature&);

    // Prohibit Assignment
    Quadrature& operator=(const Quadrature& b);

  }; // class Quadrature

  //! Print quadrature object to stream
  template <typename ordinal_type, typename value_type>
  std::ostream& operator << (std::ostream& os,
			     const Quadrature<ordinal_type, value_type>& quad) {
    return quad.print(os);
  }

} // namespace Stokhos

#endif // STOKHOS_QUADRATURE
