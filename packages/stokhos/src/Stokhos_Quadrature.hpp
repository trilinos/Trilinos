// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
