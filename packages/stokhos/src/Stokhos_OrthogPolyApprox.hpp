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

#ifndef STOKHOS_ORTHOGPOLYAPPROX_HPP
#define STOKHOS_ORTHOGPOLYAPPROX_HPP

#include "Teuchos_RCP.hpp"             // class data member
#include "Teuchos_Array.hpp"           // class data member
#include "Stokhos_OrthogPolyBasis.hpp" // class data member
#include <ostream>	               // for std::ostream

namespace Stokhos {

  /*! 
   * \brief Class to store coefficients of a projection onto an orthogonal
   * polynomial basis.
   */
  template <typename ordinal_type, typename value_type>
  class OrthogPolyApprox {
  public:

    //! Constructor with no basis
    /*!
     * Use with care!  Sets size to 1 to store constant term
     */
    OrthogPolyApprox();

    //! Constructor
    OrthogPolyApprox(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis);
    
    //! Constructor with supplied size \c sz
    /*!
     * Normally \c sz should equal the basis size, however this is not
     * enforced since other situations can arise, e.g., using a size of 1
     * for a constant expansion.
     */
    OrthogPolyApprox(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
		     ordinal_type sz);
    
    //! Copy constructor
    OrthogPolyApprox(const OrthogPolyApprox& x);
    
    //! Destructor
    ~OrthogPolyApprox();
    
    //! Assignment operator (deep copy)
    OrthogPolyApprox& operator=(const OrthogPolyApprox& x);

    //! Intialize coefficients to value
    void init(const value_type& v);

    //! Return basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > 
    basis() const;

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.  Coefficients are preserved.
     */
    void reset(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& new_basis);

    //! Return size
    ordinal_type size() const;

    //! Return coefficient array
    value_type* coeff();

    //! Return coefficient array
    const value_type* coeff() const;

    //! Array access
    value_type& operator[](ordinal_type i);

    //! Array access
    const value_type& operator[](ordinal_type i) const;

    //! Get coefficient term for given dimension and order
    value_type& term(ordinal_type dimension, ordinal_type order);

    //! Get coefficient term for given dimension and order
    const value_type& term(ordinal_type dimension, ordinal_type order) const;

    //! Evaluate polynomial approximation at a point
    value_type evaluate(const Teuchos::Array<value_type>& point) const;

    //! Evaluate polynomial approximation at a point with supplied basis values
    value_type evaluate(const Teuchos::Array<value_type>& point,
			const Teuchos::Array<value_type>& basis_vals) const;

    //! Compute mean of expansion
    value_type mean() const;

    //! Compute standard deviation of expansion
    value_type standard_deviation() const;

    //! Print approximation in basis
    std::ostream& print(std::ostream& os) const;

  protected:

    //! Basis expansion is relative to
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis_;

    //! OrthogPolyApprox coefficients
    Teuchos::Array<value_type> coeff_;
    
  }; // class OrthogPolyApprox

  //! Prints the array of coefficients (more compact than print())
  template <typename ordinal_type, typename value_type> std::ostream& 
  operator << (std::ostream& os, 
	       const OrthogPolyApprox<ordinal_type,value_type>& a);

} // namespace Stokhos

#include "Stokhos_OrthogPolyApproxImp.hpp"

#endif //  STOKHOS_ORTHOGPOLYAPPROX_HPP
