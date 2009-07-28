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

#ifndef STOKHOS_ORTHOGPOLYAPPROX_HPP
#define STOKHOS_ORTHOGPOLYAPPROX_HPP

#include <ostream>	// for std::ostream

#include "Teuchos_RCP.hpp"

#include "Stokhos_Polynomial.hpp"

namespace Stokhos {

  //! General polynomial class
  template <typename ordinal_type, typename value_type>
  class OrthogPolyApprox {
  public:

    //! Default constructor
    OrthogPolyApprox();
    
    //! Constructor with supplied value \c x
    OrthogPolyApprox(const value_type& x);
    
    //! Constructor with size \c sz and value \c x
    OrthogPolyApprox(ordinal_type sz, const value_type& x);
    
    //! Constructor with size \c sz
    OrthogPolyApprox(ordinal_type sz);
    
    //! Constructor with basis, size \c sz and length \c l
    OrthogPolyApprox(ordinal_type sz, ordinal_type l);
    
    //! Copy constructor
    OrthogPolyApprox(const OrthogPolyApprox& x);
    
    //! Destructor
    ~OrthogPolyApprox();
    
    //! Assignment operator (deep copy)
    OrthogPolyApprox& operator=(const OrthogPolyApprox& x);

    //! Resize to size \c sz
    /*!
     * Coefficients are preserved.
     */
    void resize(ordinal_type sz);

    //! Reserve space for a size \c sz expansion
    /*!
     * Coefficients are preserved.
     */
    void reserve(ordinal_type sz);

    //! Return size
    ordinal_type size() const { return coeff_.size(); }

    //! Return coefficient array
    value_type* coeff() { return &coeff_[0]; }

    //! Return coefficient array
    const value_type* coeff() const { return &coeff_[0]; }

    //! Array access
    value_type& operator[](ordinal_type i) { return coeff_[i]; }

    //! Array access
    const value_type& operator[](ordinal_type i) const { return coeff_[i]; }

    //! Get term
    template <typename BasisT>
    value_type& 
    term(const BasisT& basis, 
	 int i0 = -1, int i1 = -1, int i2 = -1, int i3 = -1, int i4 = -1,
	 int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1);

    //! Get term
    template <typename BasisT>
    const value_type& 
    term(const BasisT& basis, 
	 int i0 = -1, int i1 = -1, int i2 = -1, int i3 = -1, int i4 = -1,
	 int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1) const;

    template <typename BasisT>
    value_type& 
    term2(const BasisT& basis, ordinal_type dimension, ordinal_type order);

    template <typename BasisT>
    const value_type& 
    term2(const BasisT& basis, ordinal_type dimension, ordinal_type order) const;

    //! Write polynomial approximation in standard basis
    template <typename BasisT>
    Polynomial<value_type> toStandardBasis(const BasisT& basis) const;

    //! Evaluate polynomial approximation at a point
    template <typename BasisT>
    value_type 
    evaluate(const BasisT& basis, const std::vector<value_type>& point) const;

    //! Evaluate polynomial approximation at a point with supplied basis values
    template <typename BasisT>
    value_type 
    evaluate(const BasisT& basis, const std::vector<value_type>& point,
	     const std::vector<value_type>& basis_vals) const;

    //! Print approximation in basis
    template <typename BasisT>
    std::ostream& print(const BasisT& basis, std::ostream& os) const;

    template <typename BasisT>
    value_type mean(const BasisT& basis) const;

    template <typename BasisT>
    value_type standard_deviation(const BasisT& basis) const;

    

  protected:

    //! OrthogPolyApprox coefficients
    std::vector<value_type> coeff_;
    
  }; // class OrthogPolyApprox

  template <typename ordinal_type, typename value_type> std::ostream& 
  operator << (std::ostream& os, 
	       const OrthogPolyApprox<ordinal_type,value_type>& a);

} // namespace Stokhos

#include "Stokhos_OrthogPolyApproxImp.hpp"

#endif //  STOKHOS_ORTHOGPOLYAPPROX_HPP
