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
  template <typename T>
  class OrthogPolyApprox {
  public:

    //! Typename of values
    typedef T value_type;

    //! Default constructor
    OrthogPolyApprox();
    
    //! Constructor with supplied value \c x
    OrthogPolyApprox(const value_type& x);
    
    //! Constructor with size \c sz and value \c x
    OrthogPolyApprox(unsigned int sz, const value_type& x);
    
    //! Constructor with size \c sz
    OrthogPolyApprox(unsigned int sz);
    
    //! Constructor with basis, size \c sz and length \c l
    OrthogPolyApprox(unsigned int sz, unsigned int l);
    
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
    void resize(unsigned int sz);

    //! Reserve space for a size \c sz expansion
    /*!
     * Coefficients are preserved.
     */
    void reserve(unsigned int sz);

    //! Return size
    unsigned int size() const { return coeff_.size(); }

    //! Return coefficient array
    T* coeff() { return &coeff_[0]; }

    //! Return coefficient array
    const T* coeff() const { return &coeff_[0]; }

    //! Array access
    T& operator[](unsigned int i) { return coeff_[i]; }

    //! Array access
    const T& operator[](unsigned int i) const { return coeff_[i]; }

    //! Get term
    template <typename BasisT>
    T& 
    term(const BasisT& basis, 
	 int i0 = -1, int i1 = -1, int i2 = -1, int i3 = -1, int i4 = -1,
	 int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1);

    //! Get term
    template <typename BasisT>
    const T& 
    term(const BasisT& basis, 
	 int i0 = -1, int i1 = -1, int i2 = -1, int i3 = -1, int i4 = -1,
	 int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1) const;

    //! Write polynomial approximation in standard basis
    template <typename BasisT>
    Polynomial<T> toStandardBasis(const BasisT& basis) const;

    //! Evaluate polynomial approximation at a point
    template <typename BasisT>
    T evaluate(const BasisT& basis, const std::vector<T>& point) const;

    //! Print approximation in basis
    template <typename BasisT>
    std::ostream& print(const BasisT& basis, std::ostream& os) const;

  protected:

    //! OrthogPolyApprox coefficients
    std::vector<T> coeff_;
    
  }; // class OrthogPolyApprox

  template <typename T> std::ostream& 
  operator << (std::ostream& os, const OrthogPolyApprox<T>& a);

} // namespace Stokhos

#include "Stokhos_OrthogPolyApproxImp.hpp"

#endif //  STOKHOS_ORTHOGPOLYAPPROX_HPP
