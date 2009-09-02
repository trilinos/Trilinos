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

#ifndef STOKHOS_VECTORORTHOGPOLY_HPP
#define STOKHOS_VECTORORTHOGPOLY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"

namespace Stokhos {

  //! Base traits definition for VectorOthogPoly
  template <typename coeff_type> class VectorOrthogPolyTraits {};

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  template <typename coeff_type>
  class VectorOrthogPoly {
  public:

    //! Typename of traits
    typedef Stokhos::VectorOrthogPolyTraits<coeff_type> traits_type;

    //! Typename of values
    typedef typename traits_type::value_type value_type;

    //! Typename of ordinals
    typedef typename traits_type::ordinal_type ordinal_type;

    /*! 
     * \brief Create a polynomial for basis \c basis with empty 
     * coefficients
     */
    VectorOrthogPoly(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
		     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad = Teuchos::null);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated through a clone operation as implemented by the traits class
     * for the coefficient.
     */
    VectorOrthogPoly(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
		     const typename traits_type::cloner_type& cloner,
		     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad = Teuchos::null);

    //! Destructor
    ~VectorOrthogPoly();

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
    ordinal_type size() const;

    //! Get basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > 
    basis() const;

    //! Get quadrature
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > 
    quadrature() const;

    //! Return ref-count pointer to coefficient \c i
    Teuchos::RCP<coeff_type>
    getCoeffPtr(ordinal_type i);

    //! Return ref-count pointer to constant coefficient \c i
    Teuchos::RCP<const coeff_type>
    getCoeffPtr(ordinal_type i) const;

    //! Set coefficient \c i to \c c
    void setCoeffPtr(ordinal_type i, const Teuchos::RCP<coeff_type>& c);

    //! Array access
    coeff_type& operator[](ordinal_type i);

    //! Array access
    const coeff_type& operator[](ordinal_type i) const;

    //! Get term for dimension \c dimension and order \c order
    coeff_type& term(ordinal_type dimension, ordinal_type order);

    //! Get term for dimension \c dimension and order \c order
    const coeff_type& term(ordinal_type dimension, ordinal_type order) const;

    //! Initialize polynomial coefficients
    void init(const value_type& val);

    //! Evaluate polynomial at supplied basis values
    void evaluate(const Teuchos::Array<value_type>& basis_values, 
		  coeff_type& result) const;

    //! Evaluate polynomial at supplied basis values
    void sumIntoAllTerms(const value_type& weight,
			 const Teuchos::Array<value_type>& basis_values, 
			 const Teuchos::Array<value_type>& basis_norms,
			 const coeff_type& vec);

  private:

    //! Prohibit copying
    VectorOrthogPoly(const VectorOrthogPoly&);

    //! Prohibit copying
    VectorOrthogPoly& operator=(const VectorOrthogPoly&);

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis_;

    //! Quadrature
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > quad_;

    //! Array of polynomial coefficients
    Teuchos::Array< Teuchos::RCP<coeff_type> > coeff_;

  }; // class VectorOrthogPoly

} // end namespace Stokhos

// include template definitions
#include "Stokhos_VectorOrthogPolyImp.hpp"

#endif  // STOKHOS_VECTORORTHOGPOLY_HPP
