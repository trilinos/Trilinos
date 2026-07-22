// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_VECTORORTHOGPOLY_HPP
#define STOKHOS_VECTORORTHOGPOLY_HPP

#include "Stokhos_ProductContainer.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  template <typename coeff_type>
  class VectorOrthogPoly : public virtual ProductContainer<coeff_type> {
  public:

    //! Typename of traits
    typedef typename ProductContainer<coeff_type>::traits_type traits_type;

    //! Typename of values
    typedef typename ProductContainer<coeff_type>::value_type value_type;

    //! Typename of ordinals
    typedef typename ProductContainer<coeff_type>::ordinal_type ordinal_type;

    //! Constructor with no basis
    /*!
     * Use with care!  Generally you will want to call reset() before using
     * any of the methods on this class.
     */
    VectorOrthogPoly();

    /*! 
     * \brief Create a polynomial for basis \c basis with empty 
     * coefficients
     */
    VectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& map);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated through a clone operation as implemented by the traits class
     * for the coefficient.
     */
    VectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& map,
      const typename traits_type::cloner_type& cloner);

    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    VectorOrthogPoly(const VectorOrthogPoly&);

    //! Destructor
    virtual ~VectorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    VectorOrthogPoly& operator=(const VectorOrthogPoly&);

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.
     */
    void reset(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& new_basis,
      const Teuchos::RCP<const Epetra_BlockMap>& new_map,
      const typename traits_type::cloner_type& cloner);

    //! Get basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > 
    basis() const;

    //! Get term for dimension \c dimension and order \c order
    coeff_type& term(ordinal_type dimension, ordinal_type order);

    //! Get term for dimension \c dimension and order \c order
    const coeff_type& term(ordinal_type dimension, ordinal_type order) const;

    //! Evaluate polynomial at supplied basis values
    /*!
     * Currently doesn't work with parallel map.
     */
    void evaluate(const Teuchos::Array<value_type>& basis_values, 
		  coeff_type& result) const;

    //! Evaluate polynomial at supplied basis values
    void sumIntoAllTerms(const value_type& weight,
			 const Teuchos::Array<value_type>& basis_values, 
			 const Teuchos::Array<value_type>& basis_norms,
			 const coeff_type& vec);

    //! Print polynomial
    std::ostream& print(std::ostream& os) const;

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis_;

  }; // class VectorOrthogPoly

  template <typename coeff_type>
  std::ostream& operator << (std::ostream& os, 
			     const VectorOrthogPoly<coeff_type>& vec) {
    return vec.print(os);
  }

} // end namespace Stokhos

// include template definitions
#include "Stokhos_VectorOrthogPolyImp.hpp"

#endif  // STOKHOS_VECTORORTHOGPOLY_HPP
