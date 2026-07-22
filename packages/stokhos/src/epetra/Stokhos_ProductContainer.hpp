// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_PRODUCT_CONTAINER_HPP
#define STOKHOS_PRODUCT_CONTAINER_HPP

#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Epetra_BlockMap.h"

namespace Stokhos {

  //! Base traits definition for ProductContainer
  template <typename coeff_type> class ProductContainerTraits {};

  /*! 
   * \brief A product (in the mathematical sense) container class whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  template <typename coeff_type>
  class ProductContainer {
  public:

    //! Typename of traits
    typedef Stokhos::ProductContainerTraits<coeff_type> traits_type;

    //! Typename of values
    typedef typename traits_type::value_type value_type;

    //! Typename of ordinals
    typedef typename traits_type::ordinal_type ordinal_type;

    //! Default constructor
    /*!
     * Use with care!  Generally you will want to call reset() before using
     * any of the methods on this class.
     */
    ProductContainer();

    /*! 
     * \brief Create a container with container map \c map
     */
    ProductContainer(const Teuchos::RCP<const Epetra_BlockMap>& map);

    /*! 
     * \brief Create a container container map \c map where each coefficient is 
     * generated through a clone operation as implemented by the traits class
     * for the coefficient.
     */
    ProductContainer(const Teuchos::RCP<const Epetra_BlockMap>& map,
		     const typename traits_type::cloner_type& cloner);

    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductContainer(const ProductContainer&);

    //! Destructor
    virtual ~ProductContainer();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductContainer& operator=(const ProductContainer&);

    //! Resize to new map \c map
    /*!
     * Any previous coefficients are lost.
     */
    void reset(const Teuchos::RCP<const Epetra_BlockMap>& map);

    //! Resize to new map \c map and create coefficients from \c cloner
    /*!
     * Any previous coefficients are lost.
     */
    void reset(const Teuchos::RCP<const Epetra_BlockMap>& map,
	       const typename traits_type::cloner_type& cloner);

    //! Resize to map \c map
    /*!
     * Coefficients are preserved.
     */
    void resize(const Teuchos::RCP<const Epetra_BlockMap>& map);

    //! Reserve space for a size \c sz container
    /*!
     * Coefficients are preserved.
     */
    void reserve(ordinal_type sz);

    //! Return size
    ordinal_type size() const;

    //! Return container map
    Teuchos::RCP<const Epetra_BlockMap> map() const;

    //! Return array of coefficients
    const Teuchos::Array<Teuchos::RCP<coeff_type> >&
    getCoefficients() const;

    //! Return array of coefficients
    Teuchos::Array<Teuchos::RCP<coeff_type> >&
    getCoefficients();

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

    //! Initialize coefficients
    void init(const value_type& val);

    //! Return whether global index i resides on this processor
    bool myGID(int i) const;

    //! Print polynomial
    std::ostream& print(std::ostream& os) const;

  protected:

    //! Container map
    Teuchos::RCP<const Epetra_BlockMap> map_;

    //! Array of polynomial coefficients
    Teuchos::Array< Teuchos::RCP<coeff_type> > coeff_;

  }; // class ProductContainer

  template <typename coeff_type>
  std::ostream& operator << (std::ostream& os, 
			     const ProductContainer<coeff_type>& vec) {
    return vec.print(os);
  }

} // end namespace Stokhos

// include template definitions
#include "Stokhos_ProductContainerImp.hpp"

#endif  // STOKHOS_PRODUCT_VECTOR_HPP
