// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_PRODUCT_EPETRA_MULTIVECTOR_OPERATOR_HPP
#define STOKHOS_PRODUCT_EPETRA_MULTIVECTOR_OPERATOR_HPP

#include "Stokhos_ProductEpetraOperator.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class for products of Epetra_Vector's.  
   */
  class ProductEpetraMultiVectorOperator :
    public virtual ProductEpetraOperator {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    /*! 
     * \brief Create ProductEpetraOperator out of ProductEpetraMultiVector
     */
    ProductEpetraMultiVectorOperator(
      const Teuchos::RCP<ProductEpetraMultiVector>& product_mv,
      bool is_multi_vec_transposed);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductEpetraMultiVectorOperator(const ProductEpetraMultiVectorOperator& v);

    //! Destructor
    virtual ~ProductEpetraMultiVectorOperator();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductEpetraMultiVectorOperator& 
    operator=(const ProductEpetraMultiVectorOperator& v);

    //! Get product multi vector
    Teuchos::RCP<ProductEpetraMultiVector> productMultiVector() const;

  protected:

    //! The product multi-vector
    Teuchos::RCP<ProductEpetraMultiVector> product_mv;

  }; // class ProductEpetraMultiVectorOperator

} // end namespace Stokhos

#endif  // STOKHOS_PRODUCT_EPETRA_MULTIVECTOR_OPERATOR_HPP
