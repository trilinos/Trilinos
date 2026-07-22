// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_MULTIVECTOR_OPERATOR_ORTHOG_POLY_HPP
#define STOKHOS_EPETRA_MULTIVECTOR_OPERATOR_ORTHOG_POLY_HPP

#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "Stokhos_ProductEpetraMultiVectorOperator.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  class EpetraMultiVectorOperatorOrthogPoly : 
    public EpetraOperatorOrthogPoly,
    public ProductEpetraMultiVectorOperator {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraMultiVectorOperatorOrthogPoly(
      const Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>& sg_mv,
      bool is_multi_vec_transposed);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraMultiVectorOperatorOrthogPoly(
      const EpetraMultiVectorOperatorOrthogPoly& v);

    //! Destructor
    virtual ~EpetraMultiVectorOperatorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraMultiVectorOperatorOrthogPoly& 
    operator=(const EpetraMultiVectorOperatorOrthogPoly& v);

    //! Get multi vector orthog poly
    Teuchos::RCP<EpetraMultiVectorOrthogPoly> 
    multiVectorOrthogPoly() const;

  protected:

    //! Multivector orthog poly
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> sg_mv;

  }; // class EpetraMultiVectorOperatorOrthogPoly

} // end namespace Stokhos

#endif  // STOKHOS_EPETRA_MULTIVECTOR_OPERATOR_ORTHOG_POLY_HPP
