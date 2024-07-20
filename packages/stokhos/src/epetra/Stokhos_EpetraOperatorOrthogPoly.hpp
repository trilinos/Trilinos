// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_OPERATOR_ORTHOG_POLY_HPP
#define STOKHOS_EPETRA_OPERATOR_ORTHOG_POLY_HPP

#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_ProductEpetraOperator.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  class EpetraOperatorOrthogPoly : 
    public virtual VectorOrthogPoly<Epetra_Operator>,
    public virtual ProductEpetraOperator {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraOperatorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by a created block vector
     */
    EpetraOperatorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraOperatorOrthogPoly(const EpetraOperatorOrthogPoly& v);

    //! Destructor
    virtual ~EpetraOperatorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraOperatorOrthogPoly& operator=(const EpetraOperatorOrthogPoly& v);

  protected:

    //! Protected constructor to allow 2-stage derived setup
    EpetraOperatorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    //! Second stage of setup
    void setup(const Teuchos::RCP<const Epetra_Map>& domain_base_map,
	       const Teuchos::RCP<const Epetra_Map>& range_base_map);

  }; // class EpetraOperatorOrthogPoly

} // end namespace Stokhos

#endif  // STOKHOS_EPETRA_OPERATOR_ORTHOG_POLY_HPP
