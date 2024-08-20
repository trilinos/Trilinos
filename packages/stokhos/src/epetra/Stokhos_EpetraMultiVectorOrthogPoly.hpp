// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EPETRAMULTIVECTORORTHOGPOLY_HPP
#define STOKHOS_EPETRAMULTIVECTORORTHOGPOLY_HPP

#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  class EpetraMultiVectorOrthogPoly : 
    public VectorOrthogPoly<Epetra_MultiVector>,
    public ProductEpetraMultiVector {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Constructor with no basis
    /*!
     * Use with care!  Generally you will want to call reset() before using
     * any of the methods on this class.
     */
    EpetraMultiVectorOrthogPoly();

    /*! 
     * \brief Create a polynomial for basis \c basis with empty 
     * coefficients
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by a created block vector
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by the supplied block vector.
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      Epetra_DataAccess CV,
      const Epetra_MultiVector& block_vector);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraMultiVectorOrthogPoly(const EpetraMultiVectorOrthogPoly& v);

    //! Destructor
    virtual ~EpetraMultiVectorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraMultiVectorOrthogPoly& operator=(const EpetraMultiVectorOrthogPoly& v);

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.
     */
    void reset(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.
     */
    void reset(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    //! Compute mean
    void computeMean(Epetra_MultiVector& v) const;

    //! Compute variance
    void computeVariance(Epetra_MultiVector& v) const;

    //! Compute standard deviation
    void computeStandardDeviation(Epetra_MultiVector& v) const;

  }; // class EpetraMultiVectorOrthogPoly

} // end namespace Stokhos

#endif  // STOKHOS_EPETRAMULTIVECTORORTHOGPOLY_HPP
