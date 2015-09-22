// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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
