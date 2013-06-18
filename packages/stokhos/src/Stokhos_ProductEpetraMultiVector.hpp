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

#ifndef STOKHOS_PRODUCT_EPETRA_MULTI_VECTOR_HPP
#define STOKHOS_PRODUCT_EPETRA_MULTI_VECTOR_HPP

#include "Stokhos_ProductContainer.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Epetra_BlockMap.h"
#include "EpetraExt_MultiComm.h"
#include "EpetraExt_BlockMultiVector.h"

namespace Stokhos {

  /*! 
   * \brief A container class storing products of Epetra_MultiVector's.  
   */
  class ProductEpetraMultiVector : 
    public virtual ProductContainer<Epetra_MultiVector> {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Default constructor
    /*!
     * Use with care!  Generally you will want to call reset() before using
     * any of the methods on this class.
     */
    ProductEpetraMultiVector();

    /*! 
     * \brief Create a container with container map \c block_map
     */
    ProductEpetraMultiVector(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map);

    /*! 
     * \brief Create a container with container map \c block_map where each 
     * coefficient is generated from the supplied coefficient map \c coeff_map
     */

    ProductEpetraMultiVector(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    /*! 
     * \brief Create a container with container map \c block_map where each 
     * coefficient is generated from the supplied coefficient map \c coeff_map
     */
    /*
     * This version supplies the generated product map \c product_map
     */
    ProductEpetraMultiVector(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    /*! 
     * \brief Create a container with container map \c block_map where each 
     * coefficient is given by the supplied block vector.
     */
    ProductEpetraMultiVector(
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
    ProductEpetraMultiVector(const ProductEpetraMultiVector& v);

    //! Destructor
    virtual ~ProductEpetraMultiVector();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductEpetraMultiVector& operator=(const ProductEpetraMultiVector& v);

    //! Assignment
    ProductEpetraMultiVector& operator=(const Epetra_MultiVector& v);

    //! Assignment
    void assignToBlockMultiVector(Epetra_MultiVector& v) const;

    //! Assignment
    void assignFromBlockMultiVector(const Epetra_MultiVector& v);
      
    //! Get coefficient map
    Teuchos::RCP<const Epetra_BlockMap> coefficientMap() const;

    //! Get product map
    Teuchos::RCP<const Epetra_BlockMap> productMap() const;

    //! Get product comm
    Teuchos::RCP<const EpetraExt::MultiComm> productComm() const;

    //! Get number of vectors
    int numVectors() const;

    //! Reset to a new size
    /*!
     * This resizes array to fit new size.
     */
    void reset(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    //! Reset to a new size
    /*!
     * This resizes array to fit new size.
     */
    void reset(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      int num_vectors);

    //! Reset vector cofficients
    void resetCoefficients(Epetra_DataAccess CV,
			   const Epetra_MultiVector& block_vector);

    //! Get block vector
    Teuchos::RCP<EpetraExt::BlockMultiVector> getBlockMultiVector();

    //! Get block vector
    Teuchos::RCP<const EpetraExt::BlockMultiVector> getBlockMultiVector() const;

    //! Set block vector
    void setBlockMultiVector(
      const Teuchos::RCP<EpetraExt::BlockMultiVector>& block_vec);

  protected:

    //! Product map of block vector
    Teuchos::RCP<const Epetra_BlockMap> coeff_map;

    //! Product multi-level communicator
    Teuchos::RCP<const EpetraExt::MultiComm> product_comm;

    //! Product map of block vector
    Teuchos::RCP<const Epetra_BlockMap> product_map;

    //! Block vector storing coefficients
    Teuchos::RCP<EpetraExt::BlockMultiVector> bv;

  }; // class ProductEpetraMultiVector

} // end namespace Stokhos

#endif  // STOKHOS_PRODUCT_EPETRA_MULTI_VECTOR_HPP
