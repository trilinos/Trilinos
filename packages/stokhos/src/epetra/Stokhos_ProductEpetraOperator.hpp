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

#ifndef STOKHOS_PRODUCT_EPETRA_OPERATOR_HPP
#define STOKHOS_PRODUCT_EPETRA_OPERATOR_HPP

#include "Stokhos_ProductContainer.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Epetra_Map.h"
#include "EpetraExt_MultiComm.h"

namespace Stokhos {

  /*! 
   * \brief A container class for products of Epetra_Vector's.  
   */
  class ProductEpetraOperator :
    public virtual ProductContainer<Epetra_Operator>,
    public virtual Epetra_Operator {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    /*! 
     * \brief Create a container with container map \c block_map where each 
     * coefficient is generated from the supplied coefficient map \c coeff_map
     */
    ProductEpetraOperator(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    /*! 
     * \brief Create a container with container map \c block_map where each 
     * coefficient is generated from the supplied coefficient map \c coeff_map
     */
    /*
     * This version supplies the generated product map \c product_map
     */
    ProductEpetraOperator(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductEpetraOperator(const ProductEpetraOperator& v);

    //! Destructor
    virtual ~ProductEpetraOperator();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    ProductEpetraOperator& operator=(const ProductEpetraOperator& v);

    //! Get product comm
    Teuchos::RCP<const EpetraExt::MultiComm> productComm() const;

    /** \name Epetra_Operator methods */
    //@{
    
    //! Set to true if the transpose of the operator is requested
    virtual int SetUseTranspose(bool UseTranspose);
    
    /*! 
     * \brief Returns the result of a Epetra_Operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int Apply(const Epetra_MultiVector& Input, 
                      Epetra_MultiVector& Result) const;

    /*! 
     * \brief Returns the result of the inverse of the operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;
    
    //! Returns an approximate infinity norm of the operator matrix.
    virtual double NormInf() const;
    
    //! Returns a character string describing the operator
    virtual const char* Label () const;
  
    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;
    
    /*! 
     * \brief Returns true if the \e this object can provide an 
     * approximate Inf-norm, false otherwise.
     */
    virtual bool HasNormInf() const;

    /*! 
     * \brief Returns a reference to the Epetra_Comm communicator 
     * associated with this operator.
     */
    virtual const Epetra_Comm & Comm() const;

    /*!
     * \brief Returns the Epetra_Map object associated with the 
     * domain of this matrix operator.
     */
    virtual const Epetra_Map& OperatorDomainMap () const;

    /*! 
     * \brief Returns the Epetra_Map object associated with the 
     * range of this matrix operator.
     */
    virtual const Epetra_Map& OperatorRangeMap () const;

    //@}

  protected:

    //! Protected constructor to allow 2-stage derived setup
    ProductEpetraOperator(
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    //! Second stage of setup
    void setup(const Teuchos::RCP<const Epetra_Map>& domain_base_map,
	       const Teuchos::RCP<const Epetra_Map>& range_base_map);

  protected:

    //! Domain map of each coefficient
    Teuchos::RCP<const Epetra_Map> domain_base_map;

    //! Range map of each coefficient
    Teuchos::RCP<const Epetra_Map> range_base_map;

    //! Product range map
    Teuchos::RCP<const Epetra_Map> product_range_map;

    //! Product multi-level communicator
    Teuchos::RCP<const EpetraExt::MultiComm> product_comm;

    //! Whether to use transpose in Apply()
    bool useTranspose;

  }; // class ProductEpetraOperator

} // end namespace Stokhos

#endif  // STOKHOS_PRODUCT_EPETRA_OPERATOR_HPP
