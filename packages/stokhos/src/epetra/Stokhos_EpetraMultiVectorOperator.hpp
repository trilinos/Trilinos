// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_MULTIVECTOR_OPERATOR_HPP
#define STOKHOS_EPETRA_MULTIVECTOR_OPERATOR_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"

namespace Stokhos {
    
  /*! 
   * \brief An adaptor that supplies the operator interface to a multi-vector
   */
  /*!
   * Supplies the operator interface for the operator C = op(A)*B
   * where A, B, and C are multi-vectors, and op() is the identity or 
   * transpose operator.  The convention is when the transpose is not
   * requested, A is treated as a row-oriented multi-vector and therefore
   * the domain map is the map of A.  The corresponding range map is a local 
   * map of size given by the number of vectors  When the transpose
   * is selected, these are reversed.
   *
   * This should really live somewhere else as it has nothing to do with
   * Stokhos.
   */
  class EpetraMultiVectorOperator : public Epetra_Operator {
      
  public:

    //! Constructor 
    EpetraMultiVectorOperator(
      const Teuchos::RCP<const Epetra_MultiVector>& multi_vec,
      bool is_multi_vec_transposed);

    //! Constructor 
    EpetraMultiVectorOperator(
      const Teuchos::RCP<Epetra_MultiVector>& multi_vec,
      bool is_multi_vec_transposed);
    
    //! Destructor
    virtual ~EpetraMultiVectorOperator();

    Teuchos::RCP<const Epetra_MultiVector>
    getMultiVector() const { return multi_vec; }

    Teuchos::RCP<Epetra_MultiVector>
    getMultiVector() { return nonconst_multi_vec; }
    
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

  private:
    
    //! Private to prohibit copying
    EpetraMultiVectorOperator(const EpetraMultiVectorOperator&);
    
    //! Private to prohibit copying
    EpetraMultiVectorOperator& operator=(const EpetraMultiVectorOperator&);
    
  protected:
    
    //! Label for operator
    std::string label;
    
    //! Multi-vector
    Teuchos::RCP<const Epetra_MultiVector> multi_vec;

    //! Non-const multi-vector
    Teuchos::RCP<Epetra_MultiVector> nonconst_multi_vec;
    
    //! Whether the multivector is already transposed
    bool is_multi_vec_transposed;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Domain map ( = number of columns of multi_vec)
    Teuchos::RCP<Epetra_Map> domain_map;

  }; // class EpetraMultiVectorOperator
  
} // namespace Stokhos

#endif // STOKHOS_EPETRA_MULTIVECTOR_OPERATOR_HPP
