// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_PCE_COVARIANCE_OP_HPP
#define STOKHOS_PCE_COVARIANCE_OP_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Epetra_Operator.h"
#include "EpetraExt_BlockVector.h"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the covariance operator of a
   * polynomial chaos expansion.
   */
  /*!
   * If X is the matrix whose columns are the coefficients of a given
   * polynomial chaos expansion, starting at order 1 (not including mean term),
   * and S is a diagonal matrix whose entries are given by the norm-squared
   * of the basis polynomials, then the covariance operator is A = X*S*X^T.
   */
  class PCECovarianceOp : public Epetra_Operator {
  public:

    //! Constructor with polynomial X
    PCECovarianceOp(const Stokhos::VectorOrthogPoly<Epetra_Vector>& X_poly);

    //! Constructor with block-vector X
    PCECovarianceOp(const Teuchos::RCP<const EpetraExt::BlockVector>& X,
		    const Stokhos::OrthogPolyBasis<int,double>& basis);

    //! Constructor with multi-vector X
    PCECovarianceOp(const Teuchos::RCP<const Epetra_MultiVector>& X,
		    const Stokhos::OrthogPolyBasis<int,double>& basis);
    
    //! Destructor
    virtual ~PCECovarianceOp();
    
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
    
    //! Returns a character std::string describing the operator
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

    //! Returns PCE coefficient map
    const Epetra_BlockMap& CoeffMap() const;

  private:
    
    //! Private to prohibit copying
    PCECovarianceOp(const PCECovarianceOp&);
    
    //! Private to prohibit copying
    PCECovarianceOp& operator=(const PCECovarianceOp&);
    
  protected:
    
    //! Label for operator
    std::string label;

    //! Multivector X defining A = X*S*X^T
    Teuchos::RCP<const Epetra_MultiVector> X;

    //! Scaling vector in A = X*S*X^T
    Teuchos::Array<double> s;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Map needed for temporary vector
    Teuchos::RCP<Epetra_Map> tmp_map;

    //! Temporary vector needed for apply
    mutable Teuchos::RCP<Epetra_MultiVector> tmp;

  }; // class PCECovarianceOp
  
} // namespace Stokhos

#endif // STOKHOS_PCE_COVARIANCE_OP_HPP
