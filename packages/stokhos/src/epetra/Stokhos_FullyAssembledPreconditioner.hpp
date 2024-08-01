// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_FULLY_ASSEMBLED_PRECONDITIONER_HPP
#define STOKHOS_FULLY_ASSEMBLED_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"

#include "Stokhos_SGPreconditioner.hpp"
#include "Stokhos_AbstractPreconditionerFactory.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Stokhos {
    
  /*! 
   * \brief A stochastic preconditioner based on applying a preconditioner
   * to the fully assembled operator.
   */
  class FullyAssembledPreconditioner : public Stokhos::SGPreconditioner {
      
  public:

    //! Constructor 
    FullyAssembledPreconditioner(
      const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory,
      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    
    //! Destructor
    virtual ~FullyAssembledPreconditioner();

    /** \name Stokhos::SGPreconditioner methods */
    //@{

    //! Setup preconditioner
    virtual void 
    setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op, 
			const Epetra_Vector& x);

    //@}

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

  private:
    
    //! Private to prohibit copying
    FullyAssembledPreconditioner(const FullyAssembledPreconditioner&);
    
    //! Private to prohibit copying
    FullyAssembledPreconditioner& operator=(const FullyAssembledPreconditioner&);
    
  protected:
    
    //! Label for operator
    std::string label;

    //! Stores factory for building preconditioner
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory;

    //! Stores preconditioner
    Teuchos::RCP<Epetra_Operator> prec;

  }; // class FullyAssembledPreconditioner
  
} // namespace Stokhos

#endif // STOKHOS_FULLY_ASSEMBLED_PRECONDITIONER_HPP
