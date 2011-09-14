//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_TIMED_EPETRA_OP_HPP
#define EPETRAEXT_TIMED_EPETRA_OP_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include <Teuchos_Time.hpp>
#include <EpetraExt_ConfigDefs.h>

namespace EpetraExt {
    
  /*! 
   * \brief Class allows for timing the action and inverse action of an
   * Epetra_Opetator.
   */
  class Epetra_Timed_Operator : public Epetra_Operator {
      
  public:

    //! Constructor 
    Epetra_Timed_Operator(const Teuchos::RCP<Epetra_Operator>& A_);
    
    //! Destructor
    virtual ~Epetra_Timed_Operator();

    //! Set to true if the transpose of the operator is requested
    int SetUseTranspose(bool useTranspose);
    
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
    
    /*!
     * \brief Returns the total time applying this operator.
     */
    virtual double ApplyTime() const{return this->ApplyTimer->totalElapsedTime(false);};

    /*!
     * \brief Returns the total time applying the inverse of this operator.
     */
    virtual double ApplyInverseTime() const{return this->ApplyInverseTimer->totalElapsedTime(false);};

    /*!
     * \brief Returns a pointer to the underlying Epetra_Operator
     */
    virtual Teuchos::RCP<const Epetra_Operator> ReturnOperator() const{return this->A;};

  private:
    
    //! Private to prohibit copying
    Epetra_Timed_Operator(const Epetra_Timed_Operator&);
    
    //! Private to prohibit copying
    Epetra_Timed_Operator & operator=(const Epetra_Timed_Operator&);
    
  protected:
    
    //! Stores the base operator
    Teuchos::RCP<Epetra_Operator> A;

    //! Keeps track of the apply time 
    Teuchos::RCP<Teuchos::Time> ApplyTimer;

    //! Keeps track of the apply inverse time
    Teuchos::RCP<Teuchos::Time> ApplyInverseTimer;

  }; // class Epetra_Timed_Operator
  
} // namespace EpetraExt

#endif // EPETRAEXT_TIMED_EPETRA_OP_HPP
