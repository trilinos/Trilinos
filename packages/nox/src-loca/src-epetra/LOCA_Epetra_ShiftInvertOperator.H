
//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_SHIFTINVERTOPERATOR_H
#define LOCA_EPETRA_SHIFTINVERTOPERATOR_H

#include "Teuchos_RCP.hpp"

#include "Epetra_Operator.h"
#include "NOX_Epetra_Interface_Jacobian.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Epetra {
    class Group;
  }
}

namespace LOCA {

  namespace Epetra {

    //! Interface for LOCA::Epetra::ShifterInvertOperator
    class ShiftInvertInterface : public NOX::Epetra::Interface::Jacobian {

    public:

      //! Constructor
      ShiftInvertInterface() {};

      //! Destructor
      ~ShiftInvertInterface() {};

      //! Compute Jacobian \f$J\f$.
      bool computeJacobian(const Epetra_Vector &/* x */,
               Epetra_Operator& /* Jac */) {return true;};
    };

    //! Epetra operator for \f$(J-\sigma M)^{-1}\f$.
    class ShiftInvertOperator : public Epetra_Operator {

    public:

      //! Constructor
      ShiftInvertOperator(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const Teuchos::RCP<LOCA::Epetra::Group>& grp,
            const Teuchos::RCP<const Epetra_Operator>& jac,
            double shift);

      //! Destructor
      virtual ~ShiftInvertOperator();

      //! Set transpose
      virtual int SetUseTranspose(bool UseTranspose);

      //! Apply shifted operator
      virtual int Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y) const;

      //! Apply shifted operator inverse
      virtual int ApplyInverse(const Epetra_MultiVector& X,
                   Epetra_MultiVector& Y) const;

      //! Computing infinity norm
      virtual double NormInf() const;

      //! Label
      virtual const char * Label() const;

      //! Transpose
      virtual bool UseTranspose() const;

      //! Have norm-inf
      virtual bool HasNormInf() const;

      /*!
       * \brief Returns a pointer to the Epetra_Comm communicator associated
       * with this operator.
       */
      virtual const Epetra_Comm & Comm() const;

      /*!
       * \brief Returns the Epetra_Map object associated with the domain of
       * this operator.
       */
      virtual const Epetra_Map & OperatorDomainMap() const;

      /*!
       * \brief Returns the Epetra_Map object associated with the range of
       * this operator.
       */
      virtual const Epetra_Map & OperatorRangeMap() const;

    private:

      //! Private to prohibit copying
      ShiftInvertOperator(const ShiftInvertOperator&);

      //! Private to prohibit copying
      ShiftInvertOperator& operator=(const ShiftInvertOperator&);

    protected:

      //! Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! LOCA group
      Teuchos::RCP<const LOCA::Epetra::Group> locagrp;

      //! Jacobian operator
      Teuchos::RCP<const Epetra_Operator> jacOper;

      //! Shift value
      double shift_;

      //! Label
      const char * Label_;
    };
  } // namespace Epetra
} // namespace LOCA
#endif
