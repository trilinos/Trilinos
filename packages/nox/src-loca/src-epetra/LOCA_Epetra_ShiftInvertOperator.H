// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
