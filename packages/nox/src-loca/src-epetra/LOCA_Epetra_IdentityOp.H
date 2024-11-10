// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_IDENTITYOP_H
#define LOCA_EPETRA_IDENTITYOP_H

#include "Teuchos_RCP.hpp"

#include "Epetra_Operator.h"    // base class
#include "NOX_Common.H"         // for std::string data member

namespace LOCA {

  namespace Epetra {

    /*!
     * \brief An Epetra operator representing the identity matrix
     */
    class IdentityOp : public Epetra_Operator {

    public:

      //! Constructor
      /*!
       * \param comm [in] Comm object
       * \param map [in] Map object
       */
      IdentityOp(const Teuchos::RCP<const Epetra_Comm>& comm,
         const Teuchos::RCP<const Epetra_Map>& map);

      //! Destructor
      virtual ~IdentityOp();

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

    private:

      //! Private to prohibit copying
      IdentityOp(const IdentityOp&);

      //! Private to prohibit copying
      IdentityOp& operator=(const IdentityOp&);

    protected:

      //! Label for operator
      std::string label;

      //! Stores comm
      Teuchos::RCP<const Epetra_Comm> comm;

      //! Stores map
      Teuchos::RCP<const Epetra_Map> map;

      //! Flag indicating whether transpose was selected
      bool useTranspose;

    };
  } // namespace Epetra
} // namespace LOCA

#endif
