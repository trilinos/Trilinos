// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_INTERFACE_TIMEDEPENDENTMATRIXFREE_H
#define LOCA_EPETRA_INTERFACE_TIMEDEPENDENTMATRIXFREE_H

#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Abstract_Group.H"

// Forward declarations
class Epetra_Vector;

namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for applying the shifted matrix in a matrix-free
       * setting.
       */
      /*!
       * This interface is derived from the
       * NOX::Epetra::Interface::Required and additionally provides a
       * method for applying the shifted matrix \f$\alpha J + \beta M\f$.
       * This is needed for linear stability analysis and Hopf tracking.
       */
      class TimeDependentMatrixFree :
    public virtual LOCA::Epetra::Interface::Required {

      public:

    //! Constructor
    TimeDependentMatrixFree() {};

    //! Destructor
    virtual ~TimeDependentMatrixFree() {};

    /*!
     * \brief Call user routine for applying the shifted matrix
     * \f$\alpha J + \beta M\f$ where \f$ J\f$ is the Jacobian matrix
     * and \f$ M\f$ is the mass matrix.
     */
    virtual bool
        applyShiftedMatrix(double alpha, double beta,
               const NOX::Epetra::Vector& input,
               NOX::Epetra::Vector& result) const = 0;

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
