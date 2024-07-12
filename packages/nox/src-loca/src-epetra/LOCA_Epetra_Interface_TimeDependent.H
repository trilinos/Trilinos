// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_INTERFACE_TIMEDEPENDENT_H
#define LOCA_EPETRA_INTERFACE_TIMEDEPENDENT_H

#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Abstract_Group.H"

// Forward declarations
class Epetra_Operator;

namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for computing the shifted matrix.
       */
      /*!
       * This interface is derived from the
       * NOX::Epetra::Interface::Jacobian and additionally provides a
       * method for computing the shifted matrix \f$\alpha J + \beta M\f$ .
       * This is needed for linear stability analysis and Hopf tracking.
       */
      class TimeDependent : public virtual LOCA::Epetra::Interface::Required {

      public:

    //! Constructor
    TimeDependent() {};

    //! Destructor
    virtual ~TimeDependent() {};

    /*!
     * \brief Call user routine for computing the shifted matrix
     * \f$\alpha J + \beta M\f$ where \f$ J\f$ is the Jacobian matrix
     * and \f$ M\f$ is the mass matrix.
     */
    virtual bool
        computeShiftedMatrix(double alpha, double beta,
                 const Epetra_Vector& x,
                 Epetra_Operator& A) = 0;

        /*!
          \brief Routine used in XYZT to set x_dot and time in the interface.

          The computeF() routine for XYZT problems needs to be a function
          of x_dot, but th NOX/LOCA computeF() does not take x_dot as an
          argument. This is used to set x_dot in the application interface
          so the subsequent call to computeF has the correct x_dot value.
          The  timeStep  argument is sent so the use can set the global time,
          for cases when computeF, computeJacobian, computeMassMatrix fills
          are functions of time (nonautonomous systems).
        */
        virtual void setXdot(const Epetra_Vector& /* xdot */, const double /* time */) {
      std::cout << "WARNING: "
            << "LOCA::Epetra::Interface::TimeDependent::setXdot"
            << "\n\tempty default implementation not overloaded!" << std::endl;
        };

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
