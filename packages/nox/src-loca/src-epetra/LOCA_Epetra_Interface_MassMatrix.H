// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_INTERFACE_MASSMATRIX_H
#define LOCA_EPETRA_INTERFACE_MASSMATRIX_H

#include "NOX_Common.H"
#include "NOX_Abstract_Group.H"

// Forward declarations
class Epetra_Vector;

namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
    \brief
    Used by LOCA::Epetra::Group to provide a link to the external code
    for the MassMatrix (coefficients of time dependent terms).

    This is used for Hopf bifurcation tracking, linear stability analysis,
    and space-time solutions (xyzt).
      */
      class MassMatrix  {

      public:

    //! Constructor
    MassMatrix() {};

    //! Destructor
    virtual ~MassMatrix() {};

        /*!
      \brief
      Compute MassMatrix given the specified input vector x.
      Returns true if computation was successful.
    */
        virtual bool computeMassMatrix(const Epetra_Vector& x) = 0;

        /*!
      \brief Routines used in XYZT to set the old solution,
      the one from the previous time step.

      There is a different routine for first step,
      where the old solution is set by the user and not part of the
      solution vector. These routines are used by space-time (xyzt)
      problems, where the residual vector is a function of the previous
      solution, which is also being solved for, and where the MassMatrix
      is calculated as a function of a different solution vector then
      the Jacobian (that is, the previous time step).
      The  timeStep  argument is sent so the use can set the global time,
      for cases when computeF, computeJacobian, computeMassMatrix fills
      are functions of time (nonautonomous systems).
    */
        virtual void setOldSolution(const Epetra_Vector& x,const int timeStep) {
          std::cout << "WARNING: "
           << "LOCA::Epetra::Interface::MassMatrix::setOldSolution"
           << "\n\tempty default implementation not overloaded!" << std::endl;
    };

    //! See setOldSolution description.
        virtual void setOldSolutionFirstStep() {
          std::cout << "WARNING: "
           << "LOCA::Epetra::Interface::MassMatrix::setOldSolutionFirstStep"
           << "\n\tempty default implementation not overloaded!" << std::endl;
    };

        /*!
      \brief Provides data to application for output files.

      This routine is called from Interface::xyzt::printSolution() just
      before the call to Interface::Required::printSolution(x,param),
      and gives the application some indices that can be used for
      creating a unique name/index for the output files.
    */
        virtual void dataForPrintSolution(const int conStep, const int timeStep,
                                          const int totalTimeSteps) {};

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
