// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_INTERFACE_REQUIRED_H
#define LOCA_EPETRA_INTERFACE_REQUIRED_H

#include "NOX_Common.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Vector.H"
#include "LOCA_Abstract_Iterator.H"

// Forward declarations
class Epetra_Vector;

namespace LOCA {
  class ParameterVector;
  namespace Epetra {
    class Group;
  }
}

namespace LOCA {
  namespace Epetra {

    /*!
      \brief Provides a set of interfaces for users to provide information
      about the nonlinear problem to LOCA.

      Contains interfaces for the user to supply
      (1) the evaluation of the nonlinear equations,
      (2) any preconditioning if required,
      (3) space-time operators, and
      (4) shifted operators.
    */
    namespace Interface {

      /*!
    \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for setting problem parameters.

    This interface is derived from the NOX::Epetra::Interface::Required
    and additionally provides a method for setting problem parameters.
      */
      class Required : public virtual NOX::Epetra::Interface::Required {

      public:

    //! Constructor
    Required() {};

    //! Destructor
    virtual ~Required() {};

    /*!
      \brief Set parameters in the user's application.

      Should be called prior to calling one of the compute functions.
    */
    virtual void setParameters(const ParameterVector& p) = 0;

    //! Call user's own print routine for vector-parameter pair
    virtual void printSolution(const Epetra_Vector& /* x_ */,
                   double /* conParam */) {}

        /*!
          \brief Provides data to application for output files.

          This routine is called from Interface::xyzt::printSolution() just
          before the call to Interface::Required::printSolution(x,param),
          and gives the application some indices that can be used for
          creating a unique name/index for the output files.
        */
        virtual void dataForPrintSolution(const int /* conStep */, const int /* timeStep */,
                                          const int /* totalTimeSteps */) {};

    /*!
      \brief Set multipoint parameter in the user's application.

      Should be called prior to calling one of the compute functions.
    */
        virtual void setMultiPointParameter(const int /* stepNum */) {};


    //! Perform any preprocessing before a continuation step starts.
    /*!
     * The \c stepStatus argument indicates whether the previous step was
     * successful.  The default implementation here is empty.
     */
    virtual void
    preProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus /* stepStatus */,
                 LOCA::Epetra::Group& /* group */) {}

    //! Perform any postprocessing after a continuation step finishes.
    /*!
     * The \c stepStatus argument indicates whether the step was
     * successful. The default implementation here is empty.
     */
    virtual void
    postProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus /* stepStatus */,
                 LOCA::Epetra::Group& /* group */) {}

    /*!
      \brief Projects solution to a few scalars for
      multiparameter continuation

      Default implementation is the max norm.
    */
    virtual void projectToDraw(const NOX::Epetra::Vector& x,
                   double *px) const {
      px[0] = x.norm(NOX::Abstract::Vector::MaxNorm);
    }

    //! Returns the dimension of the projection to draw array
    virtual int projectToDrawDimension() const { return 1; }

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
