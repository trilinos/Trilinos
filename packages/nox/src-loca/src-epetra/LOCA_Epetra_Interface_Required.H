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
    virtual void printSolution(const Epetra_Vector& x_,
                   double conParam) {}

        /*!
          \brief Provides data to application for output files.

          This routine is called from Interface::xyzt::printSolution() just
          before the call to Interface::Required::printSolution(x,param),
          and gives the application some indices that can be used for
          creating a unique name/index for the output files.
        */
        virtual void dataForPrintSolution(const int conStep, const int timeStep,
                                          const int totalTimeSteps) {};

    /*!
      \brief Set multipoint parameter in the user's application.

      Should be called prior to calling one of the compute functions.
    */
        virtual void setMultiPointParameter(const int stepNum) {};


    //! Perform any preprocessing before a continuation step starts.
    /*!
     * The \c stepStatus argument indicates whether the previous step was
     * successful.  The default implementation here is empty.
     */
    virtual void
    preProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus,
                 LOCA::Epetra::Group& group) {}

    //! Perform any postprocessing after a continuation step finishes.
    /*!
     * The \c stepStatus argument indicates whether the step was
     * successful. The default implementation here is empty.
     */
    virtual void
    postProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus,
                 LOCA::Epetra::Group& group) {}

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
