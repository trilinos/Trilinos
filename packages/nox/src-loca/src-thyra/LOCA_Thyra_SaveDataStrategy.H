// $Id$
// $Source$

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

#ifndef LOCA_THYRA_SAVEDATASTRATEGY_H
#define LOCA_THYRA_SAVEDATASTRATEGY_H

#include "Teuchos_RCP.hpp"
#include "NOX_Abstract_Vector.H"
#include "LOCA_Abstract_Iterator.H"

namespace LOCA {

  namespace Thyra {

    //! Abstract interface class strategies to save continuation data
    class SaveDataStrategy {

    public:

      //! Constructor
      SaveDataStrategy() {}

      //! Destructor
      virtual ~SaveDataStrategy() {}

      //! Save solution
      /*!
       * \param x [in] Solution vector
       * \param p [in] Parameter value
       * \returns ReturnType code indicating success or failure
       */
      virtual void
      saveSolution(const NOX::Abstract::Vector& /* x */, double /* p */) {}

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.  The default implementation to empty.
       */
      virtual void
      preProcessContinuationStep(
                LOCA::Abstract::Iterator::StepStatus /* stepStatus */) {}

      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.  The default implementation to empty.
       */
      virtual void
      postProcessContinuationStep(
                LOCA::Abstract::Iterator::StepStatus /* stepStatus */) {}

      //! Projects solution to a few scalars for multiparameter continuation
      /*!
       * This method is called every time a solution is saved by the
       * multiparameter continuation code MF for later visualization
       * and should project the solution vector down to a few scalars.
       * The array \c px will be preallocated to the proper length
       * given by projectToDrawDimension().
       */
      virtual void projectToDraw(const NOX::Abstract::Vector& /* x */,
                 double */* px */) const {}

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const { return 0; }

    private:

      //! Private to prohibit copying
      SaveDataStrategy(const SaveDataStrategy&);

      //! Private to prohibit copying
      SaveDataStrategy& operator = (const SaveDataStrategy&);

    }; // Class SaveDataStrategy

  } // Namespace Thyra

} // Namespace LOCA

#endif // LOCA_THYRA_SAVEDATASTRATEGY_H
