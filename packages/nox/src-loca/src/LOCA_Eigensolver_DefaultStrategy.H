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

#ifndef LOCA_EIGENSOLVER_DEFAULTSTRATEGY_H
#define LOCA_EIGENSOLVER_DEFAULTSTRATEGY_H

#include "LOCA_Eigensolver_AbstractStrategy.H" // base class

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Eigensolver {

    //! Default eigensolver strategy
    /*!
     * This class implements a default eigensolver strategy that does
     * not compute any eigenvalues.
     */
    class DefaultStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      DefaultStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Destructor
      virtual ~DefaultStrategy();

      //! Compute eigenvalues/eigenvectors
      /*!
       * The implementation here does nothing and always returns
       * NOX::Abstract::Group::Ok.  Note that this implies the returned
       * ref-count pointers are null.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
          NOX::Abstract::Group& group,
          Teuchos::RCP< std::vector<double> >& evals_r,
          Teuchos::RCP< std::vector<double> >& evals_i,
          Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
          Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

    private:

      //! Private to prohibit copying
      DefaultStrategy(const DefaultStrategy&);

      //! Private to prohibit copying
      DefaultStrategy& operator = (const DefaultStrategy&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class DefaultStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
