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

#ifndef LOCA_BORDEREDSOLVER_FACTORY_H
#define LOCA_BORDEREDSOLVER_FACTORY_H

#include "Teuchos_RCP.hpp"

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace BorderedSolver {

    //! Factory for creating %BorderedSolver strategy objects
    /*!
     * The parameters passed to the create() through the \c solverParams
     * argument method should specify the "Bordered Solver Method" 
     * as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Bordered Solver Method" - Name of the method. Valid choices are
     *   <ul>
     *   <li> "Bordering" (LOCA::BorderedSolver::Bordering) [Default]
     *   <li> "Nested" (LOCA::BorderedSolver::Nested)
     *   </ul>
     * </ul>
     *
     * There are also %Epetra and %LAPACK specific strategies that can
     * be instantiated by the LOCA::Epetra::Factory and LOCA::LAPACK::Factory.
     * See LOCA::BorderedSolver::LAPACKDirectSolve,
     * LOCA::BorderedSolver::EpetraHouseholder and
     * LOCA::BorderedSolver::Epetra::Augmented.
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create bordered system solver strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.  
       * \param solverParams [in] Solver parameters as described above
       */
      Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>
      create(
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Return strategy name given by \c solverParams
      const std::string& strategyName(Teuchos::ParameterList& solverParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_FACTORY_H
