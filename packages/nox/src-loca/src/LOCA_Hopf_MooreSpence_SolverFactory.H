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

#ifndef LOCA_HOPF_MOORESPENCE_SOLVERFACTORY_H
#define LOCA_HOPF_MOORESPENCE_SOLVERFACTORY_H

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
  namespace Hopf {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MooreSpence {

      /*! 
       * \brief Factory for creating solver objects for solving Moore-Spence 
       * Hopf equations
       */
      /*!
       * The parameters passed to the create() through the \c solverParams
       * argument method should specify the "Solver Method" 
       * as described below, as well as any additional parameters
       * for the particular strategy.
       * <ul>
       * <li> "Solver Method" - Name of the method. Valid choices are
       *   <ul>
       *   <li> "Salinger Bordering" (LOCA::Hopf::MooreSpence::SalingerBordering) [Default]
       *   </ul>
       * </ul>
       */
      class SolverFactory {

      public:

	//! Constructor
	SolverFactory(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data);

	//! Destructor
	virtual ~SolverFactory();

	//! Create solver strategy
	/*!
	 * \param topParams [in] Parsed top-level parameter list.  
	 * \param solverParams [in] Solver parameters as described above
	 */
	Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>
	create(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

	//! Return strategy name given by \c solverParams
	const std::string& strategyName(Teuchos::ParameterList& solverParams) const;

      private:

	//! Private to prohibit copying
	SolverFactory(const SolverFactory&);

	//! Private to prohibit copying
	SolverFactory& operator = (const SolverFactory&);
	
      protected:

	//! Global data
	Teuchos::RCP<LOCA::GlobalData> globalData;

      }; // Class SolverFactory
      
    } // Namespace MooreSpence

  } // Namespace Hopf
    
} // Namespace LOCA

#endif // LOCA_HOPF_MOORESPENCE_SOLVERFACTORY_H
