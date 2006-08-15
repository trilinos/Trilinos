// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_SOLVERFACTORY_H
#define LOCA_TURNINGPOINT_MOORESPENCE_SOLVERFACTORY_H

#include "Teuchos_RefCountPtr.hpp"

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace TurningPoint {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
}

namespace LOCA {

  namespace TurningPoint {

    namespace MooreSpence {

      /*! 
       * \brief Factory for creating solver objects for solving Moore-Spence 
       * turning point equations
       */
      /*!
       * The parameters passed to the create() through the \c solverParams
       * argument method should specify the "Solver Method" 
       * as described below, as well as any additional parameters
       * for the particular strategy.
       * <ul>
       * <li> "Solver Method" - Name of the method. Valid choices are
       *   <ul>
       *   <li> "Salinger Bordering" (LOCA::TurningPoint::MooreSpence::SalingerBordering) [Default]
       *   <li> "Phipps Bordering" (LOCA::TurningPoint::MooreSpence::PhippsBordering)
       *   </ul>
       * </ul>
       */
      class SolverFactory {

      public:

	//! Constructor
	SolverFactory(
		  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

	//! Destructor
	virtual ~SolverFactory();

	//! Create solver strategy
	/*!
	 * \param topParams [in] Parsed top-level parameter list.  
	 * \param solverParams [in] Solver parameters as described above
	 */
	Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::SolverStrategy>
	create(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams);

	//! Return strategy name given by \c solverParams
	const string& strategyName(Teuchos::ParameterList& solverParams) const;

      private:

	//! Private to prohibit copying
	SolverFactory(const SolverFactory&);

	//! Private to prohibit copying
	SolverFactory& operator = (const SolverFactory&);
	
      protected:

	//! Global data
	Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      }; // Class SolverFactory
      
    } // Namespace MooreSpence

  } // Namespace TurningPoint
    
} // Namespace LOCA

#endif // LOCA_TURNINGPOINT_MOORESPENCE_SOLVERFACTORY_H
