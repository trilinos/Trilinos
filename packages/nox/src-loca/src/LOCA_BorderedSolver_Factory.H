// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSOLVER_FACTORY_H
#define LOCA_BORDEREDSOLVER_FACTORY_H

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
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create bordered system solver strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.  
       * \param solverParams [in] Solver parameters as described above
       */
      Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy>
      create(
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams);

      //! Return strategy name given by \c solverParams
      const string& strategyName(Teuchos::ParameterList& solverParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_FACTORY_H
