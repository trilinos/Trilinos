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

#ifndef LOCA_ANASAZIOPERATOR_FACTORY_H
#define LOCA_ANASAZIOPERATOR_FACTORY_H

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
  namespace AnasaziOperator {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace AnasaziOperator {

    //! Factory for creating Anasazi operator strategy objects
    /*!
     * The parameters passed to the create() through the \c eigenParams
     * argument method should specify the
     * "Operator" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Operator" - Name of the Anasazi operator. Valid choices are
     *   <ul>
     *   <li> "Jacobian Inverse" (LOCA::AnasaziOperator::JacobianInverse) 
     *        [Default]
     *   <li> "Shift-Invert" (LOCA::AnasaziOperator::ShiftInvert)
     *   <li> "Cayley" (LOCA::AnasaziOperator::Cayley)
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create Anasazi operator strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.  
       * \param eigenParams [in] %Eigensolver parameters as described above
       * \param solverParams [in] Linear solver parameters
       * \param grp [in] Group representing Jacobian/mass matrices
       */
      Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy>
      create(
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
	const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp);

      //! Return strategy name given by \c eigenParams
      const string& strategyName(Teuchos::ParameterList& eigenParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace AnasaziOperator

} // Namespace LOCA

#endif // LOCA_ANASAZIOPERATOR_FACTORY_H
