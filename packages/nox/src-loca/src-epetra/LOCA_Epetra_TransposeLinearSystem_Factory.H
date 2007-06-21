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

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_FACTORY_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_FACTORY_H

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class LinearSystem;
  }
}
namespace LOCA {
  class GlobalData;
  namespace Epetra {
    namespace TransposeLinearSystem {
      class AbstractStrategy;
    }
  }
}

namespace LOCA {

  namespace Epetra {

    namespace TransposeLinearSystem {

      //! Factory for creating transpose linear system strategy objects
      /*!
       * The parameters passed to the create() through the \c solverParams
       * argument method should specify the "Transpose Solver Method" 
       * as described below, as well as any additional parameters
       * for the particular strategy.
       * <ul>
       * <li> "Transpose Solver Method" - Name of the method. Valid choices are
       *   <ul>
       *   <li> "Tranpose Preconditioner" 
       *        (NOX::Epetra::TransposeLinearSystem::TransposePreconditioner) 
       *        [Default]
       *   <li> "Explicit Transpose" 
       *        (NOX::Epetra::TransposeLinearSystem::ExplicitTranspose) 
       *   <li> "Left Preconditioning" 
       *        (NOX::Epetra::TransposeLinearSystem::LeftPreconditioning) 
       *   </ul>
       * </ul>
       */
      class Factory {

      public:

	//! Constructor
	Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

	//! Destructor
	virtual ~Factory();

	//! Create transpose solver strategy
	/*!
	 * \param solverParams [in] Solver parameters as described above
	 * \param linsys [in] Linear system solver
	 */
	Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy>
	create(
	    const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys);

	//! Return strategy name given by \c solverParams
	const string& strategyName(Teuchos::ParameterList& solverParams) const;

      private:

	//! Private to prohibit copying
	Factory(const Factory&);
	
	//! Private to prohibit copying
	Factory& operator = (const Factory&);
	
      protected:
	
	//! Global data object
	Teuchos::RCP<LOCA::GlobalData> globalData;
	
      }; // Class Factory

    } // NameSpace TransposeLinearSystem
      
  } // Namespace Epetra

} // Namespace LOCA

#endif // LOCA_EPETRA_TRANSPOSELINEARSYSTEM_FACTORY_H
