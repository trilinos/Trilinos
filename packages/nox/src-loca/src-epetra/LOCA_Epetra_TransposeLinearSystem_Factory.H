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
	const std::string& strategyName(Teuchos::ParameterList& solverParams) const;

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
