// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef USER_APP_NOX_OBSERVER_FACTORY_HPP
#define USER_APP_NOX_OBSERVER_FACTORY_HPP

#include "Panzer_config.hpp"

#include "Panzer_STK_NOXObserverFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_Traits.hpp"
#include "NOX_PrePostOperator_Vector.H"

// Individual Observers
#include "user_app_NOXObserver_EpetraToExodus.hpp"

namespace user_app {
  
  class NOXObserverFactory_Epetra : public panzer_stk::NOXObserverFactory {
    
  public:
    NOXObserverFactory_Epetra(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary)
       : stkIOResponseLibrary_(stkIOResponseLibrary) {}
    
    Teuchos::RCP<NOX::Abstract::PrePostOperator>
    buildNOXObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
		     const Teuchos::RCP<panzer::UniqueGlobalIndexerBase>& dof_manager,
		     const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      Teuchos::RCP<NOX::PrePostOperatorVector> observer = 
	Teuchos::rcp(new NOX::PrePostOperatorVector);

      // Always register the exodus writer to output solution
      {
	Teuchos::RCP<NOX::Abstract::PrePostOperator> solution_writer = 
	  Teuchos::rcp(new user_app::NOXObserver_EpetraToExodus(mesh,dof_manager,lof,stkIOResponseLibrary_));
	observer->pushBack(solution_writer);
      }

      return observer;
    }

  private:
    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

  };

}

#endif
