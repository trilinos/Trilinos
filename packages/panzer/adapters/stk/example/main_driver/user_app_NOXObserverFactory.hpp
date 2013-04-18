// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_Traits.hpp"
#include "NOX_PrePostOperator_Vector.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// Individual Observers
#include "user_app_NOXObserver_WriteToExodus.hpp"
#include "user_app_NOXObserver_NeumannBCAnalyticSystemTest.hpp"

namespace user_app {
  
  class NOXObserverFactory :
    public panzer_stk::NOXObserverFactory,
    public Teuchos::ParameterListAcceptorDefaultBase  {
    
    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

    mutable Teuchos::RCP<Teuchos::ParameterList> valid_params_;

  public:

    NOXObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary)
       : stkIOResponseLibrary_(stkIOResponseLibrary) {}
    
    Teuchos::RCP<NOX::Abstract::PrePostOperator>
    buildNOXObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
		     const Teuchos::RCP<panzer::UniqueGlobalIndexerBase>& dof_manager,
		     const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::ParameterList;

      TEUCHOS_ASSERT(nonnull(this->getParameterList()));

      RCP<NOX::PrePostOperatorVector> observer = rcp(new NOX::PrePostOperatorVector);

      // Exodus writer to output solution
      if (this->getParameterList()->get<std::string>("Write Solution to Exodus File") == "ON") {
	Teuchos::RCP<NOX::Abstract::PrePostOperator> solution_writer = 
	  Teuchos::rcp(new user_app::NOXObserver_WriteToExodus(mesh,dof_manager,lof,stkIOResponseLibrary_));
	observer->pushBack(solution_writer);
      }

      
      // Neumann BC unit test
      if (this->getParameterList()->get<std::string>("Neumann BC Analytic System Test") == "ON") {
	Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = 
	  Teuchos::rcp(new user_app::NOXObserver_NeumannBCAnalyticSystemTest);
	observer->pushBack(ppo);
      }

      return observer;
    }

    /** \name Overridden from Teuchos::ParameterListAcceptor */
    //@{
    
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;

      paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
      setMyParamList(paramList);
    }
    
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const
    {
      if (valid_params_.is_null()) {
	
	valid_params_ = Teuchos::rcp(new Teuchos::ParameterList);

	Teuchos::setStringToIntegralParameter<int>(
          "Write Solution to Exodus File",
          "ON",
          "Enables or disables writing of solution to Exodus file at end of NOX solve",
          Teuchos::tuple<std::string>("ON","OFF"),
          valid_params_.get()
          );
	Teuchos::setStringToIntegralParameter<int>(
          "Neumann BC Analytic System Test",
          "OFF",
          "Checks solution values for Neumann BC Analytic System Test",
          Teuchos::tuple<std::string>("ON","OFF"),
          valid_params_.get()
          );

      }
      return valid_params_;
    }
 
    //@}

  };

}

#endif
