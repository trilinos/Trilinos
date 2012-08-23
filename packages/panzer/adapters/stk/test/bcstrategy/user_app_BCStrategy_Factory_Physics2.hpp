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

#ifndef USER_APP_BCSTRATEGY_FACTORY_PHYSICS2_HPP
#define USER_APP_BCSTRATEGY_FACTORY_PHYSICS2_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_GlobalData.hpp"

// Add my bcstrategies here
#include "user_app_BCStrategy_Dirichlet_Constant.hpp"

namespace user_app {
  
  /* commenting out since object is already declared in physics 1 factory 
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER("Constant 2", 
					     user_app::BCStrategy_Dirichlet_Constant,
					     BCStrategy_Dirichlet_Constant)
  */

  class BCFactoryPhysics2 : public panzer::BCStrategyFactory {

    bool m_throw_on_failure;

  public:

    BCFactoryPhysics2(bool throw_on_failure) :
      m_throw_on_failure(throw_on_failure)
    { }

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
	Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant 2", 
				      user_app::BCStrategy_Dirichlet_Constant,
				      BCStrategy_Dirichlet_Constant);

      if (m_throw_on_failure)
	TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
				   "Error - the BC Strategy Physics2 called \"" << bc.strategy() <<
				   "\" is not a valid identifier in the BCStrategyFactory.  Either add a valid implementation to your factory or fix your input file.  The relevant boundary condition is:\n\n" << bc << std::endl);
      
      if (!found)
	return Teuchos::null;
      
      return bcs_tm;
    }

  };
  
}

#endif
