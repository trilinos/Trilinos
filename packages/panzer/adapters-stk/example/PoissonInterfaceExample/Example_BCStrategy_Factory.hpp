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

#ifndef __Example_BCStrategyFactory_hpp__
#define __Example_BCStrategyFactory_hpp__

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"

// Add my bcstrategies here
#include "Example_BCStrategy_Dirichlet_Constant.hpp"
#include "Example_BCStrategy_Neumann_Constant.hpp"
#include "Example_BCStrategy_Interface_WeakDirichletMatch.hpp"
#include "Example_BCStrategy_Interface_NeumannMatch.hpp"
#include "Example_BCStrategy_Interface_Robin.hpp"

namespace Example {
  
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(BCStrategy_Dirichlet_Constant,
  BCStrategy_Dirichlet_Constant)
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(BCStrategy_Neumann_Constant,
  BCStrategy_Neumann_Constant)
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
  BCStrategy_Interface_WeakDirichletMatch,
  BCStrategy_Interface_WeakDirichletMatch)
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(BCStrategy_Interface_NeumannMatch,
  BCStrategy_Interface_NeumannMatch)
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(BCStrategy_Interface_Robin,
  BCStrategy_Interface_Robin)

struct BCStrategyFactory : public panzer::BCStrategyFactory {

   Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
   buildBCStrategy(const panzer::BC& bc,const Teuchos::RCP<panzer::GlobalData>& global_data) const
   {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
	Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant",
        BCStrategy_Dirichlet_Constant)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Constant",
        BCStrategy_Neumann_Constant)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Weak Dirichlet Match Interface",
        BCStrategy_Interface_WeakDirichletMatch)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Match Interface",
        BCStrategy_Interface_NeumannMatch)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Robin Interface",
        BCStrategy_Interface_Robin)

      TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
			 "Error - the BC Strategy called \"" << bc.strategy() <<
			 "\" is not a valid identifier in the BCStrategyFactory.  Either add a "
                         "valid implementation to your factory or fix your input file.  The "
                         "relevant boundary condition is:\n\n" << bc << std::endl);
      
      return bcs_tm;
   }

};
  
}

#endif
