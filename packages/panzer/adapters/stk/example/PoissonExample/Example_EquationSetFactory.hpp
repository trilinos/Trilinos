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

#ifndef __PoissonExample_EquationSetFactory_hpp__
#define __PoissonExample_EquationSetFactory_hpp__

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"

// Add my equation sets here
#include "Example_PoissonEquationSet.hpp"

namespace Example {

// A macro that defines a class to make construction of the equation sets easier
//   - The equation set is constructed over a list of automatic differention types
PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("Poisson", PoissonEquationSet, PoissonEquationSet)

// A user written factory that creates each equation set.  The key member here
// is buildEquationSet
class EquationSetFactory : public panzer::EquationSetFactory {
public:

   Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
   buildEquationSet(const panzer::InputEquationSet& ies,
                    const panzer::CellData& cell_data,
		    const Teuchos::RCP<panzer::GlobalData>& global_data,
                    const bool build_transient_support) const
   {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
         Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
         
      bool found = false; // this is used by PANZER_BUILD_EQSET_OBJECTS
         
      // macro checks if(ies.name=="Poisson") then an EquationSet_Energy object is constructed
      PANZER_BUILD_EQSET_OBJECTS("Poisson", PoissonEquationSet, PoissonEquationSet)
         
      // make sure your equation set has been found
      if(!found) {
         std::string msg = "Error - the \"Equation Set\" called \"" + ies.name +
                           "\" is not a valid equation set identifier. Please supply the correct factory.\n";
         TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
         
      return eq_set;
   }
    
};

}

#endif
