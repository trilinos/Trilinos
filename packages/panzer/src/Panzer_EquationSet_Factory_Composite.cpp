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

#include "Panzer_EquationSet_Factory_Composite.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_GlobalData.hpp"
#include <sstream>

namespace panzer {
  
  EquationSet_FactoryComposite::EquationSet_FactoryComposite(const std::vector<Teuchos::RCP<panzer::EquationSetFactory> >& factories) :
    m_factories(factories)
  { }
  
  Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
  EquationSet_FactoryComposite::buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& input_plist,
						 const int& default_integration_rule,
						 const panzer::CellData& cell_data,
						 const Teuchos::RCP<panzer::GlobalData>& global_data,
						 const bool build_transient_support) const
  {
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
    
    
    for (std::vector<Teuchos::RCP<panzer::EquationSetFactory> >::const_iterator factory = m_factories.begin();
	 factory != m_factories.end(); ++factory) {
      eq_set = (*factory)->buildEquationSet(input_plist,default_integration_rule,cell_data,global_data,build_transient_support);

      if (nonnull(eq_set))
	break;
    }
    
    std::ostringstream os;
    os << "The equation set factory failed to build and equation set for the following input parameter list.  Please correct the input list:\n";
    input_plist->print(os);
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(eq_set), std::logic_error,os.str());
    
    return eq_set;
  }
  
}

