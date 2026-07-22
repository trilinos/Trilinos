// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ASSEMBLY_ENGINE_TEMPLATE_BUILDER_HPP
#define PANZER_ASSEMBLY_ENGINE_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"

namespace panzer {

  class AssemblyEngine_TemplateBuilder {
    
    Teuchos::RCP<panzer::FieldManagerBuilder> m_fmb;
    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > m_lof;
    
  public:
    
    AssemblyEngine_TemplateBuilder(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                                   const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof) :
      m_fmb(fmb), m_lof(lof)  
      {}
      
    template <typename EvalT>
      Teuchos::RCP<panzer::Base> build() const {
      return Teuchos::rcp( static_cast<panzer::Base*>
			   (new panzer::AssemblyEngine<EvalT>(m_fmb,m_lof)) );
    }
    
  };
  
}

#endif 
