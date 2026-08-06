// @HEADER BEGIN
/**********************************************************************************

EMPIRE 

Copyright (c) 2015, Sandia National Laboratories
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For questions, comments or contributions contact 
Matt Bettencourt, mbetten@sandia.gov

*******************************************************************************/
// @HEADER END


#ifndef _MiniEM_SOLVERS_AuxiliaryEquationSet_MACROS_hpp_
#define _MiniEM_SOLVERS_AuxiliaryEquationSet_MACROS_hpp_
#include <iostream>
#include <string>
#include "Teuchos_Assert.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"

/** \file
  * Macros for declaring and instantiating panzer::EquationSetFactory
  * "TemplateBuilder" functors for mini-em's auxiliary equation sets
  * (see mini_em::EquationSetFactory), analogous to Panzer's own
  * PANZER_DECLARE_EQSET_TEMPLATE_BUILDER/PANZER_BUILD_EQSET_OBJECTS
  * but threading an additional panzer::GlobalEvaluationDataContainer
  * through the builder, since auxiliary equation sets need it to
  * scatter their assembled operators.
  *
  * AUX_DECLARE_EQSET_TEMPLATE_BUILDER and AUX_BUILD_EQSET_OBJECTS are
  * the macros used by mini_em::EquationSetFactory.
  */

/// \brief Declares a `<fType>_TemplateBuilder` functor that constructs one `fClass<EvalT>` auxiliary equation set per evaluation type, threading a panzer::GlobalEvaluationDataContainer through the constructor.
#undef AUX_DECLARE_EQSET_TEMPLATE_BUILDER
#define AUX_DECLARE_EQSET_TEMPLATE_BUILDER(fClass, fType)                                   \
                                                                                            \
  struct fType ## _TemplateBuilder                                                          \
  {                                                                                         \
    const Teuchos::RCP<Teuchos::ParameterList>                m_params;                     \
    const int                                                 m_default_integration_order;  \
    const panzer::CellData&                                   m_cell_data;                  \
    const Teuchos::RCP<panzer::GlobalData>                    m_global_data;                \
    const bool                                                m_build_transient_support;    \
    const Teuchos::RCP<panzer::GlobalEvaluationDataContainer> m_gedc;                       \
    fType ## _TemplateBuilder(                                                              \
      const Teuchos::RCP<Teuchos::ParameterList>&                params,                    \
      const int                                                  default_integration_order, \
      const panzer::CellData&                                    cd,                        \
      const Teuchos::RCP<panzer::GlobalData>&                    global_data,               \
      const bool                                                 build_transient_support,   \
      const Teuchos::RCP<panzer::GlobalEvaluationDataContainer>& gedc)                      \
      :                                                                                     \
      m_params(params),                                                                     \
      m_default_integration_order(default_integration_order),                               \
      m_cell_data(cd),                                                                      \
      m_global_data(global_data),                                                           \
      m_build_transient_support(build_transient_support),                                   \
      m_gedc(gedc)                                                                          \
    {                                                                                       \
    }                                                                                       \
                                                                                            \
    template<typename EvalT>                                                                \
    Teuchos::RCP<panzer::EquationSetBase> build() const                                     \
    {                                                                                       \
      fClass<EvalT>* ptr = new fClass<EvalT>(m_gedc, m_params, m_default_integration_order, \
        m_cell_data, m_global_data, m_build_transient_support);                             \
      return Teuchos::rcp(ptr);                                                             \
    }                                                                                       \
                                                                                            \
  };


/// \brief If params's "Type" equals key, builds the `fType` auxiliary equation set (via its `_TemplateBuilder`, declared by AUX_DECLARE_EQSET_TEMPLATE_BUILDER) into eq_set and sets found=true. Intended to be chained, one invocation per recognized key, inside panzer::EquationSetFactory::buildEquationSet().
#undef AUX_BUILD_EQSET_OBJECTS
#define AUX_BUILD_EQSET_OBJECTS(key, fType)                                    \
  if (params->get<std::string>("Type") == key)                                 \
  {                                                                            \
    fType ## _TemplateBuilder builder(params, default_integration_order,       \
      cell_data, global_data, build_transient_support, m_gedc); \
    eq_set->buildObjects(builder);                                             \
    found = true;                                                              \
  }

#endif
