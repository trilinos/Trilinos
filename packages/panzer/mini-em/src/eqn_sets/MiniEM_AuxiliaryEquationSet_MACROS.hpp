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


#undef AUX_BUILD_EQSET_OBJECTS
#define AUX_BUILD_EQSET_OBJECTS(key, fType)                                    \
  if (params->get<std::string>("Type") == key)                                 \
  {                                                                            \
    fType ## _TemplateBuilder builder(params, default_integration_order,       \
      cell_data, global_data, build_transient_support, m_gedc); \
    eq_set->buildObjects(builder);                                             \
    found = true;                                                              \
  }

#undef DREKAR_DECLARE_EQSET_TEMPLATE_BUILDER_GENERAL
#define DREKAR_DECLARE_EQSET_TEMPLATE_BUILDER_GENERAL(fClass, fType)	         \
  									                                                           \
  struct fType ## _TemplateBuilder {					                                 \
    const Teuchos::RCP<Teuchos::ParameterList> m_params;                       \
    const int                                  m_default_integration_order;    \
    const panzer::CellData&                    m_cell_data;                    \
    const Teuchos::RCP<panzer::GlobalData>     m_global_data;                  \
    const bool                                 m_build_transient_support;      \
    const bool                                 m_do_semi_implicit;             \
    const std::string                          m_semi_implicit_predictor_name; \
    const bool                                 m_do_lagging;                   \
    const std::string                          m_lagging_name;	               \
    const bool                                 m_do_high_order;                \
    fType ## _TemplateBuilder(                                                 \
      const Teuchos::RCP<Teuchos::ParameterList>& params,                      \
      const int default_integration_order,                                     \
      const panzer::CellData& cd,                                              \
      const Teuchos::RCP<panzer::GlobalData>& global_data,                     \
      const bool build_transient_support,                                      \
      const bool do_semi_implicit,                                             \
      const std::string semi_implicit_predictor_name,                          \
      const bool do_lagging,                                                   \
      const std::string lagging_name,                                          \
      const bool do_high_order)                                                \
      :                                                                        \
      m_params(params),                                                        \
      m_default_integration_order(default_integration_order),                  \
      m_cell_data(cd),                                                         \
      m_global_data(global_data),                                              \
      m_build_transient_support(build_transient_support),                      \
      m_do_semi_implicit(do_semi_implicit),                                    \
      m_semi_implicit_predictor_name(semi_implicit_predictor_name),            \
      m_do_lagging(do_lagging),                                                \
      m_lagging_name(lagging_name),                                            \
      m_do_high_order(do_high_order)                                           \
    {                                                                          \
    }                                                                          \
                                                                               \
    template<typename EvalT>                                                   \
    Teuchos::RCP<panzer::EquationSetBase> build() const                        \
    {           	                                                             \
      using Teuchos::RCP;                                                      \
      using Teuchos::rcp;                                                      \
      using Teuchos::rcp_dynamic_cast;                                         \
      RCP<fClass<EvalT> > ptr = rcp(new fClass<EvalT>(m_params,                \
        m_default_integration_order, m_cell_data, m_global_data,               \
        m_build_transient_support));                                           \
      Teuchos::RCP<drekar::MixIn_SemiImplicitSupport<EvalT> > semiImplicit =   \
        rcp_dynamic_cast<drekar::MixIn_SemiImplicitSupport<EvalT> >(ptr);      \
      if (semiImplicit != Teuchos::null)                                       \
        semiImplicit->setDoSemiImplicit(m_do_semi_implicit,                    \
          m_semi_implicit_predictor_name);                                     \
      Teuchos::RCP<drekar::MixIn_LaggingSupport<EvalT> > lagging =             \
        rcp_dynamic_cast<drekar::MixIn_LaggingSupport<EvalT> >(ptr);           \
      if (lagging != Teuchos::null)                                            \
        lagging->setDoLagging(m_do_lagging, m_lagging_name);                   \
      if (m_do_high_order)                                                     \
      {                                                                        \
        Teuchos::RCP<drekar::MixIn_HighOrderSupport<EvalT> > highOrder =       \
          rcp_dynamic_cast<drekar::MixIn_HighOrderSupport<EvalT> >(ptr);       \
        if (highOrder != Teuchos::null)                                        \
          highOrder->setDoHighOrder(m_do_high_order);                          \
        else                                                                   \
          std::cout << "Warning: High Order enabled, "                         \
                    << "Equation set \"EvalT\" does not support it!";          \
      }                                                                        \
      ptr->initialize(m_params, m_default_integration_order, m_cell_data,      \
        m_global_data, m_build_transient_support);                             \
      return ptr;                                                              \
    }									                                                         \
    									                                                         \
  };

#undef DREKAR_BUILD_EQSET_OBJECTS_GENERAL
#define DREKAR_BUILD_EQSET_OBJECTS_GENERAL(key, fType)                       \
  if (params->get<std::string>("Type") == key)                               \
  {                                                                          \
      fType ## _TemplateBuilder builder(params, default_integration_order,   \
        cell_data, global_data, build_transient_support, m_do_semi_implicit, \
        m_semi_implicit_predictor_name, m_do_lagging, m_lagging_name,        \
        m_do_high_order);                                                    \
      eq_set->buildObjects(builder);				                                 \
      found = true;                                                          \
  }


#endif
