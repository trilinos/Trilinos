// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"

#undef PANZER_DECLARE_EQSET_TEMPLATE_BUILDER
#define PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(fClass, fType)	               \
  									                                                         \
  struct fType ## _TemplateBuilder                                           \
  {					                                                                 \
    const Teuchos::RCP<Teuchos::ParameterList> m_params;                     \
    const int                                  m_default_integration_order;  \
    const panzer::CellData&                    m_cell_data;                  \
    const Teuchos::RCP<panzer::GlobalData>     m_global_data;                \
    const bool                                 m_build_transient_support;    \
    fType ## _TemplateBuilder(                                               \
      const Teuchos::RCP<Teuchos::ParameterList>& params,                    \
      const int                                   default_integration_order, \
      const panzer::CellData&                     cd,                        \
      const Teuchos::RCP<panzer::GlobalData>&     global_data,               \
      const bool                                  build_transient_support)   \
      :                                                                      \
      m_params(params),                                                      \
      m_default_integration_order(default_integration_order),                \
      m_cell_data(cd),                                                       \
      m_global_data(global_data),                                            \
      m_build_transient_support(build_transient_support)                     \
    {                                                                        \
    }                                                                        \
                                                                             \
    template<typename EvalT>                                                 \
    Teuchos::RCP<panzer::EquationSetBase> build() const                      \
    {          	                                                             \
      fClass<EvalT>* ptr = new fClass<EvalT>(m_params,                       \
        m_default_integration_order, m_cell_data, m_global_data,             \
        m_build_transient_support);                                          \
      return Teuchos::rcp(ptr);						                                   \
    }                                                                        \
                                                                             \
  };

#undef PANZER_BUILD_EQSET_OBJECTS
#define PANZER_BUILD_EQSET_OBJECTS(key, fType)                           \
  if (params->get<std::string>("Type") == key)                           \
  {				                                                               \
    fType ## _TemplateBuilder builder(params, default_integration_order, \
      cell_data, global_data, build_transient_support);                  \
    eq_set->buildObjects(builder);				                               \
    found = true;                                                        \
  }
