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
#include "Teuchos_RCP.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"

#undef PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER
#define PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(fClass, fType)  \
                                                                   \
  struct fType ## _TemplateBuilder                                 \
  {                                                                \
    const panzer::BC&                      m_bc;                   \
    const Teuchos::RCP<panzer::GlobalData> m_global_data;          \
    fType ## _TemplateBuilder(                                     \
      const panzer::BC&                       bc,                  \
      const Teuchos::RCP<panzer::GlobalData>& global_data)         \
      :                                                            \
      m_bc(bc),                                                    \
      m_global_data(global_data)                                   \
    {                                                              \
    }                                                              \
                                                                   \
    template<typename EvalT>                                       \
    Teuchos::RCP<panzer::BCStrategyBase> build() const             \
    {                                                              \
      fClass<EvalT>* ptr = new fClass<EvalT>(m_bc, m_global_data); \
      return Teuchos::rcp(ptr);                                    \
    }                                                              \
                                                                   \
  };

#define PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER_EXTRA(fClass, fType, \
  extraSteps)                                                           \
                                                                        \
  struct fType ## _TemplateBuilder                                      \
  {                                                                     \
    const panzer::BC&                      m_bc;                        \
    const Teuchos::RCP<panzer::GlobalData> m_global_data;               \
    fType ## _TemplateBuilder(                                          \
      const panzer::BC&                       bc,                       \
      const Teuchos::RCP<panzer::GlobalData>& global_data)              \
      :                                                                 \
      m_bc(bc),                                                         \
      m_global_data(global_data)                                        \
    {                                                                   \
    }                                                                   \
                                                                        \
    template<typename EvalT>                                            \
    Teuchos::RCP<panzer::BCStrategyBase> build() const                  \
    {                                                                   \
      fClass<EvalT>* ptr = new fClass<EvalT>(m_bc, m_global_data);      \
      {                                                                 \
        extraSteps                                                      \
      }                                                                 \
      return Teuchos::rcp(ptr);                                         \
    }                                                                   \
                                                                        \
  };

#undef PANZER_BUILD_BCSTRATEGY_OBJECTS
#define PANZER_BUILD_BCSTRATEGY_OBJECTS(key, fType)     \
  if (bc.strategy() == key)                             \
  {                                                     \
    fType ## _TemplateBuilder builder(bc, global_data); \
    bcs_tm->buildObjects(builder);                      \
    found = true;                                       \
  }
