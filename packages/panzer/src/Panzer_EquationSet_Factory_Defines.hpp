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

#include <iostream>
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"

#undef PANZER_DECLARE_EQSET_TEMPLATE_BUILDER
#define PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(key, fClass, fType)	\
  									\
  struct fType ## _TemplateBuilder {					\
    const Teuchos::RCP<Teuchos::ParameterList> m_params;                \
    const int m_default_integration_order;             			\
    const panzer::CellData& m_cell_data;                                \
    const Teuchos::RCP<panzer::GlobalData> m_global_data;               \
    const bool m_build_transient_support;                               \
    fType ## _TemplateBuilder(const Teuchos::RCP<Teuchos::ParameterList>& params, const int default_integration_order, const panzer::CellData& cd, const Teuchos::RCP<panzer::GlobalData>& global_data, const bool build_transient_support) : m_params(params), m_default_integration_order(default_integration_order), m_cell_data(cd), m_global_data(global_data), m_build_transient_support(build_transient_support) {} \
									\
    template<typename EvalT>						\
    Teuchos::RCP<panzer::EquationSetBase> build() const {           	\
      fClass <EvalT>* ptr = new fClass <EvalT>(m_params, m_default_integration_order, m_cell_data, m_global_data, m_build_transient_support); \
      return Teuchos::rcp(ptr);						\
    }									\
    									\
  };

#undef PANZER_BUILD_EQSET_OBJECTS
#define PANZER_BUILD_EQSET_OBJECTS(key, fClass, fType)                  \
  if (params->get<std::string>("Type") == key) {				\
      fType ## _TemplateBuilder builder(params, default_integration_order, cell_data, global_data, build_transient_support); \
      eq_set->buildObjects(builder);				        \
      found = true;                                                     \
    }
