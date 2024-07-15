// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step01_EquationSet_Projection_hpp__
#define __Step01_EquationSet_Projection_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace user_app {

template <typename EvalT>
class EquationSet_Projection : public panzer::EquationSet_DefaultImpl<EvalT> {

public:    

  EquationSet_Projection(const Teuchos::RCP<Teuchos::ParameterList>& params,
                         const int& default_integration_order,
                         const panzer::CellData& cell_data,
                         const Teuchos::RCP<panzer::GlobalData>& gd,
                         const bool build_transient_support);
  
  void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                             const panzer::FieldLibrary& field_library,
                                             const Teuchos::ParameterList& user_data) const;

private:

  std::string dof_name_;
};

}

#include "Step01_EquationSet_Projection_impl.hpp"

#endif
