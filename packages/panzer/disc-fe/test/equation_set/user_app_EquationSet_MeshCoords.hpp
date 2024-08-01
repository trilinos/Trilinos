// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_EQUATIONSET_MESHCOORDS_HPP
#define USER_APP_EQUATIONSET_MESHCOORDS_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace user_app {

  template <typename EvalT>
    class EquationSet_MeshCoords : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:    

    EquationSet_MeshCoords(const Teuchos::RCP<Teuchos::ParameterList>& params,
		       const int& default_integration_order,
		       const panzer::CellData& cell_data,
		       const Teuchos::RCP<panzer::GlobalData>& gd,
		       const bool build_transient_support);
    
      void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						 const panzer::FieldLibrary& field_library,
						 const Teuchos::ParameterList& user_data) const;
  private:
    int dimension_; 

  };

}

#include "user_app_EquationSet_MeshCoords_impl.hpp"

#endif
