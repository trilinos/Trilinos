// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_EquationSet_Darcy_hpp_
#define _MiniEM_EquationSet_Darcy_hpp_

#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace mini_em {

  template <typename EvalT>
    class EquationSet_Darcy : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:

    EquationSet_Darcy(const Teuchos::RCP<Teuchos::ParameterList>& params,
           const int& default_integration_order,
           const panzer::CellData& cell_data,
           const Teuchos::RCP<panzer::GlobalData>& gd,
           const bool build_transient_support);

      void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
             const panzer::FieldLibrary& field_library,
             const Teuchos::ParameterList& user_data) const;
      std::string EFieldName() const {return m_u_field_dof_name;}
      std::string BFieldName() const {return m_sigma_field_dof_name;}
  private:
      int dimension;
      std::string m_u_field_dof_name;
      std::string m_sigma_field_dof_name;
      std::string inverse_diffusivity_, forcing_;
 };

}

#include "MiniEM_EquationSet_Darcy_impl.hpp"



#endif
