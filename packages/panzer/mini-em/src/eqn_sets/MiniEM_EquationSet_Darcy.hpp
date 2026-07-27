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

  /** \brief The mixed-form Darcy flow equation set:
    * du/dt - div(sigma) = f, -grad(u) + (1/kappa)*sigma = 0, where u is
    * the scalar potential (e.g. pressure) and sigma is the vector flux.
    */
  template <typename EvalT>
    class EquationSet_Darcy : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:

    /// \brief Construct from a ParameterList specifying the basis/integration order, diffusivity, inverse diffusivity, and forcing closure model names, as documented in the class description.
    EquationSet_Darcy(const Teuchos::RCP<Teuchos::ParameterList>& params,
           const int& default_integration_order,
           const panzer::CellData& cell_data,
           const Teuchos::RCP<panzer::GlobalData>& gd,
           const bool build_transient_support);

      /** \brief Registers the evaluators implementing the Darcy weak form per the class description.
        *
        * \param[in] fm the field manager to register evaluators with.
        * \param[in] field_library field data layouts available for this physics block.
        * \param[in] user_data user-supplied parameters passed through to the equation-set evaluators.
        */
      void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
             const panzer::FieldLibrary& field_library,
             const Teuchos::ParameterList& user_data) const;
      //! Returns the DOF name of the scalar potential field u. (Named EFieldName() for consistency with EquationSet_Maxwell.)
      std::string EFieldName() const {return m_u_field_dof_name;}
      //! Returns the DOF name of the vector flux field sigma. (Named BFieldName() for consistency with EquationSet_Maxwell.)
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
