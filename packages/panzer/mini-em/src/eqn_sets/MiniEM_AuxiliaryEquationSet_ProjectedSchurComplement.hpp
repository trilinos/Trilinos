// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_AuxiliaryEquationSet_ProjectedSchurComplement_hpp_
#define _MiniEM_AuxiliaryEquationSet_ProjectedSchurComplement_hpp_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace mini_em {

  /** \brief An auxiliary equation set (not part of the physical
    * system being solved) that assembles a gradient-projected
    * approximate Schur complement operator for the Maxwell E-field
    * block (using gradient-basis integrators, rather than directly in
    * the edge-element space like mini_em::AuxiliaryEquationSet_SchurComplement)
    * and scatters it to the global evaluation data container under
    * the key "AUX_<dof_name>_ProjectedSchurComplement", for use by
    * mini_em::OperatorRequestCallback and the Maxwell block
    * preconditioner factories.
    */
  template <typename EvalT>
  class AuxiliaryEquationSet_ProjectedSchurComplement : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:

    /// \brief Construct from the global evaluation data container to scatter the projected Schur complement operator into, and a ParameterList specifying the basis/integration order, permittivity, conductivity, and inverse permeability closure model names.
    AuxiliaryEquationSet_ProjectedSchurComplement(const Teuchos::RCP<panzer::GlobalEvaluationDataContainer> & gedc,
                             const Teuchos::RCP<Teuchos::ParameterList>& params,
			     const int& default_integration_order,
			     const panzer::CellData& cell_data,
		             const Teuchos::RCP<panzer::GlobalData>& global_data,
		             const bool build_transient_support);

    /** \brief Registers the evaluators that assemble the projected approximate Schur complement operator per the class description.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] field_library field data layouts available for this physics block (unused).
      * \param[in] user_data user-supplied parameters (unused).
      */
    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					       const panzer::FieldLibrary& /* field_library */,
                                               const Teuchos::ParameterList& /* user_data */) const;

    /** \brief Registers the evaluators that scatter the assembled operator into the global evaluation data container per the class description.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] field_library field data layouts available for this physics block.
      * \param[in] lof linear object factory used to build the scatter evaluator.
      * \param[in] user_data user-supplied parameters passed through to the scatter evaluator.
      */
    void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					   const panzer::FieldLibrary& field_library,
                                           const panzer::LinearObjFactory<panzer::Traits> & lof,
                                           const Teuchos::ParameterList& user_data) const;

  protected:
    std::string dof_name;
    Teuchos::RCP<panzer::GlobalEvaluationDataContainer> m_gedc;
    Teuchos::RCP<std::vector<std::string> > m_dof_names;
    std::string permittivity_, conductivity_, inversePermeability_;
  };

/// \brief Jacobian specialization: this is the only evaluation type for which the projected Schur complement operator is actually assembled and scattered.
template < >
void AuxiliaryEquationSet_ProjectedSchurComplement<panzer::Traits::Jacobian>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				  const panzer::FieldLibrary& field_library,
                                  const panzer::LinearObjFactory<panzer::Traits> & lof,
                                  const Teuchos::ParameterList& user_data) const;

}

#endif
