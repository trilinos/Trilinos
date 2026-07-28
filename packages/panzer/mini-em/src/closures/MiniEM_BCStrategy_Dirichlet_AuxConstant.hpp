// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MiniEM_BCSTRATEGY_DIRICHLET_AUXCONSTANT_DECL_HPP
#define MiniEM_BCSTRATEGY_DIRICHLET_AUXCONSTANT_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BCStrategy.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace mini_em {

  /** \brief A Dirichlet boundary condition strategy that sets the
    * boundary DOFs of an auxiliary equation set's operator (rather
    * than a primary physics equation) to a constant value.
    */
  template <typename EvalT>
  class BCStrategy_Dirichlet_AuxConstant : public panzer::BCStrategy<EvalT> {

  public:

    /// \brief Construct from the boundary condition descriptor, as required by panzer::BCStrategy.
    BCStrategy_Dirichlet_AuxConstant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& /* global_data */);

    /// \brief Looks up the auxiliary operator's DOF basis for the given physics block. Called once before evaluator registration.
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& /* user_data */);

    /** \brief Registers the evaluators that set the boundary DOF to the constant value given in the BC descriptor.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] pb the physics block for the boundary condition's side set (unused).
      * \param[in] factory closure model factory (unused).
      * \param[in] models closure model ParameterList (unused).
      * \param[in] user_data user-supplied parameters (unused).
      */
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& /* pb */,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
				    const Teuchos::ParameterList& /* models */,
				    const Teuchos::ParameterList& /* user_data */) const;

    /** \brief Registers the evaluators that gather and scatter the auxiliary operator's boundary DOF values.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] side_pb the physics block for the boundary condition's side set.
      * \param[in] lof linear object factory used to build the gather/scatter evaluators.
      * \param[in] user_data user-supplied parameters passed through to the gather/scatter evaluators.
      */
    void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                 const panzer::PhysicsBlock& side_pb,
                                                 const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                 const Teuchos::ParameterList& user_data) const;

    /** \brief Registers the evaluators that scatter the auxiliary operator's boundary DOF values.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] side_pb the physics block for the boundary condition's side set.
      * \param[in] lof linear object factory used to build the scatter evaluators.
      * \param[in] user_data user-supplied parameters passed through to the scatter evaluators.
      */
    virtual void
    buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::PhysicsBlock& side_pb,
				      const panzer::LinearObjFactory<panzer::Traits> & lof,
				      const Teuchos::ParameterList& user_data) const;

    /** \brief Registers the evaluators that gather the auxiliary operator's boundary DOF values and apply basis orientations.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] side_pb the physics block for the boundary condition's side set.
      * \param[in] lof linear object factory used to build the gather evaluators.
      * \param[in] user_data user-supplied parameters passed through to the gather evaluators.
      */
    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					           const panzer::PhysicsBlock& side_pb,
						   const panzer::LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const;

    std::string operatorName_;
    double value_;
    std::string fieldName_;
    Teuchos::RCP<panzer::PureBasis> basis_;

  };

/// \brief Jacobian specialization of buildAndRegisterGatherScatterEvaluators(), as documented on the primary template.
template < >
void BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const;

/// \brief Jacobian specialization of buildAndRegisterScatterEvaluators(), as documented on the primary template.
template < >
void BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const;
/// \brief Jacobian specialization of buildAndRegisterGatherAndOrientationEvaluators(), as documented on the primary template.
template < >
void BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const;

}

#include "MiniEM_BCStrategy_Dirichlet_AuxConstant_impl.hpp"

#endif
