// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_BCStrategy_Dirichlet_Constant_hpp_
#define _MiniEM_BCStrategy_Dirichlet_Constant_hpp_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Phalanx_FieldManager.hpp"

namespace mini_em {

  /** \brief A Dirichlet boundary condition strategy that sets the
    * degrees of freedom on the boundary to a constant value.
    */
  template <typename EvalT>
  class BCStrategy_Dirichlet_Constant : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {

  public:

    /// \brief Construct from the boundary condition descriptor and global data, as required by panzer::BCStrategy_Dirichlet_DefaultImpl.
    BCStrategy_Dirichlet_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    /// \brief Looks up the DOF's basis and residual field name for the given physics block. Called once before evaluator registration.
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& user_data);

    /** \brief Registers the evaluators that set the boundary DOF to the constant value given in the BC descriptor.
      *
      * \param[in] fm the field manager to register evaluators with.
      * \param[in] pb the physics block for the boundary condition's side set.
      * \param[in] factory closure model factory (unused).
      * \param[in] models closure model ParameterList (unused).
      * \param[in] user_data user-supplied parameters (unused).
      */
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& pb,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				    const Teuchos::ParameterList& models,
				    const Teuchos::ParameterList& user_data) const;

    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;

  };

}

#include "MiniEM_BCStrategy_Dirichlet_Constant_impl.hpp"

#endif
