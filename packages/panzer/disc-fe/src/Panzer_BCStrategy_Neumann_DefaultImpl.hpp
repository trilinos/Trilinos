// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_NEUMANN_DEFAULT_IMPL_DECL_HPP
#define PANZER_BCSTRATEGY_NEUMANN_DEFAULT_IMPL_DECL_HPP

#include <vector>
#include <string>
#include <tuple>

#include "Teuchos_RCP.hpp"

#include "Panzer_BCStrategy.hpp"
#include "Panzer_GlobalDataAcceptor_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_Traits.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
  
  template <typename EvalT>
  class BCStrategy_Neumann_DefaultImpl : public panzer::BCStrategy<EvalT>,
					 public panzer::GlobalDataAcceptorDefaultImpl,
					 public panzer::EvaluatorWithBaseImpl<panzer::Traits>
  {
    
  public:    
    
    BCStrategy_Neumann_DefaultImpl(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);
    
    virtual ~BCStrategy_Neumann_DefaultImpl();
    
    //! \name Derived from BCStrategy
    //@{ 

    virtual void setup(const panzer::PhysicsBlock& side_pb, const Teuchos::ParameterList& user_data) = 0;
      
    virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					    const panzer::PhysicsBlock& side_pb,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					    const Teuchos::ParameterList& models,
					    const Teuchos::ParameterList& user_data) const = 0;

    virtual void 
    buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					    const panzer::PhysicsBlock& side_pb,
					    const panzer::LinearObjFactory<panzer::Traits> & lof,
					    const Teuchos::ParameterList& user_data) const;

    virtual void 
    buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::PhysicsBlock& side_pb,
				      const LinearObjFactory<panzer::Traits> & lof,
				      const Teuchos::ParameterList& user_data) const;

    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					           const panzer::PhysicsBlock& side_pb,
						   const LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const;
    //@}

    //! \name Derived from PHX::EvaluatorWithDefaultImpl
    //@{ 
    virtual void postRegistrationSetup(typename Traits::SetupData d,
				       PHX::FieldManager<Traits>& vm) = 0;

    virtual void evaluateFields(typename Traits::EvalData d) = 0;
    //@}

    //! \name User Interface methods to provide information to the default implementation to be able to build the default evaluators for a Neumann BC
    //@{ 

    //! Requires that a gather evaluator for the DOF be constructed.
    virtual void requireDOFGather(const std::string required_dof_name);

    /** \brief Adds a residual contribution for a neumann condition to a particular equation
	
      \param residual_name [in] Name of the residual field that is to be scattered to the global residual.
      \param dof_name [in] Name of the DOF residual that the Neumann contribution should be added to.
      \param flux_name [in] Name of the flux field that will be integrated over to form the residual.
      \param integration_order [in] Order of the integration rule needed to define the data layouts for the flux integration.
      \param side_pb [in] The side physics block.  Used to build the PureBasis and IntegrationRule for a residual contribution.
    */ 
    virtual void addResidualContribution(const std::string residual_name,
					 const std::string dof_name,
					 const std::string flux_name,
					 const int integration_order,
					 const panzer::PhysicsBlock& side_pb);

    //! Returns information for the residual contribution integrations associated with this Neumann BC.
    const std::vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > > getResidualContributionData() const;

    //@}

    
    //! \name Query methods for underlying data
    //@{ 

    //! Returns the boundary condition data for this object
    const panzer::BC bc() const;

    //@}

    
  private:

    //! \name Utility functions used by default implementation
    //@{ 

    //! Finds the basis for the corresponding dof_name in the physics block.
    Teuchos::RCP<panzer::PureBasis> getBasis(const std::string dof_name,
					     const panzer::PhysicsBlock& side_pb) const;

    //! Allocates and returns the integration rule associated with an integration order and side physics block.
    Teuchos::RCP<panzer::IntegrationRule> buildIntegrationRule(const int integration_order,
							       const panzer::PhysicsBlock& side_pb) const;
    //@}

  private:

    /** \brief A vector of tuples containing information for each residual contribution for a corresponding Neumann bc

        Each entry in the vector is a different contribution to an
        individual residual equation.  A single Neumann BC can
        contribute to multiple equations (e.g. all momentum
        equations).  The tuple has 6 entries, with the indices
        described below.

        \param index 0 - Name of the residual field that is to be scattered to the global residual.
	\param index 1 - Name of the DOF residual that the Neumann residual contribution should be added to.
	\param index 2 - Name of the flux field that will be integrated over to form the residual.
	\param index 3 - Order of the integration rule needed to define the data layouts for the flux integration.
	\param index 4 - RCP to the basis for the dof that the Neumann residual contribution is added to.
	\param index 5 - RCP to the Integration Rule for the side integration for the Neumann residual contribution.
    */
    std::vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > > m_residual_contributions;

    //! All DOF field names needed by this BC: this vector is used to build gather evaluators for each DOF.
    std::vector<std::string> m_required_dof_names;

  };

}

#include "Panzer_BCStrategy_Neumann_DefaultImpl_impl.hpp"

#endif
