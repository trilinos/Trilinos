// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_DECL_HPP
#define PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_DECL_HPP

#include <vector>
#include <string>
#include <map>

#include "Teuchos_RCP.hpp"

#include "Panzer_BCStrategy.hpp"
#include "Panzer_GlobalDataAcceptor_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"

#include "Phalanx_FieldManager.hpp"

namespace panzer {
  
  template <typename EvalT>
    class BCStrategy_Dirichlet_DefaultImpl : public panzer::BCStrategy<EvalT>,
					     public panzer::GlobalDataAcceptorDefaultImpl {

  public:    

    BCStrategy_Dirichlet_DefaultImpl(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data,
				     const bool check_apply_bc = false);
    
    virtual ~BCStrategy_Dirichlet_DefaultImpl();
    
    virtual void setup(const panzer::PhysicsBlock& side_pb, const Teuchos::ParameterList& user_data) = 0;
      
    virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					    const panzer::PhysicsBlock& pb,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					    const Teuchos::ParameterList& models,
					    const Teuchos::ParameterList& user_data) const = 0;

    void 
    virtual buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						    const panzer::PhysicsBlock& pb,
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

  protected:

    struct DOFDescriptor {
      DOFDescriptor()
        : dofName("")
        , coefficientResidual(false)
        , residualName(std::make_pair(false,""))
        , scatterName(std::make_pair(false,""))
        , targetName(std::make_pair(false,""))
        , timeDerivative(std::make_pair(false,"")) {}

      std::string dofName;
      bool coefficientResidual; // don't use evaluation to compute the residual
                                // use the coefficient instead. Note this is useful
                                // for vector basis functions
      std::pair<bool,std::string> residualName;
      std::pair<bool,std::string> scatterName;
      std::pair<bool,std::string> targetName;
      std::pair<bool,std::string> timeDerivative;

      void print(std::ostream & os) const {
        os << "BC DOF Desc = \"" << dofName << "\": "
           << "Coeffieint Residual = " << (coefficientResidual ? "true" : "false") << ", "
           << "Res = (" << residualName.first << ", \"" << residualName.second << ")\"), "
           << "Scatter = (" << scatterName.first << ", \"" << scatterName.second << "\"), "
           << "Scatter = (" << targetName.first << ", \"" << targetName.second << "\"), "
           << "Time = (" << timeDerivative.first << ", \"" << timeDerivative.second << "\")";
      }
    };

    /** This is to support backward compatibility and allow a migration path from the old
      * protected data approach of specifying the inputs and outputs of this class.
      */ 
    void buildDescriptorMapFromVectors() const;

    //! For convenience, declare the DOFDescriptor iterator
    typedef typename std::map<std::string,DOFDescriptor>::const_iterator DescriptorIterator;

    void addDOF(const std::string & dofName);

    /** Alert the panzer library that the DOF should be evaluated using a coefficient residual
      * as opposed to evaluating the basis and forcing the value to be equal at some point in
      * the element.
      *
      * \param[in] targetName (Required) Name of field that corresponds to the evaluated
      *                        dirichlet condition. This exists only in the PHX::FieldManager
      *                        and is required to be distinct from the <code>dofName</code>.
      * \param[in] dofName (Required) Name of field to lookup in the unique global
      *                        indexer. The 
      * \param[in] residualName (Optional) Name of field that is to be scattered associated with
      *                         this DOF.  If not supplied or an empty string used, the
      *                         default is to add the prefix "RESIDUAL_" to the dofName for
      *                         the residual field name.
      */
    void addCoefficientTarget(const std::string & targetName,
                              const std::string & dofName,
                              const std::string & residualName = "");

    /** Alert the panzer library of a DOF that is required by this boundary condition.
      * This automatically sets up the gather/scatter routines neccessary
      * to evaluate and assemble with this unknown.
      *
      * \param[in] targetName (Required) Name of field that corresponds to the evaluated
      *                        dirichlet condition. This exists only in the PHX::FieldManager
      *                        and is required to be distinct from the <code>dofName</code>.
      * \param[in] dofName (Required) Name of field to lookup in the unique global
      *                        indexer. The 
      * \param[in] residualName (Optional) Name of field that is to be scattered associated with
      *                         this DOF.  If not supplied or an empty string used, the
      *                         default is to add the prefix "RESIDUAL_" to the dofName for
      *                         the residual field name.
      */
    void addTarget(const std::string & targetName,
                   const std::string & dofName,
                   const std::string & residualName = "");

    /** Alert the panzer library that the time derivative of a DOF is required by this boundary condition.
      * This automatically sets up the gather/scatter routines neccessary
      * to evaluate and assemble with this unknown.
      *
      * \param[in] targetName (Required) Name of field that corresponds to the evaluated
      *                        dirichlet condition. This exists only in the PHX::FieldManager
      *                        and is required to be distinct from the <code>dofName</code>.
      * \param[in] dofName (Required) Name of field to lookup in the unique global
      *                        indexer. The 
      * \param[in] residualName (Optional) Name of field that is to be scattered associated with
      *                         this DOF.  If not supplied or an empty string used, the
      *                         default is to add the prefix "RESIDUAL_" to the dofName for
      *                         the residual field name.
      */
    void addDotTarget(const std::string & targetName,
                      const std::string & dofName,
                      const std::string & dotName="",
                      const std::string & residualName = "");

    std::map<std::string,DOFDescriptor> m_provided_dofs_desc;

    std::vector<std::string> required_dof_names;
    std::map<std::string,std::string> residual_to_dof_names_map;
    std::map<std::string,std::string> residual_to_target_field_map;
    bool check_apply_bc;
    mutable bool descriptor_map_built; // this is set in a const method
  };

}

#endif
