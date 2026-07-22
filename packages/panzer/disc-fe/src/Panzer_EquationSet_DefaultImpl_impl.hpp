// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EQUATIONSET_DEFAULT_IMPL_IMPL_HPP
#define PANZER_EQUATIONSET_DEFAULT_IMPL_IMPL_HPP

#include "Panzer_DOF.hpp"
#include "Panzer_DOF_PointValues.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DOFCurl.hpp"
#include "Panzer_DOFDiv.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_GatherIntegrationCoordinates.hpp"
#include "Panzer_GatherOrientation.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_ParameterList.hpp"

// For convenience we automate some evalautor registration
#include "Panzer_Sum.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(const Teuchos::RCP<Teuchos::ParameterList>& params,
                        const int& default_integration_order,
                        const panzer::CellData& cell_data,
                        const Teuchos::RCP<panzer::GlobalData>& global_data,
                        const bool build_transient_support) :
  panzer::GlobalDataAcceptorDefaultImpl(global_data),
  m_input_params(params),
  m_default_integration_order(default_integration_order),
  m_cell_data(cell_data),
  m_build_transient_support(build_transient_support)
{ 
  TEUCHOS_ASSERT(nonnull(m_input_params));

  m_eval_plist = Teuchos::rcp(new Teuchos::ParameterList);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setupDOFs()
{
  // Get defaultparameters
  m_type = m_input_params->get<std::string>("Type");

  // for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
  //     itr!=m_provided_dofs_desc.end();++itr) {
  //   itr->second.print(std::cout); std::cout << std::endl;
  // }

  this->m_provided_dofs.clear();
  this->m_int_rules.clear();

  // load up the provided dofs and unique int rules from the descriptor map
  for(typename std::map<std::string,DOFDescriptor>::iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {

    // Create the bases
    TEUCHOS_ASSERT(nonnull(itr->second.basis));
    this->m_provided_dofs.push_back(std::make_pair(itr->first, itr->second.basis));

    //{
    //  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    //  out.setOutputToRootOnly(0);
    //  itr->second.print(out);
    //  out << "Int Order = " << itr->second.integrationOrder << " (" << itr->second.intRule->order() << ")" << std::endl;
    //}

    // Create the unique integration rule map and complete descriptor objects
    TEUCHOS_ASSERT(nonnull(itr->second.intRule));
    m_int_rules[itr->second.intRule->order()] = itr->second.intRule;
    
  }

  // Setup the basis to dof mapping
  for (DescriptorIterator dof_iter = m_provided_dofs_desc.begin(); dof_iter != m_provided_dofs_desc.end(); ++dof_iter) {

    std::string basis_name = dof_iter->second.basis->name();
    Teuchos::RCP<panzer::PureBasis> basis = dof_iter->second.basis;
    std::string dof_name = dof_iter->first;

    if (is_null(m_basis_to_dofs[basis_name].first)) {
      m_basis_to_dofs[basis_name].first = basis;
      m_basis_to_dofs[basis_name].second = Teuchos::rcp(new std::vector<std::string>);
    }

    m_basis_to_dofs[basis_name].second->push_back(dof_name);
  }

  // Generate a unique list of bases
  for (DescriptorIterator dof_iter = m_provided_dofs_desc.begin(); dof_iter != m_provided_dofs_desc.end(); ++dof_iter) {
    m_unique_bases[dof_iter->second.basis->name()] = dof_iter->second.basis;
  }

  // Setup the default parameter list for closure models
  this->m_eval_plist->set("Block ID", getElementBlockId());
  this->setupDeprecatedDOFsSupport();
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::setupDeprecatedDOFsSupport()
{
  TEUCHOS_ASSERT(m_provided_dofs.size() > 0);
  TEUCHOS_ASSERT(m_int_rules.size() > 0);

  // Deprecated support assumes all equations in set use the same
  // basis and integration rule
  Teuchos::RCP<panzer::PureBasis> pure_basis = m_provided_dofs.begin()->second;
  Teuchos::RCP<panzer::IntegrationRule> int_rule = m_int_rules.begin()->second;
  Teuchos::RCP<panzer::BasisIRLayout> basis = panzer::basisIRLayout(pure_basis,*int_rule);

  this->m_eval_plist->set("Basis", basis);
  this->m_eval_plist->set("IR", int_rule);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::FieldLibrary& /* fl */,
                                               const LinearObjFactory<panzer::Traits> & lof,
                                               const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // ********************
  // DOFs (unknowns)
  // ********************

  // Gather, includes construction of orientation gathers
  for (BasisIterator basis_it = m_basis_to_dofs.begin(); basis_it != m_basis_to_dofs.end(); ++basis_it) {

    // Set tangent field names (first dimension is DOF, second is parameter)
    Teuchos::RCP< std::vector< std::vector<std::string> > > tangent_field_names;
    if (m_tangent_param_names.size() > 0) {
      tangent_field_names = rcp(new std::vector< std::vector<std::string> >(basis_it->second.second->size()));
      for (std::size_t i=0; i<basis_it->second.second->size(); ++i) {
        for (std::size_t j=0; j<m_tangent_param_names.size(); ++j) {
          const std::string tname =
            (*(basis_it->second.second))[i] + " SENSITIVITY " + m_tangent_param_names[j];
          (*tangent_field_names)[i].push_back(tname);
        }
      }
    }

    {
      ParameterList p("Gather");
      p.set("Basis", basis_it->second.first);
      p.set("DOF Names", basis_it->second.second);
      p.set("Indexer Names", basis_it->second.second);
      p.set("Sensitivities Name", "");
      p.set("First Sensitivities Available", true);
      p.set("Second Sensitivities Available", true);

      // Set tangent field names
      if (tangent_field_names != Teuchos::null)
        p.set("Tangent Names", tangent_field_names);

      RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // Create a second gather evaluator for each tangent field,
    // we never compute derivatives with respect to this field
    if (tangent_field_names != Teuchos::null) {
      for (std::size_t i=0; i<m_tangent_param_names.size(); ++i) {

        Teuchos::RCP< std::vector<std::string> > names = rcp(new std::vector<std::string>);
        for (std::size_t j=0; j<basis_it->second.second->size(); ++j)
          names->push_back((*tangent_field_names)[j][i]);

        ParameterList p(std::string("Gather Tangent ") + this->m_tangent_param_names[i]);
        p.set("Basis", basis_it->second.first);
        p.set("DOF Names", names);
        p.set("Indexer Names", basis_it->second.second);
        p.set("Sensitivities Name", "");
        p.set("First Sensitivities Available", false);
        p.set("Second Sensitivities Available", false);
        p.set("Global Data Key", "X TANGENT GATHER CONTAINER: " + this->m_tangent_param_names[i]);

        RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGatherTangent<EvalT>(p);

        this->template registerEvaluator<EvalT>(fm, op);
      }
    }
  }

  // **************************
  // Coordinates for integration points and basis functions
  // **************************
  {
    // add basis coordinates
    for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator basis =  m_unique_bases.begin();
         basis != m_unique_bases.end(); ++ basis) {
      RCP< PHX::Evaluator<panzer::Traits> > basis_op
        = rcp(new panzer::GatherBasisCoordinates<EvalT,panzer::Traits>(*basis->second));
      this->template registerEvaluator<EvalT>(fm, basis_op);
    }

    // add integration coordinates
    for (std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator ir = m_int_rules.begin();
         ir != m_int_rules.end(); ++ir)   {
      RCP< PHX::Evaluator<panzer::Traits> > quad_op
        = rcp(new panzer::GatherIntegrationCoordinates<EvalT,panzer::Traits>(*ir->second));
      this->template registerEvaluator<EvalT>(fm, quad_op);
    }

    // NOTE: You can look up the name of either coordinate field name by doing
    //       GatherBasisCoordinates<EvalT,Traits>::fieldName();
    //       GatherIntegrationCoordinates<EvalT,Traits>::fieldName();
  }

  // **************************
  // Time derivative terms
  // **************************

  // Gather of time derivative terms: One evaluator for each unique basis
  for (BasisIterator basis_it = m_basis_to_dofs.begin(); basis_it != m_basis_to_dofs.end(); ++basis_it) {

    RCP< std::vector<std::string> > t_dof_names = rcp(new std::vector<std::string>);   // time derivative indexer names
    RCP< std::vector<std::string> > t_field_names = rcp(new std::vector<std::string>); // time derivative field names
    RCP< std::vector< std::vector<std::string> > > tangent_field_names = rcp(new std::vector< std::vector<std::string> >); // tangent field names

    // determine which fields associated with this basis need time derivatives
    for (typename std::vector<std::string>::const_iterator dof_name = basis_it->second.second->begin();
         dof_name != basis_it->second.second->end(); ++dof_name) {

      DescriptorIterator desc = m_provided_dofs_desc.find(*dof_name);
      TEUCHOS_ASSERT(desc != m_provided_dofs_desc.end());

      // does this field need a time derivative?
      if(desc->second.timeDerivative.first) {
        // time derivative needed
        t_dof_names->push_back(*dof_name);
        t_field_names->push_back(desc->second.timeDerivative.second);

        // Set tangent field names (first dimension is DOF, second is parameter)
        if (m_tangent_param_names.size() > 0) {
          std::vector<std::string> tfn;
          for (std::size_t j=0; j<m_tangent_param_names.size(); ++j) {
            const std::string tname =
              desc->second.timeDerivative.second + " SENSITIVITY " + m_tangent_param_names[j];
            tfn.push_back(tname);
          }
          tangent_field_names->push_back(tfn);
        }
      }
    }

    {
      ParameterList p("Gather");
      p.set("Basis", basis_it->second.first);
      p.set("DOF Names", t_field_names);
      p.set("Indexer Names", t_dof_names);
      p.set("Use Time Derivative Solution Vector", true);

      // Set tangent field names
      if (m_tangent_param_names.size() > 0)
        p.set("Tangent Names", tangent_field_names);

      RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // Create a second gather evaluator for each tangent field,
    // we never compute derivatives with respect to this field
    if (m_tangent_param_names.size() > 0) {
      for (std::size_t i=0; i<m_tangent_param_names.size(); ++i) {

        Teuchos::RCP< std::vector<std::string> > names =
          rcp(new std::vector<std::string>);
        for (std::size_t j=0; j<tangent_field_names->size(); ++j)
          names->push_back((*tangent_field_names)[j][i]);

        ParameterList p(std::string("Gather Tangent ") + this->m_tangent_param_names[i]);
        p.set("Basis", basis_it->second.first);
        p.set("DOF Names", names);
        p.set("Indexer Names", t_dof_names);
        p.set("Use Time Derivative Solution Vector", true);
        p.set("Global Data Key", "DXDT TANGENT GATHER CONTAINER: " + this->m_tangent_param_names[i]);

        RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGatherTangent<EvalT>(p);

        this->template registerEvaluator<EvalT>(fm, op);
      }
    }
  }

  // **************************
  // Orientation terms
  // **************************

  for (BasisIterator basis_it = m_basis_to_dofs.begin(); basis_it != m_basis_to_dofs.end(); ++basis_it) {
    if(basis_it->second.first->requiresOrientations())  {
      ParameterList p("Gather Orientation");
      p.set("Basis", basis_it->second.first);
      p.set("DOF Names", basis_it->second.second);
      p.set("Indexer Names", basis_it->second.second);

      RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGatherOrientation<EvalT>(p);

      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                             const panzer::FieldLayoutLibrary& fl,
                                             const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                             const Teuchos::Ptr<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                             const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer;
  if(lof!=Teuchos::null) 
    globalIndexer = lof->getRangeGlobalIndexer();
  
  // DOFs: Scalar value @ basis --> Scalar value @ IP 
  for (DescriptorIterator dof_iter = m_provided_dofs_desc.begin(); dof_iter != m_provided_dofs_desc.end(); ++dof_iter) {

    ParameterList p;
    p.set("Name", dof_iter->first);
    p.set("Basis", fl.lookupLayout(dof_iter->first));
    p.set("IR", ir);

    if(globalIndexer!=Teuchos::null) {
      // build the offsets for this field
      int fieldNum = globalIndexer->getFieldNum(dof_iter->first);
      RCP<const std::vector<int> > offsets = 
          rcp(new std::vector<int>(globalIndexer->getGIDFieldOffsets(m_block_id,fieldNum)));
      p.set("Jacobian Offsets Vector", offsets);
    }
    // else default to the slow DOF call
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Gradients of DOFs: Scalar value @ basis --> Vector value @ IP

  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {
    
    if(itr->second.basis->supportsGrad()) {
      
      // is gradient required for this variable
      if(!itr->second.grad.first) 
        continue; // its not required, quit the loop

      const std::string dof_name =      itr->first;
      const std::string dof_grad_name = itr->second.grad.second;

      ParameterList p;
      p.set("Name", dof_name);
      p.set("Gradient Name", dof_grad_name);
      p.set("Basis", fl.lookupLayout(dof_name)); 
      p.set("IR", ir);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
        rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

  // Curl of DOFs: Vector value @ basis --> Vector value @ IP (3D) or Scalar value @ IP (2D)

  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {
    
    if(itr->second.basis->supportsCurl()) {

      // is curl required for this variable
      if(!itr->second.curl.first) 
        continue; // its not required, quit the loop

      const std::string dof_name =      itr->first;
      const std::string dof_curl_name = itr->second.curl.second;

      ParameterList p;
      p.set("Name", dof_name);
      p.set("Curl Name", dof_curl_name);
      p.set("Basis", fl.lookupLayout(dof_name)); 
      p.set("IR", ir);

      // this will help accelerate the DOFCurl evaluator when Jacobians are needed
      if(globalIndexer!=Teuchos::null) {
        // build the offsets for this field
        int fieldNum = globalIndexer->getFieldNum(dof_name);
        RCP<const std::vector<int> > offsets = 
            rcp(new std::vector<int>(globalIndexer->getGIDFieldOffsets(m_block_id,fieldNum)));
        p.set("Jacobian Offsets Vector", offsets);
      }
      // else default to the slow DOF call
    
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
        rcp(new panzer::DOFCurl<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

  }

  // Div of DOFs: Vector value @ basis --> Scalar value @ IP

  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {

    if(itr->second.basis->supportsDiv()) {

      // is div required for this variable
      if(!itr->second.div.first) 
        continue; // its not required, quit the loop

      const std::string dof_name =      itr->first;
      const std::string dof_div_name = itr->second.div.second;

      ParameterList p;
      p.set("Name", dof_name);
      p.set("Div Name", dof_div_name);
      p.set("Basis", fl.lookupLayout(dof_name)); 
      p.set("IR", ir);

      // this will help accelerate the DOFDiv evaluator when Jacobians are needed
      if(globalIndexer!=Teuchos::null) {
        // build the offsets for this field
        int fieldNum = globalIndexer->getFieldNum(dof_name);
        RCP<const std::vector<int> > offsets = 
            rcp(new std::vector<int>(globalIndexer->getGIDFieldOffsets(m_block_id,fieldNum)));
        p.set("Jacobian Offsets Vector", offsets);
      }
      // else default to the slow DOF call
    
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
        rcp(new panzer::DOFDiv<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

  // Time derivative of DOFs: Scalar value @ basis --> Scalar value @ IP 
  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {
    // is td required for this variable
    if(!itr->second.timeDerivative.first) 
      continue; // its not required, quit the loop

    const std::string td_name = itr->second.timeDerivative.second;

    ParameterList p;
    p.set("Name", td_name);
    p.set("Basis", fl.lookupLayout(itr->first)); 
    p.set("IR", ir);

    if(globalIndexer!=Teuchos::null) {
      // build the offsets for this field
      int fieldNum = globalIndexer->getFieldNum(itr->first);
      RCP<const std::vector<int> > offsets = 
          rcp(new std::vector<int>(globalIndexer->getGIDFieldOffsets(m_block_id,fieldNum)));
      p.set("Jacobian Offsets Vector", offsets);
    }
    // else default to the slow DOF call

    // set the orientiation field name explicitly if orientations are
    // required for the basis
    if(itr->second.basis->requiresOrientations())
      p.set("Orientation Field Name", itr->first+" Orientation");
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                  const panzer::FieldLibrary& /* fl */,
                                  const LinearObjFactory<panzer::Traits> & lof,
                                  const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  // this turns off the scatter contribution, and does
  // only the gather
  bool ignoreScatter = false;
  if(user_data.isParameter("Ignore Scatter")) 
    ignoreScatter = user_data.get<bool>("Ignore Scatter");
  
  if(!ignoreScatter) {
    
    for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
        itr!=m_provided_dofs_desc.end();++itr) {
      
      RCP<std::map<std::string,std::string> > names_map = rcp(new std::map<std::string,std::string>);
      RCP< std::vector<std::string> > residual_names = rcp(new std::vector<std::string>);
      
      // sanity check to make sure a residual name was registered for each provided variable
      TEUCHOS_ASSERT(itr->second.residualName.first);
      
      names_map->insert(std::make_pair(itr->second.residualName.second,itr->first));
      residual_names->push_back(itr->second.residualName.second);

      {
        ParameterList p("Scatter");
        p.set("Scatter Name", itr->second.scatterName);
        p.set("Basis", itr->second.basis.getConst());
        p.set("Dependent Names", residual_names);
        p.set("Dependent Map", names_map);
        
        RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatter<EvalT>(p);
        
        this->template registerEvaluator<EvalT>(fm, op);
      }
      
      // Require variables
      {
        PHX::Tag<typename EvalT::ScalarT> tag(itr->second.scatterName, 
                                              Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
        fm.template requireField<EvalT>(tag);
      }
    
    }

  }

}

// ***********************************************************************

template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                       const panzer::FieldLayoutLibrary& fl,
                                       const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                       const Teuchos::ParameterList& models,
                                       const Teuchos::ParameterList& user_data) const
{
  for (std::vector<std::string>::const_iterator model_name = m_closure_model_ids.begin();
       model_name != m_closure_model_ids.end(); ++model_name) {
    
    this->buildAndRegisterClosureModelEvaluators(fm,fl,ir,factory,*model_name,models,user_data);
  }
}

// ***********************************************************************

template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                       const panzer::FieldLayoutLibrary& fl,
                                       const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                       const std::string& model_name,
                                       const Teuchos::ParameterList& models,
                                       const Teuchos::ParameterList& user_data) const
{
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
    factory.getAsObject<EvalT>()->buildClosureModels(model_name,
                                                     models,
                                                     fl,
                                                     ir,
                                                     *(this->m_eval_plist),
                                                     user_data,
                                                     this->getGlobalData(),
                                                     fm);
  
  for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
    this->template registerEvaluator<EvalT>(fm, (*evaluators)[i]);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                           const panzer::FieldLibrary& fl,
                                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                           const std::string& model_name,
                                           const Teuchos::ParameterList& models,
                                           const LinearObjFactory<panzer::Traits> & lof,
                                           const Teuchos::ParameterList& user_data) const
{
  // add basis coordinates
  for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator basis =  m_unique_bases.begin();
       basis != m_unique_bases.end(); ++ basis) {
    Teuchos::RCP< PHX::Evaluator<panzer::Traits> > basis_op
      = Teuchos::rcp(new panzer::GatherBasisCoordinates<EvalT,panzer::Traits>(*basis->second));
    this->template registerEvaluator<EvalT>(fm, basis_op);
  }

  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {
    
    Teuchos::ParameterList p("Scatter");
    p.set("Scatter Name", itr->second.scatterName);
    p.set("Basis", itr->second.basis.getConst());
    Teuchos::RCP<std::vector<std::string> > name = Teuchos::rcp(new std::vector<std::string>);
    name->push_back(itr->first);
    p.set("Dependent Names", name);

    // Create an identity map
    Teuchos::RCP<std::map<std::string,std::string> > names_map = Teuchos::rcp(new std::map<std::string,std::string>);
    names_map->insert(std::make_pair(itr->first,itr->first));
    p.set("Dependent Map", names_map);

    // Set flag for ScatterDirichlet evaluators
    p.set("Scatter Initial Condition", true);

    // Use ScatterDirichlet to scatter the initial condition
    Teuchos::RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatterDirichlet<EvalT>(p);
    
    this->template registerEvaluator<EvalT>(fm, op);


    // Require field
    PHX::Tag<typename EvalT::ScalarT> tag(itr->second.scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
    fm.template requireField<EvalT>(tag);
  }

  // Add in closure models.  This is a hack that we should revisit.
  {
    // Closure models are normally evaluated at integration points,
    // but some evaluator models support evaluation at both basis and
    // integration points.  For initial guesses, we should only
    // evaluate at basis points, so integration rule is meaningless.
    // We use this to build all closure model evaluators in model
    // (including integration point based ones that will never be
    // used).  In the future we may need ir for using L2 projection to
    // basis points for initial guesses (for non-nodal bases).
    Teuchos::RCP<panzer::IntegrationRule> dummy_ir;
    if (m_int_rules.size() > 0)
      dummy_ir = m_int_rules.begin()->second;
    Teuchos::RCP<const panzer::FieldLayoutLibrary> fll = fl.buildFieldLayoutLibrary(*dummy_ir);

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
      factory.getAsObject<EvalT>()->buildClosureModels(model_name, models, *fll, dummy_ir, *(this->m_eval_plist), user_data, this->getGlobalData(), fm);
    
    for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
      this->template registerEvaluator<EvalT>(fm, (*evaluators)[i]);
  }

  // **************************
  // Add Orientation terms
  // **************************

  for (BasisIterator basis_it = m_basis_to_dofs.begin(); basis_it != m_basis_to_dofs.end(); ++basis_it) {
    if(basis_it->second.first->requiresOrientations())  {
      Teuchos::ParameterList p("Gather Orientation");
      p.set("Basis", basis_it->second.first);
      p.set("DOF Names", basis_it->second.second);
      p.set("Indexer Names", basis_it->second.second);
     
      Teuchos::RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGatherOrientation<EvalT>(p);
      
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

}

// ***********************************************************************
template <typename EvalT>
const Teuchos::RCP<Teuchos::ParameterList>
panzer::EquationSet_DefaultImpl<EvalT>::getEvaluatorParameterList() const
{
  return m_eval_plist;
}

// ***********************************************************************
template <typename EvalT>
const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >&
panzer::EquationSet_DefaultImpl<EvalT>::getProvidedDOFs() const
{
  return m_provided_dofs;
}

// ***********************************************************************
template <typename EvalT>
const std::vector<std::vector<std::string> > &
panzer::EquationSet_DefaultImpl<EvalT>::getCoordinateDOFs() const
{
  return m_coordinate_dofs;
}

// ***********************************************************************
template <typename EvalT>
const std::map<int,Teuchos::RCP<panzer::IntegrationRule> > &
panzer::EquationSet_DefaultImpl<EvalT>::getIntegrationRules() const
{
  return m_int_rules;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setElementBlockId(const std::string & blockId)
{
   TEUCHOS_ASSERT(m_block_id=="");
   m_block_id = blockId;
   this->m_eval_plist->set("Block ID", getElementBlockId());  // set the value in parameter list
                                                              // used by closure model factory
}

// ***********************************************************************
template <typename EvalT>
std::string panzer::EquationSet_DefaultImpl<EvalT>::
getElementBlockId() const
{
   return m_block_id;
}

// ***********************************************************************
template <typename EvalT>
std::string panzer::EquationSet_DefaultImpl<EvalT>::getType() const
{
  return m_type;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::setTangentParamNames(const std::vector<std::string>& tangent_param_names)
{
  m_tangent_param_names = tangent_param_names;
}

// ***********************************************************************
template <typename EvalT>
bool panzer::EquationSet_DefaultImpl<EvalT>::buildTransientSupport() const
{
  return m_build_transient_support;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
getAddedDOFs(std::vector<std::string> & dofNames) const
{
  dofNames.clear();
  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr)
    dofNames.push_back(itr->first);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
updateDOF(const std::string & dofName,
          int basisOrder,
          int integrationOrder)
{
  typename std::map<std::string,DOFDescriptor>::const_iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::updateDOF: DOF \"" << dofName << "\" has not been specified "
                             "by derived equation set \"" << this->getType() << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  desc.basisOrder = basisOrder;
  desc.basis = Teuchos::rcp(new panzer::PureBasis(desc.basisType,desc.basisOrder,m_cell_data));

  if (integrationOrder == -1)
    desc.integrationOrder = m_default_integration_order;
  else
    desc.integrationOrder = integrationOrder;

  desc.intRule = Teuchos::rcp(new panzer::IntegrationRule(desc.integrationOrder,m_cell_data));
}

// ***********************************************************************
template <typename EvalT>
int panzer::EquationSet_DefaultImpl<EvalT>::
getBasisOrder(const std::string & dofName) const
{
  typename std::map<std::string,DOFDescriptor>::const_iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::getBasisOrder: DOF \"" << dofName << "\" has not been specified "
                             "by derived equation set \"" << this->getType() << "\".");

  return itr->second.basisOrder;
}

// ***********************************************************************
template <typename EvalT>
int panzer::EquationSet_DefaultImpl<EvalT>::
getIntegrationOrder(const std::string & dofName) const
{
  typename std::map<std::string,DOFDescriptor>::const_iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::getIntegrationOrder: DOF \"" << dofName << "\" has not been specified "
                             "by derived equation set \"" << this->getType() << "\".");

  return itr->second.integrationOrder;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addDOF(const std::string & dofName,
       const std::string & basisType,
       const int & basisOrder,
       const int integrationOrder,
       const std::string residualName,
       const std::string scatterName)
{
  typename std::map<std::string,DOFDescriptor>::const_iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr!=m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addProvidedDOF: DOF \"" << dofName << "\" was previously specified "
                             "by derived equation set \"" << this->getType() << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  desc.dofName = dofName;
  desc.basisType = basisType;
  desc.basisOrder = basisOrder;
  desc.basis = Teuchos::rcp(new panzer::PureBasis(desc.basisType,desc.basisOrder,m_cell_data));

  if (integrationOrder == -1)
    desc.integrationOrder = m_default_integration_order;
  else
    desc.integrationOrder = integrationOrder;

  desc.intRule = Teuchos::rcp(new panzer::IntegrationRule(desc.integrationOrder,m_cell_data));

  // this function always creates a residual and scatter
  desc.residualName.first = true;

  if (residualName == "")
    desc.residualName.second = "RESIDUAL_" + dofName;
  else
    desc.residualName.second = residualName;

  if (scatterName == "")
    desc.scatterName = "SCATTER_" + dofName;
  else
    desc.scatterName = scatterName;

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addDOFGrad(const std::string & dofName,
           const std::string & gradName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addDOFGrad: DOF \"" << dofName << "\" has not been specified as a DOF "
                             "by derived equation set \"" << this->getType() << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check

  if (gradName == "")
    desc.grad = std::make_pair(true,std::string("GRAD_")+dofName);
  else
    desc.grad = std::make_pair(true,gradName);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addDOFCurl(const std::string & dofName,
           const std::string & curlName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addDOFCurl: DOF \"" << dofName << "\" has not been specified as a DOF "
                             "by derived equation set \"" << this->getType() << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check

  if (curlName == "")
    desc.curl = std::make_pair(true,std::string("CURL_")+dofName);
  else
    desc.curl = std::make_pair(true,curlName);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addDOFDiv(const std::string & dofName,
          const std::string & divName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addDOFDiv: DOF \"" << dofName << "\" has not been specified as a DOF "
                             "by derived equation set \"" << this->getType() << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check

  if (divName == "")
    desc.div = std::make_pair(true,std::string("DIV_")+dofName);
  else
    desc.div = std::make_pair(true,divName);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addDOFTimeDerivative(const std::string & dofName,
                              const std::string & dotName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addDOFTimeDerivative: DOF \"" << dofName << "\" has not been specified as a DOF "
                             "by derived equation set \"" << this->getType() << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check

  if (dotName == "")
    desc.timeDerivative = std::make_pair(true,std::string("DXDT_")+dofName);
  else
    desc.timeDerivative = std::make_pair(true,dotName);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setCoordinateDOFs(const std::vector<std::string> & dofNames)
{
  TEUCHOS_TEST_FOR_EXCEPTION(m_cell_data.baseCellDimension()!=Teuchos::as<int>(dofNames.size()),std::invalid_argument,
                             "EquationSet_DefaultImpl::setCoordinateDOFs: Size of vector is not equal to the "
                             "spatial dimension.");

  for(std::size_t d=0;d<dofNames.size();d++) {
    typename std::map<std::string,DOFDescriptor>::const_iterator desc_it = m_provided_dofs_desc.find(dofNames[d]);
    TEUCHOS_TEST_FOR_EXCEPTION(desc_it == m_provided_dofs_desc.end(),std::invalid_argument,
                             "EquationSet_DefaultImpl::setCoordinateDOFs: DOF of name \"" + dofNames[d] + "\" "
                             "has not been added, thus cannot be set as a coordinate DOF.");
  }

  m_coordinate_dofs.push_back(dofNames);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setDefaultValidParameters(Teuchos::ParameterList& valid_parameters)
{
  std::string default_type = "";
  valid_parameters.set("Type",default_type,"The equation set type. This must corespond to the type keyword used to build the equation set in the equation set factory.");
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<panzer::PureBasis> 
panzer::EquationSet_DefaultImpl<EvalT>::getBasisForDOF(const std::string& dof_name) const
{
  typename std::map<std::string,DOFDescriptor>::const_iterator desc_it = m_provided_dofs_desc.find(dof_name);
  TEUCHOS_ASSERT(desc_it != m_provided_dofs_desc.end());
  return desc_it->second.basis;
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<panzer::IntegrationRule> 
panzer::EquationSet_DefaultImpl<EvalT>::getIntRuleForDOF(const std::string& dof_name) const
{
  typename std::map<std::string,DOFDescriptor>::const_iterator desc_it = m_provided_dofs_desc.find(dof_name);
  TEUCHOS_TEST_FOR_EXCEPTION(desc_it == m_provided_dofs_desc.end(),std::logic_error,
                             "EquationSet_DefaultImpl::getIntRuleForDOF: Failed to find degree of freedom "
                             "with name \"" << dof_name << "\".");
  return desc_it->second.intRule;
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<panzer::BasisIRLayout> 
panzer::EquationSet_DefaultImpl<EvalT>::getBasisIRLayoutForDOF(const std::string& dof_name) const
{
  typename std::map<std::string,DOFDescriptor>::const_iterator desc_it = m_provided_dofs_desc.find(dof_name);
  TEUCHOS_TEST_FOR_EXCEPTION(desc_it == m_provided_dofs_desc.end(),std::logic_error,
                             "EquationSet_DefaultImpl::getBasisIRLayoutForDOF: Failed to find degree of freedom "
                             "with name \"" << dof_name << "\".");

  return panzer::basisIRLayout(desc_it->second.basis,*(desc_it->second.intRule));
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterResidualSummationEvaluator(PHX::FieldManager<panzer::Traits>& fm,
                                           const std::string dof_name,
                                           const std::vector<std::string>& residual_contributions,
                                           const std::string residual_field_name) const
{
  using Teuchos::rcp;
  using Teuchos::RCP;

  Teuchos::ParameterList p;

  if (residual_field_name != "")
    p.set("Sum Name", residual_field_name);
  else
    p.set("Sum Name", "RESIDUAL_" + dof_name);
  
  RCP<std::vector<std::string> > rcp_residual_contributions = rcp(new std::vector<std::string>);
  *rcp_residual_contributions = residual_contributions;

  p.set("Values Names", rcp_residual_contributions);
  
  DescriptorIterator desc_it = m_provided_dofs_desc.find(dof_name);
  TEUCHOS_ASSERT(desc_it != m_provided_dofs_desc.end());
  
  p.set("Data Layout", desc_it->second.basis->functional);

  RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::SumStatic<EvalT,panzer::Traits,panzer::Cell,panzer::BASIS>(p));
  
  this->template registerEvaluator<EvalT>(fm, op);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterResidualSummationEvaluator(PHX::FieldManager<panzer::Traits>& fm,
                                           const std::string dof_name,
                                           const std::vector<std::string>& residual_contributions,
                                           const std::vector<double>& scale_contributions,
                                           const std::string residual_field_name) const
{
  using Teuchos::rcp;
  using Teuchos::RCP;

  Teuchos::ParameterList p;

  if (residual_field_name != "")
    p.set("Sum Name", residual_field_name);
  else
    p.set("Sum Name", "RESIDUAL_" + dof_name);
  
  RCP<std::vector<std::string> > rcp_residual_contributions = rcp(new std::vector<std::string>);
  *rcp_residual_contributions = residual_contributions;
  p.set("Values Names", rcp_residual_contributions);

  RCP<const std::vector<double> > rcp_scale_contributions = rcp(new std::vector<double>(scale_contributions));
  p.set("Scalars", rcp_scale_contributions);
  
  DescriptorIterator desc_it = m_provided_dofs_desc.find(dof_name);
  TEUCHOS_ASSERT(desc_it != m_provided_dofs_desc.end());
  
  p.set("Data Layout", desc_it->second.basis->functional);

  RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::SumStatic<EvalT,panzer::Traits,panzer::Cell,panzer::BASIS>(p));
  
  this->template registerEvaluator<EvalT>(fm, op);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::addClosureModel(const std::string& closure_model)
{
  m_closure_model_ids.push_back(closure_model);
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<Teuchos::ParameterList> 
panzer::EquationSet_DefaultImpl<EvalT>::getEquationSetParameterList() const
{
  return m_input_params;
}

// ***********************************************************************
#endif
