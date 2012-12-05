// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EQUATIONSET_DEFAULT_IMPL_IMPL_HPP
#define PANZER_EQUATIONSET_DEFAULT_IMPL_IMPL_HPP

#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DOFCurl.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_GatherIntegrationCoordinates.hpp"
#include "Panzer_GatherOrientation.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_RCP.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(const panzer::InputEquationSet& ies,
			const panzer::CellData& cell_data,
			const Teuchos::RCP<panzer::GlobalData>& global_data,
			const bool build_transient_support) :
  panzer::GlobalDataAcceptorDefaultImpl(global_data),
  m_input_eq_set(ies),
  m_cell_data(cell_data),
  m_build_transient_support(build_transient_support),
  m_eqset_prefix("")
{ 
  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_gradient_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_curl_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_time_derivative_names = Teuchos::rcp(new std::vector<std::string>);
  m_residual_names = Teuchos::rcp(new std::vector<std::string>);
  m_eval_plist = Teuchos::rcp(new Teuchos::ParameterList);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setupDOFs(int equation_dimension)
{
  if(m_provided_dofs_desc.size()>0) {

    // this conditional and following assert is needed to check that the new
    // approach is not interfering with the old approach of specifying DOFs.
    bool backward_compatibility_check
           = (m_dof_names->size()==0) && 
             (m_dof_gradient_names->size()==0) && 
             (m_dof_curl_names->size()==0) && 
             (m_dof_time_derivative_names->size()==0) && 
             (m_residual_names->size()==0);
 
    TEUCHOS_TEST_FOR_EXCEPTION(!backward_compatibility_check,std::runtime_error,
                               "EquationSet_DefaultImpl::setupDOFs: You cannot call both addProvidedDOFs (and it's ilk) and set m_dof_names (for example). "
                               "The later approach to using the default implementation of the equation set is now deprecated."); 
  }
  else {
    // register with the "new" equation set way the various fields
    for(std::size_t i=0;i<m_dof_names->size();i++) {
      if(m_residual_names->size()>0)        
        addProvidedDOF((*m_dof_names)[i],(*m_residual_names)[i]);
      else
        addProvidedDOF((*m_dof_names)[i]);

      if(m_dof_gradient_names->size()>0)        
        addDOFGrad((*m_dof_names)[i],(*m_dof_gradient_names)[i]);

      if(m_dof_curl_names->size()>0)            
        addDOFCurl((*m_dof_names)[i],(*m_dof_curl_names)[i]);

      if(m_dof_time_derivative_names->size()>0) 
        addDOFTimeDerivative((*m_dof_names)[i],(*m_dof_time_derivative_names)[i]);
    }
  }

  // for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
  //     itr!=m_provided_dofs_desc.end();++itr) {
  //   itr->second.print(std::cout); std::cout << std::endl;
  // }

  this->m_int_rule = 
    Teuchos::rcp(new panzer::IntegrationRule(m_input_eq_set.integration_order,
					     m_cell_data));
  
  this->m_pure_basis = Teuchos::rcp(new panzer::PureBasis(m_input_eq_set.basis,m_cell_data));
  this->m_basis = Teuchos::rcp(new panzer::BasisIRLayout(m_pure_basis,
						 *(this->m_int_rule)));
  
  this->m_provided_dofs.clear();
  this->m_dof_names->clear();

  // load up the provided dofs from the descriptor map
  for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end();++itr) {
    this->m_dof_names->push_back(itr->first); // this contains all the dof names provided by this equation set
    this->m_provided_dofs.push_back(std::make_pair(itr->first, this->m_pure_basis));
  }

  this->m_eval_plist->set("IR", this->m_int_rule);
  this->m_eval_plist->set("Basis", this->m_basis);
  this->m_eval_plist->set("Equation Dimension", equation_dimension);
  this->m_eval_plist->set("DOF Names", this->m_dof_names);  
  this->m_eval_plist->set("Block ID", getElementBlockId());  
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
					       const LinearObjFactory<panzer::Traits> & lof,
					       const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // ********************
  // DOFs (unknowns)
  // ********************

  // Gather, includes construction of orientation gathers
  {
    ParameterList p("Gather");
    p.set("Basis", m_pure_basis);
    p.set("DOF Names", m_dof_names);
    p.set("Indexer Names", m_dof_names);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // **************************
  // Coordinates for integration points and basis functions
  // **************************
  {
    // add basis coordinates
    RCP< PHX::Evaluator<panzer::Traits> > basis_op
       = rcp(new panzer::GatherBasisCoordinates<EvalT,panzer::Traits>(*this->m_basis->getBasis()));
    fm.template registerEvaluator<EvalT>(basis_op);

    // add integration coordinates
    RCP< PHX::Evaluator<panzer::Traits> > quad_op
       = rcp(new panzer::GatherIntegrationCoordinates<EvalT,panzer::Traits>(*this->m_int_rule));
    fm.template registerEvaluator<EvalT>(quad_op);

    // NOTE: You can look up the name of either coordinate field name by doing
    //       GatherBasisCoordinates<EvalT,Traits>::fieldName();
    //       GatherIntegrationCoordinates<EvalT,Traits>::fieldName();
  }

  // **************************
  // Time derivative terms
  // **************************

  // Gather of time derivative terms
  {
    RCP< std::vector<std::string> > dof_names = rcp(new std::vector<std::string>);   // time derivative indexer names
    RCP< std::vector<std::string> > field_names = rcp(new std::vector<std::string>); // time derivative field names

    // determine which fields need time derivatives
    for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
        itr!=m_provided_dofs_desc.end();++itr) {
      std::string dofName = itr->first;
      const DOFDescriptor & desc = itr->second;

      // does this field need a time derivative?
      if(desc.timeDerivative.first) {
        // time derivaitive needed
        dof_names->push_back(dofName);
        field_names->push_back(desc.timeDerivative.second);
      }
    }

    ParameterList p("Gather");
    p.set("Basis", m_pure_basis);
    p.set("DOF Names", field_names);
    p.set("Indexer Names", dof_names);
    p.set("Use Time Derivative Solution Vector", true);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // **************************
  // Orientation terms
  // **************************

  if(this->m_pure_basis->requiresOrientations())  {
    ParameterList p("Gather Orientation");
    p.set("Basis", m_pure_basis);
    p.set("DOF Names", m_dof_names);
    p.set("Indexer Names", m_dof_names);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGatherOrientation<EvalT>(p);

    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					     const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
					     const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  // DOFs: Scalar value @ basis --> Scalar value @ IP 
  for (std::vector<std::string>::const_iterator dof_name = this->m_dof_names->begin();
       dof_name != this->m_dof_names->end(); ++dof_name) {
    ParameterList p;
    p.set("Name", *dof_name);
    p.set("Basis", this->m_basis); 
    p.set("IR", this->m_int_rule);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Gradients of DOFs: Scalar value @ basis --> Vector value @ IP
  if(this->m_pure_basis->supportsGrad()) {

    for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
        itr!=m_provided_dofs_desc.end();++itr) {
      // is gradient required for this variable
      if(!itr->second.grad.first) 
        continue; // its not required, quit the loop

      const std::string dof_name =      itr->first;
      const std::string dof_grad_name = itr->second.grad.second;

      ParameterList p;
      p.set("Name", dof_name);
      p.set("Gradient Name", dof_grad_name);
      p.set("Basis", this->m_basis); 
      p.set("IR", this->m_int_rule);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Curl of DOFs: Vector value @ basis --> Vector value @ IP (3D) or Scalar value @ IP (2D)
  if(this->m_pure_basis->supportsCurl()) {

    for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
        itr!=m_provided_dofs_desc.end();++itr) {
      // is curl required for this variable
      if(!itr->second.curl.first) 
        continue; // its not required, quit the loop

      const std::string dof_name =      itr->first;
      const std::string dof_curl_name = itr->second.curl.second;

      ParameterList p;
      p.set("Name", dof_name);
      p.set("Curl Name", dof_curl_name);
      p.set("Basis", this->m_basis); 
      p.set("IR", this->m_int_rule);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::DOFCurl<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
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
    p.set("Basis", this->m_basis); 
    p.set("IR", this->m_int_rule);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				  const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
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
     // Scatter
     RCP<std::map<std::string,std::string> > names_map = rcp(new std::map<std::string,std::string>);
     RCP< std::vector<std::string> > residual_names = rcp(new std::vector<std::string>);
   
     for(typename std::map<std::string,DOFDescriptor>::const_iterator itr=m_provided_dofs_desc.begin();
         itr!=m_provided_dofs_desc.end();++itr) {
       // sanity check to make sure a residual name was registered for each provided variable
       TEUCHOS_ASSERT(itr->second.residualName.first);

       names_map->insert(std::make_pair(itr->second.residualName.second,itr->first));
       residual_names->push_back(itr->second.residualName.second);
     }
    
   
     {
       ParameterList p("Scatter");
       p.set("Scatter Name", this->m_scatter_name);
       p.set("Basis", this->m_pure_basis.getConst());
       p.set("Dependent Names", residual_names);
       p.set("Dependent Map", names_map);
   
       RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatter<EvalT>(p);
         // rcp(new panzer::ScatterResidual_Epetra<EvalT,panzer::Traits>(p));
       
       fm.template registerEvaluator<EvalT>(op);
     }
     
     // Require variables
     {
       PHX::Tag<typename EvalT::ScalarT> tag(this->m_scatter_name, 
   					  Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
       fm.template requireField<EvalT>(tag);
     }
  }

}

// ***********************************************************************

template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
				       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				       const Teuchos::ParameterList& models,
				       const Teuchos::ParameterList& user_data) const
{
  buildAndRegisterClosureModelEvaluators(fm,dofs,factory,this->m_input_eq_set.model_id,models,user_data);
}

// ***********************************************************************

template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				   const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
				   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				   const std::string& model_name,
				   const Teuchos::ParameterList& models,
				   const Teuchos::ParameterList& user_data) const
{
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
    factory.getAsObject<EvalT>()->buildClosureModels(model_name, this->m_input_eq_set, models, 
                                                     *(this->m_eval_plist), user_data, this->getGlobalData(), fm);
    
  for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
    fm.template registerEvaluator<EvalT>((*evaluators)[i]);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					   const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
					   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					   const std::string& model_name,
					   const Teuchos::ParameterList& models,
					   const LinearObjFactory<panzer::Traits> & lof,
					   const Teuchos::ParameterList& user_data) const
{
  // **************************
  // Coordinates for integration points and basis functions
  // **************************
  {
    // add basis coordinates
    Teuchos::RCP< PHX::Evaluator<panzer::Traits> > basis_op
       = Teuchos::rcp(new panzer::GatherBasisCoordinates<EvalT,panzer::Traits>(*this->m_basis->getBasis()));
    fm.template registerEvaluator<EvalT>(basis_op);
  }

  {
    Teuchos::RCP<const PureBasis> basis = m_field_layout_lib->lookupBasis((*m_dof_names)[0]);

    Teuchos::ParameterList p("Scatter");
    p.set("Scatter Name", this->m_scatter_name);
    p.set("Basis", basis);
    p.set("Dependent Names", this->m_dof_names);

    Teuchos::RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatterInitialCondition<EvalT>(p);
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Require variables
  {
    PHX::Tag<typename EvalT::ScalarT> tag(this->m_scatter_name, 
					  Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
    fm.template requireField<EvalT>(tag);
  }  

  // Add in closure models
  {
    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
      factory.getAsObject<EvalT>()->buildClosureModels(model_name, this->m_input_eq_set, models, *(this->m_eval_plist), user_data, this->getGlobalData(), fm);
    
    for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
      fm.template registerEvaluator<EvalT>((*evaluators)[i]);
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
const std::vector<std::string>&
panzer::EquationSet_DefaultImpl<EvalT>::getDOFNames() const
{
  return *m_dof_names;
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
Teuchos::RCP<panzer::IntegrationRule> panzer::EquationSet_DefaultImpl<EvalT>::
getIntegrationRule() const
{
   return m_int_rule;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setFieldLayoutLibrary(const panzer::FieldLibrary & fieldLibrary)
{
   m_field_layout_lib = fieldLibrary.buildFieldLayoutLibrary(*m_int_rule);

   this->m_eval_plist->set("Field Layout Library", m_field_layout_lib);  // set the value in parameter list
                                                                         // used by closure model factory
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<const panzer::FieldLayoutLibrary> panzer::EquationSet_DefaultImpl<EvalT>::
getFieldLayoutLibrary() const
{
   return m_field_layout_lib;
}

// ***********************************************************************
template <typename EvalT>
const panzer::InputEquationSet& 
panzer::EquationSet_DefaultImpl<EvalT>::getInputEquationSet() const
{
  return m_input_eq_set;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addProvidedDOF(const std::string & dofName,
               const std::string & residualName)
{
  typename std::map<std::string,DOFDescriptor>::const_iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr!=m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addProvidedDOF: DOF \"" << dofName << "\" was previously specified "
                             "by derived equation set \"" << m_scatter_name << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  desc.dofName = dofName;
  desc.residualName.first = true;
  desc.residualName.second = residualName;
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
addProvidedDOF(const std::string & dofName)
{
  typename std::map<std::string,DOFDescriptor>::const_iterator itr = m_provided_dofs_desc.find(dofName);

  TEUCHOS_TEST_FOR_EXCEPTION(itr!=m_provided_dofs_desc.end(),std::runtime_error,
                             "EquationSet_DefaultImpl::addProvidedDOF: DOF \"" << dofName << "\" was previously specified "
                             "by derived equation set \"" << m_scatter_name << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  desc.dofName = dofName;
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
                             "by derived equation set \"" << m_scatter_name << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check
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
                             "by derived equation set \"" << m_scatter_name << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check
  desc.curl = std::make_pair(true,curlName);
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
                             "by derived equation set \"" << m_scatter_name << "\".");

  // allocate and populate a dof descriptor associated with the field "dofName"
  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  TEUCHOS_ASSERT(desc.dofName==dofName); // safety check
  desc.timeDerivative = std::make_pair(true,dotName);
}

// ***********************************************************************

#endif
