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

#ifndef PANZER_EQUATION_SET_DEFAULTIMPL_DECL_HPP
#define PANZER_EQUATION_SET_DEFAULTIMPL_DECL_HPP

#include "Panzer_EquationSet.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_GlobalDataAcceptor_DefaultImpl.hpp"
#include "Panzer_CellData.hpp"

namespace PHX {
  template<typename Traits> class FieldManager;
}

namespace panzer {

  template <typename EvalT>
  class EquationSet_DefaultImpl : public panzer::EquationSet<EvalT>,
				  public panzer::GlobalDataAcceptorDefaultImpl {
    
  public:    
    
    EquationSet_DefaultImpl(const panzer::InputEquationSet& ies, const panzer::CellData& cell_data, const Teuchos::RCP<panzer::GlobalData>& global_data, const bool build_transient_support);
    
    virtual ~EquationSet_DefaultImpl() {}
    
    //! Builds the integration rule, basis, DOFs, and default parameter list
    virtual void setupDOFs(int equation_dimension);

    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
						       const Teuchos::ParameterList& user_data) const = 0;

    virtual void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
								const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
								const LinearObjFactory<panzer::Traits> & lof,
								const Teuchos::ParameterList& user_data) const;
    
    virtual void buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							      const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
							      const Teuchos::ParameterList& user_data) const;
    
    virtual void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						   const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
						   const LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
							const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							const Teuchos::ParameterList& models,
							const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
							const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                        const std::string & model_name,
							const Teuchos::ParameterList& models,
							const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							    const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
							    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							    const std::string& model_name,
							    const Teuchos::ParameterList& models,
							    const LinearObjFactory<panzer::Traits> & lof,
							    const Teuchos::ParameterList& user_data) const;

    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const;
    
    virtual const std::vector<std::string> & getDOFNames() const;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > > & getProvidedDOFs() const;

    void setElementBlockId(const std::string & blockId);

    std::string getElementBlockId() const;

    virtual Teuchos::RCP<panzer::IntegrationRule> getIntegrationRule() const;

    virtual void setFieldLayoutLibrary(const FieldLibrary & fieldLibrary);
    Teuchos::RCP<const FieldLayoutLibrary> getFieldLayoutLibrary() const;

    const panzer::InputEquationSet& getInputEquationSet() const;

  protected:
    
    const panzer::InputEquationSet m_input_eq_set;
    const panzer::CellData m_cell_data;
    const bool m_build_transient_support;
    
    std::string m_eqset_prefix;
    
    Teuchos::RCP<panzer::IntegrationRule> m_int_rule;
    Teuchos::RCP<panzer::PureBasis> m_pure_basis;
    Teuchos::RCP<panzer::BasisIRLayout> m_basis;
    
    std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >  m_provided_dofs;
    Teuchos::RCP< std::vector<std::string> > m_dof_names;
    Teuchos::RCP< std::vector<std::string> > m_dof_gradient_names;
    Teuchos::RCP< std::vector<std::string> > m_dof_curl_names;
    Teuchos::RCP< std::vector<std::string> > m_dof_time_derivative_names;
    Teuchos::RCP< std::vector<std::string> > m_residual_names;
    std::string m_scatter_name;
    Teuchos::RCP<Teuchos::ParameterList> m_eval_plist;
    Teuchos::RCP<const FieldLayoutLibrary> m_field_layout_lib;

    std::string m_block_id;
  };
  
}

#endif
