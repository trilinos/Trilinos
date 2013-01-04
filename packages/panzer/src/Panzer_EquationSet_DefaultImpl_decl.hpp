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

#include <map>

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

    //! Returns true if transient support should be enabled in the equation set
    bool buildTransientSupport() const;

    // The set of functions below are for use by derived classes to specify the 
    // provided degree of freedom (and associated residual name), in addition
    // to enabling the, gradient, curl and time derivative for those.

    /** Alert the panzer library of a DOF provided by this equation set.
      * This automatically sets up the gather/scatter routines neccessary
      * to evaluate and assemble with this unknown.
      *
      * \param[in] dofName Name of field to lookup in the unique global
      *                        indexer. This also serves as a key for the remaining
      *                        <code>addDOF*</code> methods.
      * \param[in] residualName Name of field that is to be scattered associated with
      *                         this DOF.
      */
    void addProvidedDOF(const std::string & dofName,
                        const std::string & residualName);

    /** Alert the panzer library of a DOF provided by this equation set.
      * This version of the method does not sets up the scatter routines
      * because there is no specified residual name.
      *
      * \param[in] dofName Name of field to lookup in the unique global
      *                        indexer. This also serves as a key for the remaining
      *                        <code>addDOF*</code> methods.
      */
    void addProvidedDOF(const std::string & dofName);

    /** Alert the panzer library that a gradient of particular a DOF is needed.
      *
      * \param[in] dofName Name of field to lookup in the unique global indexer. 
      * \param[in] gradName Name of the gradient field associated with
      *                         this DOF.
      */
    void addDOFGrad(const std::string & dofName,
                    const std::string & gradName);

    /** Alert the panzer library that a curl of particular a DOF is needed.
      *
      * \param[in] dofName Name of field to lookup in the unique global indexer. 
      * \param[in] curlName Name of the curl field associated with
      *                         this DOF.
      */
    void addDOFCurl(const std::string & dofName,
                    const std::string & curlName);

    /** Alert the panzer library that a time derivative of particular a DOF is needed.
      *
      * \param[in] dofName Name of field to lookup in the unique global indexer. 
      * \param[in] dotName Name of the time derivative field associated with
      *                         this DOF.
      */
    void addDOFTimeDerivative(const std::string & dofName,
                              const std::string & dotName);

    struct DOFDescriptor {
      DOFDescriptor() 
        : dofName("")
        , residualName(std::make_pair(false,""))
        , grad(std::make_pair(false,""))
        , curl(std::make_pair(false,""))
        , timeDerivative(std::make_pair(false,"")) {}

      std::string dofName;
      std::pair<bool,std::string> residualName;
      std::pair<bool,std::string> grad;
      std::pair<bool,std::string> curl;
      std::pair<bool,std::string> timeDerivative;

      void print(std::ostream & os) const {
        os << "DOF Desc = \"" << dofName << "\": "
           << "Res = (" << residualName.first << ", \"" << residualName.second << "\"), "
           << "Grad = (" << grad.first << ", \"" << grad.second << "\"), "
           << "Curl = (" << curl.first << ", \"" << curl.second << "\"), "
           << "Time = (" << timeDerivative.first << ", \"" << timeDerivative.second << "\")";
      }
    };

    std::map<std::string,DOFDescriptor> m_provided_dofs_desc;
    

    const panzer::InputEquationSet m_input_eq_set;
    const panzer::CellData m_cell_data;
    const bool m_build_transient_support;
    
    std::string m_eqset_prefix;
    
    Teuchos::RCP<panzer::IntegrationRule> m_int_rule;
    Teuchos::RCP<panzer::PureBasis> m_pure_basis;
    Teuchos::RCP<panzer::BasisIRLayout> m_basis;
    
    std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >  m_provided_dofs;
    Teuchos::RCP< std::vector<std::string> > m_dof_names;
    std::string m_scatter_name;
    Teuchos::RCP<Teuchos::ParameterList> m_eval_plist;
    Teuchos::RCP<const FieldLayoutLibrary> m_field_layout_lib;

    std::string m_block_id;
  };
  
}

#endif
