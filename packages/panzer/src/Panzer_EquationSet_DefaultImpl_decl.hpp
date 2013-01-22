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
#include "Panzer_GlobalDataAcceptor_DefaultImpl.hpp"
#include "Panzer_CellData.hpp"

#include <map>

namespace Teuchos {
  class ParameterList;
}

namespace PHX {
  template<typename Traits> class FieldManager;
}

namespace panzer {

  template <typename EvalT>
  class EquationSet_DefaultImpl : public panzer::EquationSet<EvalT>,
				  public panzer::GlobalDataAcceptorDefaultImpl {
    
  public:    
    
    EquationSet_DefaultImpl(const Teuchos::RCP<Teuchos::ParameterList>& params,
			    const int& default_integration_order,
			    const panzer::CellData& cell_data,
			    const Teuchos::RCP<panzer::GlobalData>& global_data,
			    const bool build_transient_support);
    
    virtual ~EquationSet_DefaultImpl() {}
    
    virtual void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
								const panzer::FieldLibrary& fl,
								const LinearObjFactory<panzer::Traits> & lof,
								const Teuchos::ParameterList& user_data) const;
    
    virtual void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						   const panzer::FieldLibrary& fl,
						   const LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						       const panzer::FieldLibrary& fl,
						       const Teuchos::ParameterList& user_data) const = 0;

    virtual void buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							      const panzer::FieldLayoutLibrary& fl,
							      const Teuchos::RCP<panzer::IntegrationRule>& ir,
							      const Teuchos::ParameterList& user_data) const;
    
    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const panzer::FieldLayoutLibrary& fl,
							const Teuchos::RCP<panzer::IntegrationRule>& ir,
							const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							const Teuchos::ParameterList& models,
							const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const panzer::FieldLayoutLibrary& fl,
							const Teuchos::RCP<panzer::IntegrationRule>& ir,
							const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                        const std::string & model_name,
							const Teuchos::ParameterList& models,
							const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							    const panzer::FieldLibrary& fl,
							    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							    const std::string& model_name,
							    const Teuchos::ParameterList& models,
							    const LinearObjFactory<panzer::Traits> & lof,
							    const Teuchos::ParameterList& user_data) const;

    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > > & getProvidedDOFs() const;

    virtual const std::map<int,Teuchos::RCP<panzer::IntegrationRule> > & getIntegrationRules() const;

    void setElementBlockId(const std::string & blockId);

    std::string getElementBlockId() const;

    virtual std::string getType() const;

    virtual std::string getKey() const;

  protected:

    //! Builds the integration rule, basis, DOFs, and default parameter list
    virtual void setupDOFs();

    //! Returns true if transient support should be enabled in the equation set
    bool buildTransientSupport() const;

    //! If the concrete equation set will be used multiple times in the same physics block, then a unique key must be set using this call from the user derived equation set constructor. If not called by the user, it defaults to the "type" and this equation set can only be instantiated once per physics block. 
    void setKey(const std::string& key);

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
      * \param[in] basisType Name of the basis type for this DOF.
      * \param[in] basisOrder Polynomial order for the basis for this DOF.
      * \param[in] integrationOrder Order of the integration rule associated with this
      *                             DOF.  If set to -1 (default), it will use the default
      *                             integration order.
      * \param[in] residualName Name of field that is to be scattered associated with
      *                         this DOF.  If not supplied or an empty string used, the
      *                         default is to add the prefix "RESIDUAL_" to the dofName for
      *                         the residual field name.
      * \param[in] scatterName Name of the required scatter field associated with
      *                        this DOF.  If not supplied or an empty string used,
      *                        the default is to add the prefix "SCATTER_"
      *                        to the dofName for the scatter field name.
      */
    void addDOF(const std::string & dofName,
		const std::string & basisType,
		const int & basisOrder,
		const int integrationOrder = -1,
		const std::string residualName = "",
		const std::string scatterName = "");

    /** Alert the panzer library that a gradient of particular a DOF is needed.
      *
      * \param[in] dofName Name of field to lookup in the unique global indexer. 
      * \param[in] gradName Name of the gradient field associated with
      *                     this DOF.  If not supplied or an empty string used,
      *                     the default is to add the prefix "GRAD_"
      *                     to the dofName for the name of the gradient field.
      */
    void addDOFGrad(const std::string & dofName,
                    const std::string & gradName = "");

    /** Alert the panzer library that a curl of particular a DOF is needed.
      *
      * \param[in] dofName Name of field to lookup in the unique global indexer. 
      * \param[in] curlName Name of the curl field associated with
      *                     this DOF.  If not supplied or an empty string used,
      *                     the default is to add the prefix "CURL_"
      *                     to the dofName for the naem of the curl field.
      */
    void addDOFCurl(const std::string & dofName,
                    const std::string & curlName = "");

    /** Alert the panzer library that a time derivative of particular a DOF is needed.
      *
      * \param[in] dofName Name of field to lookup in the unique global indexer. 
      * \param[in] dotName Name of the time derivative field associated with
      *                    this DOF.  If not supplied or an empty string used,
      *                    the default is to add the prefix "DXDT_"
      *                    to the dofName for the name of the time derivative field.
      */
    void addDOFTimeDerivative(const std::string & dofName,
                              const std::string & dotName = "");

    struct DOFDescriptor {
      DOFDescriptor() 
        : dofName("")
        , residualName(std::make_pair(false,""))
        , grad(std::make_pair(false,""))
        , curl(std::make_pair(false,""))
        , timeDerivative(std::make_pair(false,"")) {}

      std::string dofName;
      std::string basisType;
      int basisOrder;
      Teuchos::RCP<panzer::PureBasis> basis;
      int integrationOrder;
      Teuchos::RCP<panzer::IntegrationRule> intRule;
      std::pair<bool,std::string> residualName;
      std::string scatterName;
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

    //! Maps the dof name into a DOFDescriptor
    std::map<std::string,DOFDescriptor> m_provided_dofs_desc;

    //! For convenience, declare the DOFDescriptor iterator
    typedef typename std::map<std::string,DOFDescriptor>::const_iterator DescriptorIterator;


    /** \brief Map that links a common basis to a vector of dof names.  Key is the unique basis name, the value is a pair that contains an RCP to a basis and an RCP to a vector of dof names that share the basis.

        Some of our evaluators are vectorized to work on a block of
        dofs as long as they share a common basis.  We can minimize
        the evaluators built below by grouping dofs with a common
        basis.  This struct is for grouping dofs with a common basis.
    */
    std::map<std::string,std::pair<Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<std::vector<std::string> > > > m_basis_to_dofs;

    //! For convenience, declare a basis iterator
    typedef typename std::map<std::string,std::pair<Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<std::vector<std::string> > > >::const_iterator BasisIterator;
    

    const Teuchos::RCP<Teuchos::ParameterList> m_input_params;
    int m_default_integration_order;
    const panzer::CellData m_cell_data;
    const bool m_build_transient_support;
    
    //! Key is the dof name and the value is the corresponding basis
    std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >  m_provided_dofs;

    //! Key is the integration rule order and the value is the corresponding integration rule
    std::map<int,Teuchos::RCP<panzer::IntegrationRule> > m_int_rules;

    //! Key is the basis name from panzer::PureBasis::name() and value is the corresponding PureBasis
    std::map<std::string,Teuchos::RCP<panzer::PureBasis> > m_unique_bases;

    Teuchos::RCP<Teuchos::ParameterList> m_eval_plist;

    std::string m_block_id;
    std::string m_type;
    std::string m_key;
    std::string m_model_id;

    // Deprecated code support, NOTE: this assumes the same basis and inte rule are used for all dofs in the phsyics block!!!  We are setting these to avoid having to change closure model factories for all physics right away.
    void setupDeprecatedDOFsSupport();

  };
  
}

#endif
