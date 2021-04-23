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
// Questions? Contact Graham B. Harper (gbharpe@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __muelu_region_poisson_hpp__
#define __muelu_region_poisson_hpp__

// Standard headers
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

// Teuchos headers
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// Panzer headers
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_ClosureModel_Factory.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

// Phalanx headers
#include "Phalanx_config.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_MDField.hpp"

// Sacado headers
#include "Sacado_mpl_apply.hpp"

void perceptrenumbertest(const unsigned int num_levels_refinement)
{
  std::cout << "Starting perceptrenumber..." << std::endl;

  const unsigned int d = 2; // spatial dimension (only 2 or 3 makes sense)... I only really care about the 2D case since it's where the ccw rotation lies
  const unsigned int r = num_levels_refinement; // levels of recursion/refinement
  const unsigned int cells_per_dim = (1 << r); // number of cells per dim is 2^(r-1), i.e. 1<<r
  const unsigned int points_per_dim = cells_per_dim + 1; // number of points per spatial dimension is just 1 higher
  const unsigned int num_cells = (d<3)? cells_per_dim*cells_per_dim : cells_per_dim*cells_per_dim*cells_per_dim; // evil ternary op code
  const unsigned int num_points = (d<3)? points_per_dim*points_per_dim : points_per_dim*points_per_dim*points_per_dim; // evil ternary op code

  std::cout << "Parameters: " << std::endl;
  std::cout << "   Dimension = " << d << std::endl;
  std::cout << "      Levels = " << r << std::endl;
  std::cout << "   Cells/dim = " << cells_per_dim << std::endl;
  std::cout << "       Cells = " << num_cells << std::endl;
  std::cout << "  Points/dim = " << points_per_dim << std::endl;
  std::cout << "      Points = " << num_points << std::endl;
  std::cout << std::endl;
  std::cout << "Running re-ordering scheme..." << std::endl;

  unsigned int percept[2][2][2];

  // notice the numbering is ccw in the xy plane
  // assuming the indexing is [x][y][z]... this looks backwards as far as access patterns go
  percept[0][0][0] = 0;
  percept[1][0][0] = 1;
  percept[1][1][0] = 2;
  percept[0][1][0] = 3;
  percept[0][0][1] = 4;
  percept[1][0][1] = 5;
  percept[1][1][1] = 6;
  percept[0][1][1] = 7;

  // scope in case I decide to copypasta
  {
    // initialize to 0
    unsigned int outputorder[cells_per_dim][cells_per_dim][cells_per_dim];
    for(unsigned int i=0; i<cells_per_dim; ++i)
      for(unsigned int j=0; j<cells_per_dim; ++j)
        for(unsigned int k=0; k<cells_per_dim; ++k)
          outputorder[i][j][k] = 0;

    // do this for the first levels
    for(unsigned int level=0; level<(r-1); ++level)
      for(unsigned int i=0; i<cells_per_dim; ++i)
        for(unsigned int j=0; j<cells_per_dim; ++j)
          for(unsigned int k=0; k<cells_per_dim; ++k)
            outputorder[i][j][k] = outputorder[i][j][k] + (1<<(d*level))*percept[(i%(2<<level))/(1<<level)][(j%(2<<level))/(1<<level)][(k%(2<<level))/(1<<level)];

    // do this for the top level
    for(unsigned int i=0; i<cells_per_dim; ++i)
      for(unsigned int j=0; j<cells_per_dim; ++j)
        for(unsigned int k=0; k<cells_per_dim; ++k)
          outputorder[i][j][k] = outputorder[i][j][k] + (1<<(d*(r-1)))*percept[i/(cells_per_dim/2)][j/(cells_per_dim/2)][k/(cells_per_dim/2)];


    std::cout << "Outputting Percept's order in 2D... " << std::endl;
    for(int j=cells_per_dim-1; j>=0; --j)
    {
      std::cout << outputorder[0][j][0];
      for(unsigned int i=1; i<cells_per_dim; ++i)
      {
        std::cout << " " << outputorder[i][j][0];
      }
      std::cout << std::endl;
    }

    std::cout << "Outputting the lexicographic reordering..." << std::endl;
    for(unsigned int j=0; j<cells_per_dim; ++j)
      for(unsigned int i=0; i<cells_per_dim; ++i)
        std::cout << " " << outputorder[i][j][0];
    std::cout << std::endl;
  }

  std::cout << "Done!" << std::endl;
  std::cout << std::endl;
}


namespace panzer {
  class InputEquationSet;
}

namespace Example {
  
  // ******************************************************* //
  // ******************** BC STRATEGIES ******************** //
  // ******************************************************* //
  template <typename EvalT>
  class BCStrategy_Dirichlet_Constant : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {
  public:    
    
    BCStrategy_Dirichlet_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);
    
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& user_data);
    
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& pb,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				    const Teuchos::ParameterList& models,
				    const Teuchos::ParameterList& user_data) const;

    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;
  };

  // ***********************************************************************
  template <typename EvalT>
  BCStrategy_Dirichlet_Constant<EvalT>::
  BCStrategy_Dirichlet_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
    panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc,global_data)
  {
    TEUCHOS_ASSERT(this->m_bc.strategy() == "Constant");
  }

  // ***********************************************************************
  template <typename EvalT>
  void BCStrategy_Dirichlet_Constant<EvalT>::
  setup(const panzer::PhysicsBlock& side_pb,
	const Teuchos::ParameterList& /* user_data */)
  {
    using Teuchos::RCP;
    using std::vector;
    using std::string;
    using std::pair;

    // need the dof value to form the residual
    this->required_dof_names.push_back(this->m_bc.equationSetName());

    // unique residual name
    this->residual_name = "Residual_" + this->m_bc.identifier();

    // map residual to dof 
    this->residual_to_dof_names_map[residual_name] = this->m_bc.equationSetName();

    // map residual to target field
    this->residual_to_target_field_map[residual_name] = "Constant_" + this->m_bc.equationSetName();

    // find the basis for this dof 
    const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it = 
	   dofs.begin(); dof_it != dofs.end(); ++dof_it) {
      if (dof_it->first == this->m_bc.equationSetName())
	this->basis = dof_it->second;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
			       "Error the name \"" << this->m_bc.equationSetName()
			       << "\" is not a valid DOF for the boundary condition:\n"
			       << this->m_bc << "\n");

  }

  // ***********************************************************************
  template <typename EvalT>
  void BCStrategy_Dirichlet_Constant<EvalT>::
  buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			     const panzer::PhysicsBlock& /* pb */,
			     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			     const Teuchos::ParameterList& /* models */,
			     const Teuchos::ParameterList& /* user_data */) const
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // provide a constant target value to map into residual
    {
      ParameterList p("BC Constant Dirichlet");
      p.set("Name", "Constant_" + this->m_bc.equationSetName());
      p.set("Data Layout", basis->functional);
      p.set("Value", this->m_bc.params()->template get<double>("Value"));
    
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
    
      this->template registerEvaluator<EvalT>(fm, op);
    }

  }

  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(BCStrategy_Dirichlet_Constant,
					     BCStrategy_Dirichlet_Constant)

  struct BCStrategyFactory : public panzer::BCStrategyFactory {

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc,const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
	Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant",
				      BCStrategy_Dirichlet_Constant)

	TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
				   "Error - the BC Strategy called \"" << bc.strategy() <<
				    "\" is not a valid identifier in the BCStrategyFactory.  Either add a "
                         "valid implementation to your factory or fix your input file.  The "
				   "relevant boundary condition is:\n\n" << bc << std::endl);
      
      return bcs_tm;
    }

  };

  // ******************************************************* //
  // ******************** MODEL FACTORY ******************** //
  // ******************************************************* //
  template<typename EvalT>
  class ModelFactory : public panzer::ClosureModelFactory<EvalT> {

  public:

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  };

  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

  /** A source for the curl Laplacian that results in the solution
   */
  template<typename EvalT, typename Traits>
  class SimpleSource : public panzer::EvaluatorWithBaseImpl<Traits>,
		       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    SimpleSource(const std::string & name,
		 const panzer::IntegrationRule & ir);
                                                                        
    void postRegistrationSetup(typename Traits::SetupData d,           
                               PHX::FieldManager<Traits>& fm);        
                                                                     
    void evaluateFields(typename Traits::EvalData d);               


  private:
    typedef typename EvalT::ScalarT ScalarT;

    // Simulation source
    PHX::MDField<ScalarT,Cell,Point> source;
    int ir_degree, ir_index;
  };

  //**********************************************************************
  template <typename EvalT,typename Traits>
  SimpleSource<EvalT,Traits>::SimpleSource(const std::string & name,
					   const panzer::IntegrationRule & ir)
  {
    using Teuchos::RCP;

    Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
    ir_degree = ir.cubature_degree;

    source = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);

    this->addEvaluatedField(source);
  
    std::string n = "Simple Source";
    this->setName(n);
  }

  //**********************************************************************
  template <typename EvalT,typename Traits>
  void SimpleSource<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
							 PHX::FieldManager<Traits>& /* fm */)
  {
    ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
  }

  //**********************************************************************
  template <typename EvalT,typename Traits>
  void SimpleSource<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
  { 
    using panzer::index_t;
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < source.extent_int(1); ++point) {
	const double& x = workset.int_rules[ir_index]->ip_coordinates(cell,point,0);
	const double& y = workset.int_rules[ir_index]->ip_coordinates(cell,point,1);

	source(cell,point) = 8.0*M_PI*M_PI*std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y);
      }
    }
  }

  //**********************************************************************

  /** The analytic solution to the mixed poisson equation for the sine source.
   */
  template<typename EvalT, typename Traits>
  class SimpleSolution : public panzer::EvaluatorWithBaseImpl<Traits>,
			 public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    SimpleSolution(const std::string & name,
		   const panzer::IntegrationRule & ir);
                                                                        
    void postRegistrationSetup(typename Traits::SetupData d,           
                               PHX::FieldManager<Traits>& fm);        
                                                                     
    void evaluateFields(typename Traits::EvalData d);               


  private:
    typedef typename EvalT::ScalarT ScalarT;

    // Simulation solution
    PHX::MDField<ScalarT,Cell,Point> solution;
    PHX::MDField<ScalarT,Cell,Point,Dim> solution_grad;
    int ir_degree, ir_index;
  };

  //**********************************************************************
  template <typename EvalT,typename Traits>
  SimpleSolution<EvalT,Traits>::SimpleSolution(const std::string & name,
					       const panzer::IntegrationRule & ir)
  {
    using Teuchos::RCP;

    Teuchos::RCP<PHX::DataLayout> data_layout_scalar = ir.dl_scalar;
    Teuchos::RCP<PHX::DataLayout> data_layout_vector = ir.dl_vector;
    ir_degree = ir.cubature_degree;

    solution = PHX::MDField<ScalarT,Cell,Point>(name, data_layout_scalar);
    solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_"+name, data_layout_vector);

    this->addEvaluatedField(solution);
    this->addEvaluatedField(solution_grad);
  
    std::string n = "Simple Solution";
    this->setName(n);
  }

  //**********************************************************************
  template <typename EvalT,typename Traits>
  void SimpleSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
							   PHX::FieldManager<Traits>& /* fm */)
  {
    ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
  }

  //**********************************************************************
  template <typename EvalT,typename Traits>
  void SimpleSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
  { 
    using panzer::index_t;
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < solution.extent_int(1); ++point) {

	const double & x = this->wda(workset).int_rules[ir_index]->ip_coordinates(cell,point,0);
	const double & y = this->wda(workset).int_rules[ir_index]->ip_coordinates(cell,point,1);

	solution(cell,point) = std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
	solution_grad(cell,point,0) = 2.0*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
	solution_grad(cell,point,1) = 2.0*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
      }
    }
  }

  //**********************************************************************


  // ********************************************************************
  // ********************************************************************
  template<typename EvalT>
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
  ModelFactory<EvalT>::
  buildClosureModels(const std::string& model_id,
		     const Teuchos::ParameterList& models, 
		     const panzer::FieldLayoutLibrary& fl,
		     const Teuchos::RCP<panzer::IntegrationRule>& ir,
		     const Teuchos::ParameterList& /* default_params */,
		     const Teuchos::ParameterList& /* user_data */,
		     const Teuchos::RCP<panzer::GlobalData>& /* global_data */,
		     PHX::FieldManager<panzer::Traits>& /* fm */) const
  {
    using std::string;
    using std::vector;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
    using PHX::Evaluator;

    RCP< vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
      rcp(new vector< RCP<Evaluator<panzer::Traits> > > );

    if (!models.isSublist(model_id)) {
      models.print(std::cout);
      std::stringstream msg;
      msg << "Falied to find requested model, \"" << model_id 
	  << "\", for equation set:\n" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(!models.isSublist(model_id), std::logic_error, msg.str());
    }

    std::vector<Teuchos::RCP<const panzer::PureBasis> > bases;
    fl.uniqueBases(bases);

    const ParameterList& my_models = models.sublist(model_id);

    for (ParameterList::ConstIterator model_it = my_models.begin(); 
	 model_it != my_models.end(); ++model_it) {
    
      bool found = false;
    
      const std::string key = model_it->first;
      ParameterList input;
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);

      if (plist.isType<double>("Value")) {
	{ // at IP
	  input.set("Name", key);
	  input.set("Value", plist.get<double>("Value"));
	  input.set("Data Layout", ir->dl_scalar);
	  RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
	}
      
	for (std::vector<Teuchos::RCP<const panzer::PureBasis> >::const_iterator basis_itr = bases.begin();
	     basis_itr != bases.end(); ++basis_itr) { // at BASIS
	  input.set("Name", key);
	  input.set("Value", plist.get<double>("Value"));
	  Teuchos::RCP<const panzer::BasisIRLayout> basis = basisIRLayout(*basis_itr,*ir);
	  input.set("Data Layout", basis->functional);
	  RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
	}
	found = true;
      }

      if (plist.isType<std::string>("Type")) {
	std::string type = plist.get<std::string>("Type");
	if(type=="SIMPLE SOURCE") {
	  RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new Example::SimpleSource<EvalT,panzer::Traits>(key,*ir));
	  evaluators->push_back(e);

	  found = true;
	}
	else if(type=="TEMPERATURE_EXACT") {
	  RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new Example::SimpleSolution<EvalT,panzer::Traits>(key,*ir));
	  evaluators->push_back(e);

	  found = true;
	}
	else if(type=="L2 ERROR_CALC") {
	  {
	    std::vector<std::string> values(2);
	    values[0] = plist.get<std::string>("Field A");
	    values[1] = plist.get<std::string>("Field B");
  
	    std::vector<double> scalars(2); 
	    scalars[0] = 1.0; 
	    scalars[1] = -1.0;
  
	    Teuchos::ParameterList p;
	    p.set("Sum Name",key+"_DIFF"); // Name of sum
	    p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
	    p.set<RCP<const std::vector<double> > >("Scalars",Teuchos::rcpFromRef(scalars));
	    p.set("Data Layout",ir->dl_scalar);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
	    evaluators->push_back(e);
	  }

	  {
	    Teuchos::RCP<const panzer::PointRule> pr = ir;
	    std::vector<std::string> values(2);
	    values[0] = key+"_DIFF";
	    values[1] = key+"_DIFF";
  
	    Teuchos::ParameterList p;
	    p.set("Product Name",key);
	    p.set("Values Names",Teuchos::rcpFromRef(values));
	    p.set("Data Layout",ir->dl_scalar);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      rcp(new panzer::Product<EvalT,panzer::Traits>(p));
  
	    evaluators->push_back(e);
	  }

	  found = true;
	}
	else if(type=="H1 ERROR_CALC") {
	  // Compute L2 contribution
	  {
	    std::vector<std::string> values(2);
	    values[0] = plist.get<std::string>("Field A");
	    values[1] = plist.get<std::string>("Field B");
  
	    std::vector<double> scalars(2); 
	    scalars[0] = 1.0; 
	    scalars[1] = -1.0;
  
	    Teuchos::ParameterList p;
	    p.set("Sum Name",key+"_H1_L2DIFF"); // Name of sum
	    p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
	    p.set<RCP<const std::vector<double> > >("Scalars",Teuchos::rcpFromRef(scalars));
	    p.set("Data Layout",ir->dl_scalar);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
	    evaluators->push_back(e);
	  }

	  {
	    std::vector<std::string> values(2);
	    values[0] = "GRAD_"+plist.get<std::string>("Field A");
	    values[1] = "GRAD_"+plist.get<std::string>("Field B");
  
	    std::vector<double> scalars(2); 
	    scalars[0] = 1.0; 
	    scalars[1] = -1.0;
  
	    Teuchos::ParameterList p;
	    p.set("Sum Name","GRAD_"+key+"_DIFF"); // Name of sum
	    p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
	    p.set<RCP<const std::vector<double> > >("Scalars",Teuchos::rcpFromRef(scalars));
	    p.set("Data Layout",ir->dl_vector);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
	    evaluators->push_back(e);
	  }

	  {
	    Teuchos::RCP<const panzer::PointRule> pr = ir;
	    std::vector<std::string> values(2);
	    values[0] = "GRAD_"+key+"_DIFF";
	    values[1] = "GRAD_"+key+"_DIFF";
  
	    Teuchos::ParameterList p;
	    p.set("Product Name",key+"_H1Semi");
	    p.set("Values Names",Teuchos::rcpFromRef(values));
	    p.set("Data Layout",ir->dl_vector);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      panzer::buildEvaluator_DotProduct<EvalT,panzer::Traits>(key,*ir,values[0],values[1]);
  
	    evaluators->push_back(e);
	  }

	  {
	    Teuchos::RCP<const panzer::PointRule> pr = ir;
	    std::vector<std::string> values(2);
	    values[0] = key+"_H1_L2DIFF";
	    values[1] = key+"_H1_L2DIFF";
  
	    Teuchos::ParameterList p;
	    p.set("Product Name",key+"_L2");
	    p.set("Values Names",Teuchos::rcpFromRef(values));
	    p.set("Data Layout",ir->dl_scalar);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      rcp(new panzer::Product<EvalT,panzer::Traits>(p));
  
	    evaluators->push_back(e);
	  }

	  {
	    std::vector<std::string> values(2);
	    values[0] = key+"_L2";
	    values[1] = key+"_H1Semi";
  
	    Teuchos::ParameterList p;
	    p.set("Sum Name",key); // Name of sum
	    p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
	    p.set("Data Layout",ir->dl_scalar);
  
	    RCP< Evaluator<panzer::Traits> > e = 
	      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
	    evaluators->push_back(e);
	  }

	  found = true;
	}
      }

      if (!found) {
	std::stringstream msg;
	msg << "ClosureModelFactory failed to build evaluator for key \"" << key 
	    << "\"\nin model \"" << model_id 
	    << "\".  Please correct the type or add support to the \nfactory." <<std::endl;
	TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, msg.str());
      }

    }

    return evaluators;
  }

  class ClosureModelFactory_TemplateBuilder {
  public:
      
     template <typename EvalT>
     Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const 
     {
        return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>(new Example::ModelFactory<EvalT>) );
     }
      
  };
  
  /** The equation set serves two roles. The first is to let the panzer library
   * know which fields this equation set defines and their names. It registers
   * the evaluators required for a particular equation set. The level of the
   * granularity is largely up to a user. For instance this could be the momentum
   * or continuity equation in Navier-Stokes, or it could simply be the Navier-Stokes
   * equations. 
   *
   * Generally, this inherits from the panzer::EquationSet_DefaultImpl which takes
   * care of adding the gather (extract basis coefficients from solution vector) and 
   * scatter (using element matrices and vectors distribute and sum their values
   * to a global linear system) evaluators. These use data members that must be set by
   * the user.
   */
  template <typename EvalT>
  class PoissonEquationSet : public panzer::EquationSet_DefaultImpl<EvalT> {
  public:    

    /** In the constructor you set all the fields provided by this
     * equation set. 
     */
    PoissonEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		       const int& default_integration_order,
		       const panzer::CellData& cell_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       const bool build_transient_support);
    
    /** The specific evaluators are registered with the field manager argument.
     */
    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					       const panzer::FieldLibrary& field_library,
					       const Teuchos::ParameterList& user_data) const;

  };

  // ***********************************************************************
  template <typename EvalT>
  PoissonEquationSet<EvalT>::
  PoissonEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		     const int& default_integration_order,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) :
    panzer::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support )
  {
    // ********************
    // Validate and parse parameter list
    // ********************
    {    
      Teuchos::ParameterList valid_parameters;
      this->setDefaultValidParameters(valid_parameters);
    
      valid_parameters.set("Model ID","","Closure model id associated with this equaiton set");
      valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
      valid_parameters.set("Basis Order",1,"Order of the basis");
      valid_parameters.set("Integration Order",-1,"Order of the integration rule");
    
      params->validateParametersAndSetDefaults(valid_parameters);
    }
  
    std::string basis_type = params->get<std::string>("Basis Type");
    int basis_order = params->get<int>("Basis Order");
    int integration_order = params->get<int>("Integration Order");
    std::string model_id = params->get<std::string>("Model ID");

    // ********************
    // Panzer uses strings to match fields. In this section we define the
    // name of the fields provided by this equation set. This is a bit strange
    // in that this is not the fields necessarily required by this equation set.
    // For instance for the momentum equations in Navier-Stokes only the velocity
    // fields are added, the pressure field is added by continuity.
    //
    // In this case "TEMPERATURE" is the lone field.  We also name the gradient
    // for this field. These names automatically generate evaluators for "TEMPERATURE"
    // and "GRAD_TEMPERATURE" gathering the basis coefficients of "TEMPERATURE" and
    // the values of the TEMPERATURE and GRAD_TEMPERATURE fields at quadrature points.
    //
    // After all the equation set evaluators are added to a given field manager, the
    // panzer code adds in appropriate scatter evaluators to distribute the
    // entries into the residual and the Jacobian operator. These operators will be
    // "required" by the field manager and will serve as roots of evaluation tree.
    // The leaves of this tree will include the gather evaluators whose job it is to
    // gather the solution from a vector.
    // ********************

    // ********************
    // Assemble DOF names and Residual names
    // ********************

    this->addDOF("TEMPERATURE",basis_type,basis_order,integration_order);
    this->addDOFGrad("TEMPERATURE");
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative("TEMPERATURE");

    // ********************
    // Build Basis Functions and Integration Rules
    // ********************
   
    this->addClosureModel(model_id);

    this->setupDOFs();
  }

  // ***********************************************************************
  template <typename EvalT>
  void PoissonEquationSet<EvalT>::
  buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					const panzer::FieldLibrary& /* fl */,
					const Teuchos::ParameterList& /* user_data */) const
  {
    using panzer::BasisIRLayout;
    using panzer::EvaluatorStyle;
    using panzer::IntegrationRule;
    using panzer::Integrator_BasisTimesScalar;
    using panzer::Integrator_GradBasisDotVector;
    using panzer::Traits;
    using PHX::Evaluator;
    using std::string;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
  
    RCP<IntegrationRule> ir = this->getIntRuleForDOF("TEMPERATURE");
    RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF("TEMPERATURE"); 

    // ********************
    // Energy Equation
    // ********************

    // Transient Operator: Assembles \int \dot{T} v
    if (this->buildTransientSupport())
      {
	string resName("RESIDUAL_TEMPERATURE"), valName("DXDT_TEMPERATURE");
	double multiplier(1);
	RCP<Evaluator<Traits>> op = rcp(new
					Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
										   resName, valName, *basis, *ir, multiplier));
	this->template registerEvaluator<EvalT>(fm, op);
      }

    // Diffusion Operator: Assembles \int \nabla T \cdot \nabla v
    {
      double thermal_conductivity = 1.0;

      ParameterList p("Diffusion Residual");
      p.set("Residual Name", "RESIDUAL_TEMPERATURE");
      p.set("Flux Name", "GRAD_TEMPERATURE"); // this field is constructed by the panzer library
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", thermal_conductivity);
    
      RCP<Evaluator<Traits>> op = rcp(new
				      Integrator_GradBasisDotVector<EvalT, Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }
  
    // Source Operator
    {   
      string resName("RESIDUAL_TEMPERATURE"),
	valName("SOURCE_TEMPERATURE");
      double multiplier(-1);
      RCP<Evaluator<Traits>> op = rcp(new
				      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
										 resName, valName, *basis, *ir, multiplier));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

  // ***********************************************************************

  // A macro that defines a class to make construction of the equation sets easier
  //   - The equation set is constructed over a list of automatic differention types
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(PoissonEquationSet, PoissonEquationSet)

  // A user written factory that creates each equation set.  The key member here
  // is buildEquationSet
  class EquationSetFactory : public panzer::EquationSetFactory {
  public:

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		     const int& default_integration_order,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) const
    {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
	Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
         
      bool found = false; // this is used by PANZER_BUILD_EQSET_OBJECTS
         
      // macro checks if(ies.name=="Poisson") then an EquationSet_Energy object is constructed
      PANZER_BUILD_EQSET_OBJECTS("Poisson", PoissonEquationSet)
         
	// make sure your equation set has been found
	if(!found) {
	  std::string msg = "Error - the \"Equation Set\" called \"" + params->get<std::string>("Type") +
	    "\" is not a valid equation set identifier. Please supply the correct factory.\n";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
	}
         
      return eq_set;
    }
    
  };

}

#endif
