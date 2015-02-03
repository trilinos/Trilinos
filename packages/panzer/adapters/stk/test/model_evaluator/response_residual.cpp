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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_Response_Residual.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Thyra_LinearOpTester.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  struct AssemblyPieces {
    RCP<panzer::FieldManagerBuilder> fmb;  
    RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    RCP<panzer::GlobalData> gd;
    RCP<panzer::LinearObjFactory<panzer::Traits> > lof;
    RCP<panzer::LinearObjFactory<panzer::Traits> > param_lof;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager;
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer;
    Teuchos::ParameterList user_data;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::EquationSetFactory> eqset_factory;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    Teuchos::ParameterList closure_models;
    std::vector<panzer::BC> bcs;
    Teuchos::RCP<EpetraVector_ReadOnly_GlobalEvaluationData> param_ged;
  };

  void buildAssemblyPieces(bool parameter_on,bool distr_parameter_on,
                           AssemblyPieces & ap);

  bool testEqualityOfVectorValues(const Thyra::VectorBase<double> & a, 
                                  const Thyra::VectorBase<double> & b, 
                                  double tolerance, bool write_to_cout=false);

  // Test that the response library can build the correct residual and jacobian
  TEUCHOS_UNIT_TEST(response_residual, basic)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::LinearOpBase<RealType> OperatorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    bool parameter_on = true;
    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(parameter_on,distr_param_on,ap);

    // Use the model evaluator to setup a bunch of comparisons
    ///////////////////////////////////////////////////////////////////////
    RCP<VectorType> x;
    RCP<VectorType> x_dot;
    RCP<VectorType> f_me;
    RCP<OperatorType> J_me;
    double alpha = 1.3, beta = 0.2;

    RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);
    RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraRangeSpace());
    Thyra::assign(param_density.ptr(),3.7);
    {
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef panzer::ModelEvaluator<double> PME;

      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = true;

      RCP<PME> me = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,Teuchos::null,ap.gd,build_transient_support,0.0));
      me->addDistributedParameter("DENSITY",th_param_lof->getThyraRangeSpace(),ap.param_ged,param_density);

      x = Thyra::createMember(*me->get_x_space());
      x_dot = Thyra::createMember(*me->get_x_space());

      Thyra::randomize(-1.0,1.0,x.ptr());
      Thyra::randomize(-1.0,1.0,x_dot.ptr());

      InArgs in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_x_dot(x_dot);
      in_args.set_alpha(alpha);
      in_args.set_beta(beta);
      
      OutArgs out_args = me->createOutArgs();
      f_me = Thyra::createMember(*me->get_f_space());
      J_me = me->create_W_op();
      out_args.set_f(f_me);
      out_args.set_W_op(J_me);

      me->evalModel(in_args, out_args);
    }

    bool residualType = true;
    RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>(ap.wkstContainer,ap.dofManager,ap.lof,residualType)); 
 
    // verify the residual type
    TEST_ASSERT(rLibrary->isResidualType());

    // test that responses are the right type
    {
      RCP<ResponseBase> response_residual = rLibrary->getResponse<Traits::Residual>("RESIDUAL");
      RCP<ResponseBase> response_jacobian = rLibrary->getResponse<Traits::Jacobian>("RESIDUAL");

      TEST_ASSERT(response_residual!=Teuchos::null);
      TEST_ASSERT(response_jacobian!=Teuchos::null);

      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Residual<Traits::Residual> >(response_residual,true));
      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(response_jacobian,true));
    }

    // build residual response evaluators
    {
      user_app::BCFactory bc_factory;

      TEST_ASSERT(!rLibrary->responseEvaluatorsBuilt());

      rLibrary->buildResidualResponseEvaluators(ap.physicsBlocks,*ap.eqset_factory,ap.bcs,bc_factory,
                                                ap.cm_factory,ap.closure_models,ap.user_data);

      TEST_ASSERT(rLibrary->responseEvaluatorsBuilt());
    }

    RCP<Response_Residual<Traits::Residual> > response_residual = 
      rcp_dynamic_cast<Response_Residual<Traits::Residual> >(rLibrary->getResponse<Traits::Residual>("RESIDUAL"));
    RCP<Response_Residual<Traits::Jacobian> > response_jacobian = 
      rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(rLibrary->getResponse<Traits::Jacobian>("RESIDUAL"));

    Teuchos::RCP<EpetraVector_ReadOnly_GlobalEvaluationData> resp_param_ged =
      Teuchos::rcp(new EpetraVector_ReadOnly_GlobalEvaluationData(*ap.param_ged));
    resp_param_ged->setUniqueVector(param_density);

    // evaluate residual responses
    {
      // set solution vectors
      RCP<LinearObjContainer> loc = ap.lof->buildLinearObjContainer();
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_x_th(x);
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_dxdt_th(x_dot);

      // allocate fill vectors
      RCP<LinearObjContainer> gloc = ap.lof->buildGhostedLinearObjContainer();
      ap.lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::DxDt,*gloc);

      // setup output arguments for the residual response 
      response_residual->setResidual(response_residual->allocateResidualVector());

      // setup in arguments
      AssemblyEngineInArgs ae_inargs(gloc,loc);
      ae_inargs.addGlobalEvaluationData("DENSITY",resp_param_ged);
      ae_inargs.alpha = alpha;
      ae_inargs.beta = beta;
      ae_inargs.evaluate_transient_terms = true;
      rLibrary->addResponsesToInArgs<Traits::Residual>(ae_inargs);

      // evaluate
      rLibrary->evaluate<Traits::Residual>(ae_inargs);

      TEST_ASSERT(testEqualityOfVectorValues(*response_residual->getResidual(),*f_me,1e-16));
    }

    // evaluate jacobian responses
    {
      // set solution vectors
      RCP<LinearObjContainer> loc = ap.lof->buildLinearObjContainer();
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_x_th(x);
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_dxdt_th(x_dot);

      // allocate fill vectors
      RCP<LinearObjContainer> gloc = ap.lof->buildGhostedLinearObjContainer();
      ap.lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::DxDt,*gloc);

      // setup output arguments for the residual response 
      response_jacobian->setJacobian(response_jacobian->allocateJacobian());

      // setup in arguments
      AssemblyEngineInArgs ae_inargs(gloc,loc);
      ae_inargs.addGlobalEvaluationData("DENSITY",resp_param_ged);
      ae_inargs.alpha = alpha;
      ae_inargs.beta = beta;
      ae_inargs.evaluate_transient_terms = true;
      rLibrary->addResponsesToInArgs<Traits::Jacobian>(ae_inargs);

      // evaluate
      rLibrary->evaluate<Traits::Jacobian>(ae_inargs);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
      tester.set_all_error_tol(1e-16);      
      tester.num_random_vectors(20);      

      Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
      const bool op_cmp = tester.compare( *J_me, *response_jacobian->getJacobian(), Teuchos::ptrFromRef(fout));
      TEST_ASSERT(op_cmp);
    }
  }

  bool comparison(double a,double b)
  { 
    double exact = std::fabs(a+b)/2.0;

    if(exact>1e-12)
      return (std::fabs(a-b)/exact) < 1e-12;
    else 
      return (std::fabs(a-b)) < 1e-12;
  }

  // Test that the response library can build the correct residual and jacobian
  TEUCHOS_UNIT_TEST(response_residual, distr_parameter)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::LinearOpBase<RealType> OperatorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    bool parameter_on = true;
    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(parameter_on,distr_param_on,ap);

    RCP<ThyraObjFactory<double> > th_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.lof);
    RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

    // Use the model evaluator to setup a bunch of comparisons
    ///////////////////////////////////////////////////////////////////////
    RCP<VectorType> x = Thyra::createMember(th_lof->getThyraDomainSpace());
    RCP<VectorType> x_dot = Thyra::createMember(th_lof->getThyraDomainSpace());
    RCP<VectorType> f = Thyra::createMember(th_lof->getThyraRangeSpace());

    RCP<OperatorType> J = th_param_lof->getThyraMatrix();

    Thyra::randomize(-1.0,1.0,x.ptr());
    Thyra::put_scalar(1.0,x_dot.ptr());

    double alpha = Teuchos::ScalarTraits<double>::nan();
    double beta = Teuchos::ScalarTraits<double>::nan();

    bool residualType = true;
    RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>(ap.wkstContainer,ap.dofManager,ap.param_lof,residualType)); 
 
    // verify the residual type
    TEST_ASSERT(rLibrary->isResidualType());

    // test that responses are the right type
    {
      RCP<ResponseBase> response_jacobian = rLibrary->getResponse<Traits::Jacobian>("RESIDUAL");

      TEST_ASSERT(response_jacobian!=Teuchos::null);

      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(response_jacobian,true));
    }

    // build residual response evaluators
    {
      user_app::BCFactory bc_factory;

      TEST_ASSERT(!rLibrary->responseEvaluatorsBuilt());

      rLibrary->buildResidualResponseEvaluators(ap.physicsBlocks,*ap.eqset_factory,ap.bcs,bc_factory,
                                                ap.cm_factory,ap.closure_models,ap.user_data);

      TEST_ASSERT(rLibrary->responseEvaluatorsBuilt());
    }

    RCP<Response_Residual<Traits::Jacobian> > response_jacobian = 
      rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(rLibrary->getResponse<Traits::Jacobian>("RESIDUAL"));

    Teuchos::RCP<EpetraVector_ReadOnly_GlobalEvaluationData> resp_param_ged =
      Teuchos::rcp(new EpetraVector_ReadOnly_GlobalEvaluationData(*ap.param_ged));

    RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraRangeSpace());
    Thyra::assign(param_density.ptr(),3.7);
    resp_param_ged->setUniqueVector(param_density);

    // evaluate residual responses
    {
      // set solution vectors
      RCP<LinearObjContainer> loc = ap.lof->buildLinearObjContainer();
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_x_th(x);
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_dxdt_th(x_dot);

      // allocate fill vectors
      RCP<LinearObjContainer> gloc = ap.param_lof->buildGhostedLinearObjContainer();
      ap.lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::DxDt,*gloc);

      // setup output arguments for the residual response 
      response_jacobian->setJacobian(response_jacobian->allocateJacobian());

      // setup in arguments
      AssemblyEngineInArgs ae_inargs(gloc,loc);
      ae_inargs.addGlobalEvaluationData("DENSITY",resp_param_ged);
      ae_inargs.alpha = alpha;
      ae_inargs.beta = beta;
      ae_inargs.evaluate_transient_terms = true;
      ae_inargs.sensitivities_name = "DENSITY";
      ae_inargs.gather_seeds.push_back(1.0); // gather seed index 0 (see closure model)
      rLibrary->addResponsesToInArgs<Traits::Jacobian>(ae_inargs);

      // evaluate
      rLibrary->evaluate<Traits::Jacobian>(ae_inargs);

      // std::cout << Teuchos::describe(*response_jacobian->getJacobian(),Teuchos::VERB_EXTREME) << std::endl;
    }

    // test the jacobian for correctness
    {
 
      // build vectors for each type of node
      std::vector<double> corner(4);
      corner[0] = 1./1152.; corner[1] = 1./576.;
      corner[2] = 1./288.;  corner[3] = 1./128.;

      std::vector<double> edge(6);
      edge[0] = 1./1152.; edge[1] = 1./1152.; edge[2] = 1./576.;
      edge[3] = 1./576.;  edge[4] = 1./288.;  edge[5] = 1./144.;

      std::vector<double> volume(9); 
      volume[0] = 1./1152.; volume[1] = 1./1152.; volume[2] = 1./1152.;
      volume[3] = 1./1152.; volume[4] = 1./288.;  volume[5] = 1./288.;
      volume[6] = 1./288.;  volume[7] = 1./288.;  volume[8] = 1./72.;
      
      RCP<const Thyra::LinearOpBase<double> > thJac = response_jacobian->getJacobian();
      RCP<const Epetra_CrsMatrix> jac = rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*thJac));

      TEUCHOS_ASSERT(jac!=Teuchos::null);

      for(int i=0;i<jac->NumMyRows();i++) {
        int numEntries = -1;
        int * indices = 0;
        double * values = 0;

        // get a view of the row entries
        jac->ExtractMyRowView(i,numEntries,values,indices);

        TEUCHOS_ASSERT(numEntries>0);

        // sort the row entries
        std::vector<double> sorted_values(numEntries);
        for(int j=0;j<numEntries;j++)
          sorted_values.push_back(values[j]);
        std::sort(sorted_values.begin(),sorted_values.end());

        if(sorted_values[0]!=sorted_values[sorted_values.size()-1]) {
          std::vector<double>::const_iterator found_itr;

          ///////////////////////////////////////////////////////////////////////

          found_itr = std::find_end(sorted_values.begin(),sorted_values.end(),
                                    corner.begin(),corner.end(),comparison);

          // test passed this row corresponds to a corner
          if(found_itr!=corner.end()) 
            continue;

          ///////////////////////////////////////////////////////////////////////

          found_itr = std::find_end(sorted_values.begin(),sorted_values.end(),
                                    edge.begin(),edge.end(),comparison);

          // test passed this row corresponds to a edge
          if(found_itr!=edge.end()) 
            continue;

          ///////////////////////////////////////////////////////////////////////

          found_itr = std::find_end(sorted_values.begin(),sorted_values.end(),
                                    volume.begin(),volume.end(),comparison);

          // test passed this row corresponds to a volume
          if(found_itr!=volume.end()) 
            continue;

          ///////////////////////////////////////////////////////////////////////

          TEST_ASSERT(false); // non of the required row types were found

          out << "Row didn't match expectation " << i << ": ";
          for(std::size_t i=0;i<sorted_values.size();i++)
            out << sorted_values[i] << " ";
          out << std::endl;
        }
        else {
          TEST_ASSERT(sorted_values[0]==0.0); 
          TEST_ASSERT(sorted_values[sorted_values.size()-1]==0.0);
        }
                      
      }
    }
  }

  bool testEqualityOfVectorValues(const Thyra::VectorBase<double> & a, 
                                  const Thyra::VectorBase<double> & b, 
                                  double tolerance, bool write_to_cout)
  {  
    bool is_equal = true;

    TEUCHOS_ASSERT(a.space()->dim() == b.space()->dim());

    Teuchos::ArrayRCP<const double> a_data,b_data;
    dynamic_cast<const Thyra::SpmdVectorBase<double> &>(a).getLocalData(Teuchos::ptrFromRef(a_data));
    dynamic_cast<const Thyra::SpmdVectorBase<double> &>(b).getLocalData(Teuchos::ptrFromRef(b_data));

    for (int i = 0; i < a_data.size(); ++i) {
      
      std::string output = "    equal!: ";
      
      if (std::fabs(a_data[i] - b_data[i]) > tolerance) {
	is_equal = false;
	output = "NOT equal!: ";
      }
      
      if (write_to_cout)
	std::cout << output << a_data[i] << " - " << b_data[i] << " = " << (a_data[i] - b_data[i]) << std::endl;
    }

    int globalSuccess = -1;
    Teuchos::RCP<const Teuchos::Comm<Teuchos::Ordinal> > comm = Teuchos::DefaultComm<Teuchos::Ordinal>::getComm();
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, is_equal ? 0 : 1, Teuchos::outArg(globalSuccess) );
    return (globalSuccess==0);
  }

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",2);
    }

    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",2);
    }

    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }    
    {
      std::size_t bc_id = 1;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }   
    {
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }
  }
  
  void buildAssemblyPieces(bool parameter_on,bool distr_parameter_on,
                           AssemblyPieces & ap)
  {
    using Teuchos::RCP;
  
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",4);
    pl->set("Y Elements",4);
    
    panzer_stk_classic::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk_classic::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Comm);

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> & bcs = ap.bcs;
    testInitialzation(ipb, bcs);
    

    ap.fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    ap.eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    ap.gd = panzer::createGlobalData();
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
      
      int default_integration_order = 1;
      
      bool build_transient_support = true;
      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
                                 ap.eqset_factory,
				 ap.gd,
			         build_transient_support,
                                 ap.physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    // build WorksetContainer
    Teuchos::RCP<panzer_stk_classic::WorksetFactory> wkstFactory 
       = Teuchos::rcp(new panzer_stk_classic::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,ap.physicsBlocks,workset_size));
    ap.wkstContainer = wkstContainer;

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk_classic::STKConnManager<int>(mesh));

    // build the state dof manager and LOF
    {
      panzer::DOFManagerFactory<int,int> globalIndexerFactory;
      RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
           = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),ap.physicsBlocks,conn_manager);
      ap.dofManager = dofManager;

      Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
          // = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64>(Comm.getConst(),dofManager));
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpiComm,dofManager));
      ap.lof = linObjFactory;
    }

    // build the dof manager and LOF for DENSITY control
    if(distr_parameter_on) {
      Teuchos::RCP<panzer::DOFManager<int,int> > dofManager 
          = Teuchos::rcp(new panzer::DOFManager<int,int>(conn_manager,MPI_COMM_WORLD));

      Teuchos::RCP<IntrepidFieldPattern> fp 
          = Teuchos::rcp(new IntrepidFieldPattern(panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >("HGrad",1,mesh->getCellTopology("eblock-0_0"))));
      dofManager->addField("eblock-0_0","DENSITY",fp);
      dofManager->addField("eblock-1_0","DENSITY",fp);

      dofManager->setOrientationsRequired(false);
      dofManager->buildGlobalUnknowns();

      // build a nonsquare LOF for the parameter vector
      Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > linObjFactory
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpiComm,ap.dofManager,dofManager));

      Teuchos::RCP<Epetra_Map> uniqueMap = linObjFactory->getColMap();
      Teuchos::RCP<Epetra_Map> ghostedMap = linObjFactory->getGhostedColMap();
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*ghostedMap,*uniqueMap));

      ap.param_ged = Teuchos::rcp(new EpetraVector_ReadOnly_GlobalEvaluationData(importer,ghostedMap,uniqueMap));
      ap.param_lof = linObjFactory;
    }

    ap.rLibrary = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,ap.dofManager,ap.lof)); 

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    // Add in the application specific closure model factory
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    cm_builder.setDistributedParameterLOF(ap.param_lof);
    ap.cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    if(parameter_on)
       closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<std::string>("Type","Parameter");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    if(distr_parameter_on)
      closure_models.sublist("solid").sublist("DENSITY").set("Type","Distributed Parameter");
    else
      closure_models.sublist("solid").sublist("DENSITY").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("HEAT_CAPACITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_HEAT_CAPACITY").set<double>("Value",1.0);
    ap.closure_models = closure_models;

    ap.user_data = Teuchos::ParameterList("User Data");

    ap.fmb->setWorksetContainer(wkstContainer);
    ap.fmb->setupVolumeFieldManagers(ap.physicsBlocks,ap.cm_factory,closure_models,*ap.lof,ap.user_data);
    ap.fmb->setupBCFieldManagers(bcs,ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,bc_factory,closure_models,*ap.lof,ap.user_data);
  }


}
