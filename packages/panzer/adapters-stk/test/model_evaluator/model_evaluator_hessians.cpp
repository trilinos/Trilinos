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

#include "Panzer_NodeType.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Phalanx_KokkosUtilities.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
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
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include <cstdio> // for get char
#include <fstream>

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
    RCP<panzer::UniqueGlobalIndexer<int,int> > param_dofManager;
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer;
    Teuchos::ParameterList user_data;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::EquationSetFactory> eqset_factory;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    Teuchos::ParameterList closure_models;
    Teuchos::RCP<EpetraVector_ReadOnly_GlobalEvaluationData> param_ged;
    std::vector<panzer::BC> bcs;
    Teuchos::RCP<panzer::BCStrategyFactory> bc_factory;
  };

  void buildAssemblyPieces(bool distr_parameter_on,
                           AssemblyPieces & ap,
                           const std::vector<std::string>& tangentParamNames = std::vector<std::string>());

  bool testEqualityOfVectorValues(const Thyra::VectorBase<double> & a, 
                                  const Thyra::VectorBase<double> & b, 
                                  double tolerance, bool write_to_cout=false);

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_hessians, d2g_dx2)
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::SpmdVectorBase<RealType> SpmdVectorType;
    typedef Thyra::LinearOpBase<RealType> OperatorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
    typedef panzer::ModelEvaluator<double> PME;


    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(distr_param_on,ap);

    int pIndex = -1;
    int rIndex = -1;

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    bool build_transient_support = true;
    RCP<PME> me 
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,p_values,Teuchos::null,ap.gd,build_transient_support,0.0));

    const double DENSITY_VALUE = 3.7;
    const double TEMPERATURE_VALUE = 2.0;
    const double PERTURB_VALUE = 0.1;

    // add in a flexible response
    {
      Teuchos::RCP<panzer::FunctionalResponse_Builder<int,int> > builder
        = Teuchos::rcp(new panzer::FunctionalResponse_Builder<int,int>);
 
      builder->comm = MPI_COMM_WORLD; // good enough
      builder->cubatureDegree = 1;
      builder->requiresCellIntegral = true;
      builder->quadPointField = "";

      std::vector<panzer::WorksetDescriptor> blocks;
      blocks.push_back(panzer::blockDescriptor("eblock-0_0"));
      blocks.push_back(panzer::blockDescriptor("eblock-1_0"));
      rIndex = me->addFlexibleResponse("INTEGRAND",blocks,builder); // integrate the density
    }

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),DENSITY_VALUE);
      pIndex = me->addDistributedParameter("DENSITY",th_param_lof->getThyraDomainSpace(),
                                           ap.param_ged,param_density,ap.param_dofManager);
    }

    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                   ap.closure_models,
                   ap.user_data,false,"");

    // test value
    {
      RCP<VectorType> g = Thyra::createMember(*me->get_g_space(rIndex)); 

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      OutArgs out_args = me->createOutArgs();
      out_args.set_g(0,g);

      me->evalModel(in_args, out_args);

      out << "RESPONSE = " << std::endl;
      out << Teuchos::describe(*g,Teuchos::VERB_EXTREME) << std::endl;

      Teuchos::ArrayRCP<const double> g_data;
      rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(g)->getLocalData(Teuchos::ptrFromRef(g_data));
      TEST_FLOATING_EQUALITY(g_data[0],TEMPERATURE_VALUE*TEMPERATURE_VALUE*DENSITY_VALUE*DENSITY_VALUE*DENSITY_VALUE,1e-12);
        // Volume of integral is 1.0, then just compute the integrand
    }

    RCP<VectorType> D2gDx2 = Thyra::createMember(*me->get_x_space()); 

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dx = Thyra::createMember(*me->get_x_space());
      Thyra::assign(dx.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      me->evalModel_D2gDx2(rIndex,in_args,dx,D2gDx2);

      out << "D2gDx2 = \n" << Teuchos::describe(*D2gDx2,Teuchos::VERB_EXTREME) << std::endl;
    }

    Teuchos::ArrayRCP<const double> D2gDx2_data;
    RCP<const SpmdVectorType> spmd_D2gDx2 = rcp_dynamic_cast<const SpmdVectorType>(D2gDx2);
    rcp_dynamic_cast<const SpmdVectorType>(D2gDx2)->getLocalData(Teuchos::ptrFromRef(D2gDx2_data));
    double scale = 2.0*DENSITY_VALUE*DENSITY_VALUE*DENSITY_VALUE*PERTURB_VALUE;
    double int_phi = 1.0/192.0;
    for(int i=0;i<D2gDx2_data.size();i++) {
      out << D2gDx2_data[i]  << " " << D2gDx2_data[i]/(scale*int_phi) << std::endl;
      bool a = std::fabs(D2gDx2_data[i]-scale*int_phi)/(scale*int_phi)          <= 1e-14;
      bool b = std::fabs(D2gDx2_data[i]-2.0*scale*int_phi)/(2.0*scale*int_phi)  <= 1e-14;
      bool c = std::fabs(D2gDx2_data[i]-4.0*scale*int_phi)/(4.0*scale*int_phi)  <= 1e-14;
      bool d = (D2gDx2_data[i]==0.0);

      TEST_ASSERT(a || b || c || d);
    }
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_hessians, d2g_dp2)
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::SpmdVectorBase<RealType> SpmdVectorType;
    typedef Thyra::LinearOpBase<RealType> OperatorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
    typedef panzer::ModelEvaluator<double> PME;


    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(distr_param_on,ap);

    int pIndex = -1;
    int rIndex = -1;

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    bool build_transient_support = true;
    RCP<PME> me 
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,p_values,Teuchos::null,ap.gd,build_transient_support,0.0));

    const double DENSITY_VALUE = 3.7;
    const double TEMPERATURE_VALUE = 2.0;
    const double PERTURB_VALUE = 0.1;

    // add in a flexible response
    {
      Teuchos::RCP<panzer::FunctionalResponse_Builder<int,int> > builder
        = Teuchos::rcp(new panzer::FunctionalResponse_Builder<int,int>);
 
      builder->comm = MPI_COMM_WORLD; // good enough
      builder->cubatureDegree = 1;
      builder->requiresCellIntegral = true;
      builder->quadPointField = "";

      std::vector<panzer::WorksetDescriptor> blocks;
      blocks.push_back(panzer::blockDescriptor("eblock-0_0"));
      blocks.push_back(panzer::blockDescriptor("eblock-1_0"));
      rIndex = me->addFlexibleResponse("INTEGRAND",blocks,builder); // integrate the density
    }

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),DENSITY_VALUE);
      pIndex = me->addDistributedParameter("DENSITY",th_param_lof->getThyraDomainSpace(),
                                           ap.param_ged,param_density,ap.param_dofManager);
    }

    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                   ap.closure_models,
                   ap.user_data,false,"");

    // test value
    {
      RCP<VectorType> g = Thyra::createMember(*me->get_g_space(rIndex)); 

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      OutArgs out_args = me->createOutArgs();
      out_args.set_g(0,g);

      me->evalModel(in_args, out_args);

      out << "RESPONSE = " << std::endl;
      out << Teuchos::describe(*g,Teuchos::VERB_EXTREME) << std::endl;

      Teuchos::ArrayRCP<const double> g_data;
      rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(g)->getLocalData(Teuchos::ptrFromRef(g_data));
      TEST_FLOATING_EQUALITY(g_data[0],TEMPERATURE_VALUE*TEMPERATURE_VALUE*DENSITY_VALUE*DENSITY_VALUE*DENSITY_VALUE,1e-12);
        // Volume of integral is 1.0, then just compute the integrand
    }

    RCP<VectorType> D2gDp2 = Thyra::createMember(*me->get_p_space(pIndex)); 

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dx = Thyra::createMember(*me->get_x_space());
      Thyra::assign(dx.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      me->evalModel_D2gDp2(rIndex,pIndex,in_args,dx,D2gDp2);

      out << "D2gDp2 = \n" << Teuchos::describe(*D2gDp2,Teuchos::VERB_EXTREME) << std::endl;
    }

    Teuchos::ArrayRCP<const double> D2gDp2_data;
    RCP<const SpmdVectorType> spmd_D2gDp2 = rcp_dynamic_cast<const SpmdVectorType>(D2gDp2);
    rcp_dynamic_cast<const SpmdVectorType>(D2gDp2)->getLocalData(Teuchos::ptrFromRef(D2gDp2_data));
    double scale = 6.0*TEMPERATURE_VALUE*TEMPERATURE_VALUE*DENSITY_VALUE*PERTURB_VALUE;
    double int_phi = 1.0/192.0;
    for(int i=0;i<D2gDp2_data.size();i++) {
      out << D2gDp2_data[i]  << " " << D2gDp2_data[i]/(scale*int_phi) << std::endl;
      bool a = std::fabs(D2gDp2_data[i]-scale*int_phi)/(scale*int_phi)          <= 1e-14;
      bool b = std::fabs(D2gDp2_data[i]-2.0*scale*int_phi)/(2.0*scale*int_phi)  <= 1e-14;
      bool c = std::fabs(D2gDp2_data[i]-4.0*scale*int_phi)/(4.0*scale*int_phi)  <= 1e-14;
      bool d = (D2gDp2_data[i]==0.0);

      TEST_ASSERT(a || b || c || d);
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
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
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
  
  void buildAssemblyPieces(bool distr_parameter_on,
                           AssemblyPieces & ap,
                           const std::vector<std::string>& tangentParamNames)
  {
    using Teuchos::RCP;
  
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);
    
    panzer_stk_classic::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk_classic::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm 
       = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Comm);

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    ap.fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    ap.eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    ap.bc_factory = Teuchos::rcp(new user_app::BCFactory);
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
                                 ap.physicsBlocks,
                                 tangentParamNames);
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
        = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpiComm,dofManager));
      ap.lof = linObjFactory;
    }

    // build the dof manager and LOF for DENSITY control
    if(distr_parameter_on) {
      Teuchos::RCP<panzer::DOFManager<int,int> > dofManager 
          = Teuchos::rcp(new panzer::DOFManager<int,int>(conn_manager,MPI_COMM_WORLD));

      Teuchos::RCP<Intrepid2FieldPattern> fp 
          = Teuchos::rcp(new Intrepid2FieldPattern(panzer::createIntrepid2Basis<double,Kokkos::DynRankView<double,PHX::Device> >("HGrad",1,mesh->getCellTopology("eblock-0_0"))));
      dofManager->addField("eblock-0_0","DENSITY",fp);
      dofManager->addField("eblock-1_0","DENSITY",fp);

      dofManager->setOrientationsRequired(false);
      dofManager->buildGlobalUnknowns();

      // build a nonsquare LOF for the parameter vector
      Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > linObjFactory
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpiComm,ap.dofManager,dofManager));

      Teuchos::RCP<Epetra_Map> uniqueMap = linObjFactory->getColMap(0);
      Teuchos::RCP<Epetra_Map> ghostedMap = linObjFactory->getGhostedColMap(0);
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*ghostedMap,*uniqueMap));

      ap.param_dofManager = dofManager;
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

    closure_models.sublist("solid").sublist("INTEGRAND").set<std::string>("Type","Product");
    closure_models.sublist("solid").sublist("INTEGRAND").set<std::string>("Term Names","TEMPERATURE,TEMPERATURE,DENSITY,DENSITY,DENSITY");

    ap.closure_models = closure_models;

    ap.user_data = Teuchos::ParameterList("User Data");

    ap.fmb->setWorksetContainer(wkstContainer);
    ap.fmb->setupVolumeFieldManagers(ap.physicsBlocks,ap.cm_factory,closure_models,*ap.lof,ap.user_data);
    ap.fmb->setupBCFieldManagers(bcs,ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,*ap.bc_factory,
                                 closure_models,*ap.lof,ap.user_data);
  }


}
