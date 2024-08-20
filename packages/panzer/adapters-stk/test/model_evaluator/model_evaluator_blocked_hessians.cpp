// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_NodeType.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

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
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
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
#include "Panzer_GlobalIndexer_Utilities.hpp"

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
    Teuchos::RCP<panzer::FieldManagerBuilder> fmb;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    Teuchos::RCP<panzer::GlobalData> gd;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > param_lof;
    Teuchos::RCP<panzer::GlobalIndexer> dofManager;
    Teuchos::RCP<panzer::GlobalIndexer> param_dofManager;
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
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2g_dx2)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::SpmdVectorBase<RealType> SpmdVectorType;
    typedef Thyra::ProductVectorBase<RealType> ProductVectorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
    typedef panzer::ModelEvaluator<double> PME;


    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(distr_param_on,ap);

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
      me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
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

    for(int b=0;b<2;b++) {
      out << "checking vector: " << b << std::endl;
      Teuchos::ArrayRCP<const double> D2gDx2_data;

      RCP<const ProductVectorType> product_D2gDx2 = rcp_dynamic_cast<const ProductVectorType>(D2gDx2,true);
      RCP<const SpmdVectorType> spmd_D2gDx2 = rcp_dynamic_cast<const SpmdVectorType>(product_D2gDx2->getVectorBlock(b),true);
      spmd_D2gDx2->getLocalData(Teuchos::ptrFromRef(D2gDx2_data));
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
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2g_dp2)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::SpmdVectorBase<RealType> SpmdVectorType;

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
      pIndex = me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
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

      RCP<VectorType> dp = Thyra::createMember(*me->get_p_space(pIndex));
      Thyra::assign(dp.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      me->evalModel_D2gDp2(rIndex,pIndex,in_args,dp,D2gDp2);

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

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2g_dpdx)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::SpmdVectorBase<RealType> SpmdVectorType;

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
      pIndex = me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
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

    RCP<VectorType> D2gDpDx = Thyra::createMember(*me->get_p_space(pIndex));

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dx = Thyra::createMember(*me->get_x_space());
      Thyra::assign(dx.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      me->evalModel_D2gDpDx(rIndex,pIndex,in_args,dx,D2gDpDx);

      out << "D2gDpDx = \n" << Teuchos::describe(*D2gDpDx,Teuchos::VERB_EXTREME) << std::endl;
    }

    Teuchos::ArrayRCP<const double> D2gDpDx_data;

    RCP<const SpmdVectorType> spmd_D2gDpDx = rcp_dynamic_cast<const SpmdVectorType>(D2gDpDx);
    spmd_D2gDpDx->getLocalData(Teuchos::ptrFromRef(D2gDpDx_data));
    double scale = 6.0*TEMPERATURE_VALUE*DENSITY_VALUE*DENSITY_VALUE*PERTURB_VALUE;
    double int_phi = 1.0/192.0;
    for(int i=0;i<D2gDpDx_data.size();i++) {
      out << D2gDpDx_data[i]  << " " << D2gDpDx_data[i]/(scale*int_phi) << std::endl;
      bool a = std::fabs(D2gDpDx_data[i]-scale*int_phi)/(scale*int_phi)          <= 1e-14;
      bool b = std::fabs(D2gDpDx_data[i]-2.0*scale*int_phi)/(2.0*scale*int_phi)  <= 1e-14;
      bool c = std::fabs(D2gDpDx_data[i]-4.0*scale*int_phi)/(4.0*scale*int_phi)  <= 1e-14;
      bool d = (D2gDpDx_data[i]==0.0);

      TEST_ASSERT(a || b || c || d);
    }
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2g_dxdp)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::SpmdVectorBase<RealType> SpmdVectorType;
    typedef Thyra::ProductVectorBase<RealType> ProductVectorType;

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
      pIndex = me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
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

    RCP<VectorType> D2gDxDp = Thyra::createMember(*me->get_x_space());

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dp = Thyra::createMember(*me->get_p_space(pIndex));
      Thyra::assign(dp.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);

      me->evalModel_D2gDxDp(rIndex,pIndex,in_args,dp,D2gDxDp);

      out << "D2gDxDp = \n" << Teuchos::describe(*D2gDxDp,Teuchos::VERB_EXTREME) << std::endl;
    }

    for(int b=0;b<2;b++) {
      out << "checking vector: " << b << std::endl;
      Teuchos::ArrayRCP<const double> D2gDxDp_data;

      RCP<const ProductVectorType> product_D2gDxDp = rcp_dynamic_cast<const ProductVectorType>(D2gDxDp,true);
      RCP<const SpmdVectorType> spmd_D2gDxDp = rcp_dynamic_cast<const SpmdVectorType>(product_D2gDxDp->getVectorBlock(b));
      spmd_D2gDxDp->getLocalData(Teuchos::ptrFromRef(D2gDxDp_data));
      double scale = 6.0*TEMPERATURE_VALUE*DENSITY_VALUE*DENSITY_VALUE*PERTURB_VALUE;
      double int_phi = 1.0/192.0;
      for(int i=0;i<D2gDxDp_data.size();i++) {
        out << D2gDxDp_data[i]  << " " << D2gDxDp_data[i]/(scale*int_phi) << std::endl;
        bool a = std::fabs(D2gDxDp_data[i]-scale*int_phi)/(scale*int_phi)          <= 1e-14;
        bool b = std::fabs(D2gDxDp_data[i]-2.0*scale*int_phi)/(2.0*scale*int_phi)  <= 1e-14;
        bool c = std::fabs(D2gDxDp_data[i]-4.0*scale*int_phi)/(4.0*scale*int_phi)  <= 1e-14;
        bool d = (D2gDxDp_data[i]==0.0);

        TEST_ASSERT(a || b || c || d);
      }
    }
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2f_dx2)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::LinearOpBase<RealType> OperatorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
    typedef panzer::ModelEvaluator<double> PME;


    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(distr_param_on,ap);

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    bool build_transient_support = true;
    RCP<PME> me
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,p_values,Teuchos::null,ap.gd,build_transient_support,0.0));

    const double DENSITY_VALUE = 3.7;
    const double TEMPERATURE_VALUE = 2.0;
    const double PERTURB_VALUE = 0.1;

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),DENSITY_VALUE);
      me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
                                  ap.param_ged,param_density,ap.param_dofManager);
    }

    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                   ap.closure_models,
                   ap.user_data,false,"");

    // panzer::printMeshTopology(out,*ap.dofManager);


    RCP<OperatorType> D2fDx2 = me->create_W_op();

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dx = Thyra::createMember(*me->get_x_space());
      Thyra::assign(dx.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0/0.1);
      in_args.set_beta(1.0);

      me->evalModel_D2fDx2(in_args,dx,D2fDx2);

      out << "D2fDx2 = \n" << Teuchos::describe(*D2fDx2,Teuchos::VERB_EXTREME) << std::endl;
    }

    RCP<OperatorType> W = me->create_W_op();
    {
      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      InArgs in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0);
      in_args.set_beta(0.0);

      OutArgs out_args = me->createOutArgs();
      out_args.set_W_op(W);

      me->evalModel(in_args,out_args);

      out << "W = \n" << Teuchos::describe(*W,Teuchos::VERB_EXTREME) << std::endl;
    }

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-15);
    tester.num_random_vectors(20);

    double scaling_of_mass_matrix = -2.0*PERTURB_VALUE*DENSITY_VALUE*DENSITY_VALUE;

    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(out));
    const bool op_cmp = tester.compare( *Thyra::scaleAndAdjoint(scaling_of_mass_matrix,Thyra::NOTRANS,W.getConst()), *D2fDx2, Teuchos::ptrFromRef(fout));
    TEST_ASSERT(op_cmp);
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2f_dxdp)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
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

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    bool build_transient_support = true;
    RCP<PME> me
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,p_values,Teuchos::null,ap.gd,build_transient_support,0.0));

    const double DENSITY_VALUE = 3.7;
    const double TEMPERATURE_VALUE = 2.0;
    const double PERTURB_VALUE = 0.1;

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),DENSITY_VALUE);
      pIndex = me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
                                           ap.param_ged,param_density,ap.param_dofManager);
    }

    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                   ap.closure_models,
                   ap.user_data,false,"");

    // panzer::printMeshTopology(out,*ap.dofManager);


    RCP<OperatorType> D2fDxDp = me->create_W_op();

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dp = Thyra::createMember(*me->get_p_space(0));
      Thyra::assign(dp.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0/0.1);
      in_args.set_beta(1.0);

      me->evalModel_D2fDxDp(pIndex,in_args,dp,D2fDxDp);

      out << "D2fDxDp = \n" << Teuchos::describe(*D2fDxDp,Teuchos::VERB_EXTREME) << std::endl;
    }

    RCP<OperatorType> W = me->create_W_op();
    {
      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      InArgs in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0);
      in_args.set_beta(0.0);

      OutArgs out_args = me->createOutArgs();
      out_args.set_W_op(W);

      me->evalModel(in_args,out_args);

      out << "W = \n" << Teuchos::describe(*W,Teuchos::VERB_EXTREME) << std::endl;
    }

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-15);
    tester.num_random_vectors(20);

    double scaling_of_mass_matrix = -2.0*TEMPERATURE_VALUE*2.0*DENSITY_VALUE*PERTURB_VALUE;

    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
    const bool op_cmp = tester.compare( *Thyra::scaleAndAdjoint(scaling_of_mass_matrix,Thyra::NOTRANS,W.getConst()), *D2fDxDp, Teuchos::ptrFromRef(fout));
    TEST_ASSERT(op_cmp);
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2f_dpdx)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
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

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    bool build_transient_support = true;
    RCP<PME> me
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,p_values,Teuchos::null,ap.gd,build_transient_support,0.0));

    const double DENSITY_VALUE = 3.7;
    const double TEMPERATURE_VALUE = 2.0;
    const double PERTURB_VALUE = 0.1;

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),DENSITY_VALUE);
      pIndex = me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
                                           ap.param_ged,param_density,ap.param_dofManager);
    }

    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                   ap.closure_models,
                   ap.user_data,false,"");

    // panzer::printMeshTopology(out,*ap.dofManager);


    RCP<OperatorType> D2fDpDx = me->create_DfDp_op(pIndex);

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dx = Thyra::createMember(*me->get_x_space());
      Thyra::assign(dx.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0/0.1);
      in_args.set_beta(1.0);

      me->evalModel_D2fDpDx(pIndex,in_args,dx,D2fDpDx);

      out << "D2fDpDx = \n" << Teuchos::describe(*D2fDpDx,Teuchos::VERB_EXTREME) << std::endl;
    }

    RCP<OperatorType> W = me->create_DfDp_op(pIndex);
    {
      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      InArgs in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0);
      in_args.set_beta(0.0);

      OutArgs out_args = me->createOutArgs();
      out_args.set_DfDp(pIndex,W);

      me->evalModel(in_args,out_args);

      out << "W = \n" << Teuchos::describe(*W,Teuchos::VERB_EXTREME) << std::endl;
    }

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-15);
    tester.num_random_vectors(20);

    double scaling_of_dfdp = 2.0*PERTURB_VALUE/TEMPERATURE_VALUE;

    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
    const bool op_cmp = tester.compare( *Thyra::scaleAndAdjoint(scaling_of_dfdp,Thyra::NOTRANS,W.getConst()), *D2fDpDx, Teuchos::ptrFromRef(fout));
    TEST_ASSERT(op_cmp);
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator_blocked_hessians, d2f_dp2)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
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

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    bool build_transient_support = true;
    RCP<PME> me
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,p_values,Teuchos::null,ap.gd,build_transient_support,0.0));

    const double DENSITY_VALUE = 3.7;
    const double TEMPERATURE_VALUE = 2.0;
    const double PERTURB_VALUE = 0.1;

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),DENSITY_VALUE);
      pIndex = me->addDistributedParameter("DENSITY_P",th_param_lof->getThyraDomainSpace(),
                                           ap.param_ged,param_density,ap.param_dofManager);
    }

    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                   ap.closure_models,
                   ap.user_data,false,"");

    // panzer::printMeshTopology(out,*ap.dofManager);


    RCP<OperatorType> D2fDp2 = me->create_DfDp_op(pIndex);

    // test hessian
    {

      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      RCP<VectorType> dp = Thyra::createMember(*me->get_p_space(0));
      Thyra::assign(dp.ptr(),PERTURB_VALUE);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0/0.1);
      in_args.set_beta(1.0);

      me->evalModel_D2fDp2(pIndex,in_args,dp,D2fDp2);

      out << "D2fDp2 = \n" << Teuchos::describe(*D2fDp2,Teuchos::VERB_EXTREME) << std::endl;
    }

    RCP<OperatorType> W = me->create_DfDp_op(pIndex);
    {
      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),TEMPERATURE_VALUE);

      InArgs in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_alpha(1.0);
      in_args.set_beta(0.0);

      OutArgs out_args = me->createOutArgs();
      out_args.set_DfDp(pIndex,W);

      me->evalModel(in_args,out_args);

      out << "W = \n" << Teuchos::describe(*W,Teuchos::VERB_EXTREME) << std::endl;
    }

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-15);
    tester.num_random_vectors(20);

    double scaling_of_dfdp = PERTURB_VALUE/DENSITY_VALUE;

    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
    const bool op_cmp = tester.compare( *Thyra::scaleAndAdjoint(scaling_of_dfdp,Thyra::NOTRANS,W.getConst()), *D2fDp2, Teuchos::ptrFromRef(fout));
    TEST_ASSERT(op_cmp);
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

    RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm
       = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Comm);

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    ap.fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
    ap.bcs = bcs;

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
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<ap.physicsBlocks.size();i++)
      wkstContainer->setNeeds(ap.physicsBlocks[i]->elementBlockID(),ap.physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);
    ap.wkstContainer = wkstContainer;

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // build the state dof manager and LOF
    {
      panzer::BlockedDOFManagerFactory globalIndexerFactory;
      auto dofManager = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),ap.physicsBlocks,
                                                                      conn_manager,"blocked: TEMPERATURE ION_TEMPERATURE");
      ap.dofManager = dofManager;

      Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
        = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(mpiComm,dofManager));
      ap.lof = linObjFactory;
    }

    // build the dof manager and LOF for DENSITY control
    if(distr_parameter_on) {
      Teuchos::RCP<panzer::DOFManager> dofManager
          = Teuchos::rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

      Teuchos::RCP<Intrepid2FieldPattern> fp
        = Teuchos::rcp(new Intrepid2FieldPattern(panzer::createIntrepid2Basis<PHX::exec_space,double,double>("HGrad",1,mesh->getCellTopology("eblock-0_0"))));
      dofManager->addField("eblock-0_0","DENSITY_P",fp);
      dofManager->addField("eblock-1_0","DENSITY_P",fp);

      dofManager->setOrientationsRequired(false);
      dofManager->buildGlobalUnknowns();

      // build a nonsquare LOF for the parameter vector
      Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > linObjFactory
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(mpiComm,ap.dofManager,dofManager));

      Teuchos::RCP<Epetra_Map> ownedMap = linObjFactory->getColMap(0);
      Teuchos::RCP<Epetra_Map> ghostedMap = linObjFactory->getGhostedColMap(0);
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*ghostedMap,*ownedMap));

      ap.param_dofManager = dofManager;
      ap.param_ged = Teuchos::rcp(new EpetraVector_ReadOnly_GlobalEvaluationData(importer,ghostedMap,ownedMap));
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
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<std::string>("Type","Product");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<std::string>("Term Names","TEMPERATURE,TEMPERATURE,DENSITY_P,DENSITY_P");
    closure_models.sublist("solid").sublist("DENSITY").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("HEAT_CAPACITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<std::string>("Type","Product");
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<std::string>("Term Names","ION_TEMPERATURE,ION_TEMPERATURE,DENSITY_P,DENSITY_P");
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_HEAT_CAPACITY").set<double>("Value",1.0);

    closure_models.sublist("solid").sublist("DENSITY_P").set("Type","Distributed Parameter");
    closure_models.sublist("solid").sublist("INTEGRAND").set<std::string>("Type","Product");
    closure_models.sublist("solid").sublist("INTEGRAND").set<std::string>("Term Names","TEMPERATURE,TEMPERATURE,DENSITY_P,DENSITY_P,DENSITY_P");

    ap.closure_models = closure_models;

    ap.user_data = Teuchos::ParameterList("User Data");

    ap.fmb->setWorksetContainer(wkstContainer);
    ap.fmb->setupVolumeFieldManagers(ap.physicsBlocks,ap.cm_factory,closure_models,*ap.lof,ap.user_data);
    ap.fmb->setupBCFieldManagers(bcs,ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,*ap.bc_factory,
                                 closure_models,*ap.lof,ap.user_data);
  }


}
