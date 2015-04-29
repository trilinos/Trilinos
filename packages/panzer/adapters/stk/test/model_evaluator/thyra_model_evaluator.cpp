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

#include "Phalanx_KokkosUtilities.hpp"

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

  struct RespFactoryFunc_Builder {
    MPI_Comm comm;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
    Teuchos::RCP<const panzer::UniqueGlobalIndexer<int,int> > globalIndexer;

    template <typename T>
    Teuchos::RCP<ResponseEvaluatorFactoryBase> build() const
    { return Teuchos::rcp(new ResponseEvaluatorFactory_Functional<T,int,int>(comm,1,true,"",linearObjFactory,globalIndexer)); }
  };

  void buildAssemblyPieces(bool parameter_on,bool distr_parameter_on,
                           AssemblyPieces & ap);

  bool testEqualityOfVectorValues(const Thyra::VectorBase<double> & a, 
                                  const Thyra::VectorBase<double> & b, 
                                  double tolerance, bool write_to_cout=false);

  TEUCHOS_UNIT_TEST(thyra_model_evaluator, basic)
  {
    using Teuchos::RCP;

    PHX::KokkosDeviceSession session;

    bool parameter_on = true;
    AssemblyPieces ap;
  
    buildAssemblyPieces(parameter_on,false,ap);

    // Test a transient me
    {
      typedef Thyra::ModelEvaluatorBase MEB;
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef Thyra::VectorBase<double> VectorType;
      typedef Thyra::LinearOpBase<double> OperatorType;
      typedef panzer::ModelEvaluator<double> PME;

      //std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = true;
    
      // RCP<PME> me = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,Teuchos::null,ap.gd,build_transient_support,0.0));
      RCP<PME> me = Teuchos::rcp(new PME(ap.lof,Teuchos::null,ap.gd,build_transient_support,0.0));
      me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                     *ap.eqset_factory,
                     *ap.bc_factory,
                     ap.cm_factory,
                     ap.cm_factory,
                     ap.closure_models,
                     ap.user_data,false,"");

      InArgs in_args = me->createInArgs();
      OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(MEB::IN_ARG_x));
      TEST_ASSERT(in_args.supports(MEB::IN_ARG_x_dot));
      TEST_ASSERT(in_args.supports(MEB::IN_ARG_alpha));
      TEST_ASSERT(in_args.supports(MEB::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(MEB::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(MEB::OUT_ARG_W_op));

      InArgs nomValues = me->getNominalValues();
      RCP<const VectorType> x = nomValues.get_x();
      RCP<VectorType> x_dot = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x_dot.ptr(),0.0);
      in_args.set_x(x);
      in_args.set_x_dot(x_dot);
      in_args.set_alpha(0.0);
      in_args.set_beta(1.0);
      
      RCP<VectorType> f = Thyra::createMember(*me->get_f_space());
      RCP<OperatorType> J_tmp = me->create_W_op();
      out_args.set_f(f);
      out_args.set_W_op(J_tmp);

      me->evalModel(in_args, out_args);
    }
  }

  TEUCHOS_UNIT_TEST(thyra_model_evaluator, response)
  {
    PHX::KokkosDeviceSession session;

    bool parameter_on = true;
    AssemblyPieces ap;
  
    buildAssemblyPieces(parameter_on,false,ap);

    {
      typedef Thyra::ModelEvaluatorBase MEB;
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef Thyra::VectorBase<double> VectorType;
      typedef Thyra::LinearOpBase<double> OperatorType;
      typedef panzer::ModelEvaluator<double> PME;

      // std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = false;
      // RCP<PME> me = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,Teuchos::null,ap.gd,build_transient_support,0.0));
      RCP<PME> me = Teuchos::rcp(new PME(ap.lof,Teuchos::null,ap.gd,build_transient_support,0.0));
 
      // parameterize the builder
      Teuchos::RCP<panzer::FunctionalResponse_Builder<int,int> > builder
        = Teuchos::rcp(new panzer::FunctionalResponse_Builder<int,int>);
 
      builder->comm = MPI_COMM_WORLD; // good enough
      builder->cubatureDegree = 1;
      builder->requiresCellIntegral = true;
      builder->quadPointField = "";

      std::vector<panzer::WorksetDescriptor> blocks;
      blocks.push_back(panzer::blockDescriptor("eblock-0_0"));
      blocks.push_back(panzer::blockDescriptor("eblock-1_0"));
      me->addFlexibleResponse("TEMPERATURE",blocks,builder);
      me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                     *ap.eqset_factory,
                     *ap.bc_factory,
                     ap.cm_factory,
                     ap.cm_factory,
                     ap.closure_models,
                     ap.user_data,false,"");

      InArgs nomValues = me->getNominalValues();
      RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
      Thyra::assign(x.ptr(),0.0);

      InArgs  in_args = me->createInArgs();
      in_args.set_x(x);
      
      RCP<VectorType> f = Thyra::createMember(*me->get_f_space());
      RCP<VectorType> g = Thyra::createMember(*me->get_g_space(0)); 
      RCP<VectorType> DgDx = Thyra::createMember(*me->get_x_space());
      OutArgs out_args = me->createOutArgs();
      out_args.set_f(f);
      out_args.set_g(0,g);
      out_args.set_DgDx(0,MEB::Derivative<double>(DgDx,MEB::DERIV_MV_GRADIENT_FORM));

      me->evalModel(in_args, out_args);

      out << "RESIDUAL = " << std::endl;
      out << Teuchos::describe(*f,Teuchos::VERB_EXTREME) << std::endl;

      out << "RESPONSE = " << std::endl;
      out << Teuchos::describe(*g,Teuchos::VERB_EXTREME) << std::endl;

      out << "RESPONSE DERIVATIVE = " << std::endl;
      out << Teuchos::describe(*DgDx,Teuchos::VERB_EXTREME) << std::endl;
    }
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(thyra_model_evaluator, scalar_parameters)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    PHX::KokkosDeviceSession session;

    bool parameter_on = true;
    AssemblyPieces ap;
  
    buildAssemblyPieces(parameter_on,false,ap);

    panzer::registerScalarParameter("DUMMY",*ap.gd->pl,3.0);
    panzer::registerScalarParameter("DUMMY_A",*ap.gd->pl,4.0);
    panzer::registerScalarParameter("DUMMY_B",*ap.gd->pl,5.0);

    {
      typedef Thyra::ModelEvaluatorBase MEB;
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef Thyra::VectorBase<double> VectorType;
      typedef Thyra::LinearOpBase<double> OperatorType;
      typedef Thyra::SpmdVectorBase<double> SpmdVector;
      typedef panzer::ModelEvaluator<double> PME;

      bool build_transient_support = false;
      RCP<PME> me = Teuchos::rcp(new PME(ap.lof,Teuchos::null,ap.gd,build_transient_support,0.0));
      me->addParameter("SOURCE_TEMPERATURE");
      me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                     *ap.eqset_factory,
                     *ap.bc_factory,
                     ap.cm_factory,
                     ap.cm_factory,
                     ap.closure_models,
                     ap.user_data,false,"");

      Teuchos::Array<std::string> params;
      params.push_back("DUMMY_A");
      params.push_back("DUMMY_B");

      int index_dummy = me->addParameter("DUMMY");
      int index_dummy_pair = me->addParameter(params);

      InArgs inArgs = me->createInArgs();

      TEST_EQUALITY(inArgs.Np(),3);

      RCP<const Teuchos::Array<std::string> > names_0 = me->get_p_names(0); 
      RCP<const Teuchos::Array<std::string> > names_1 = me->get_p_names(index_dummy); 
      RCP<const Teuchos::Array<std::string> > names_2 = me->get_p_names(index_dummy_pair); 

      TEST_ASSERT(names_0!=Teuchos::null);
      TEST_ASSERT(names_1!=Teuchos::null);
      TEST_ASSERT(names_2!=Teuchos::null);

      TEST_EQUALITY(names_0->size(),1);
      TEST_EQUALITY(names_1->size(),1);
      TEST_EQUALITY(names_2->size(),2);

      TEST_EQUALITY((*names_0)[0],"SOURCE_TEMPERATURE");
      TEST_EQUALITY((*names_1)[0],"DUMMY");
      TEST_EQUALITY((*names_2)[0],"DUMMY_A");
      TEST_EQUALITY((*names_2)[1],"DUMMY_B");

      RCP<const Thyra::VectorSpaceBase<double> > vs_0 = me->get_p_space(0);
      RCP<const Thyra::VectorSpaceBase<double> > vs_1 = me->get_p_space(index_dummy);
      RCP<const Thyra::VectorSpaceBase<double> > vs_2 = me->get_p_space(index_dummy_pair);

      TEST_THROW(me->get_p_space(-1),std::runtime_error);
      TEST_THROW(me->get_p_space(3),std::runtime_error);

      TEST_THROW(me->get_p_names(-1),std::runtime_error);
      TEST_THROW(me->get_p_names(3),std::runtime_error);

      TEST_ASSERT(vs_0!=Teuchos::null);
      TEST_ASSERT(vs_1!=Teuchos::null);
      TEST_ASSERT(vs_2!=Teuchos::null);

      TEST_EQUALITY(vs_0->dim(),1);
      TEST_EQUALITY(vs_1->dim(),1);
      TEST_EQUALITY(vs_2->dim(),2);


      RCP<const Thyra::VectorBase<double> > param_0 = inArgs.get_p(0);
      RCP<const Thyra::VectorBase<double> > param_1 = inArgs.get_p(index_dummy);
      RCP<const Thyra::VectorBase<double> > param_2 = inArgs.get_p(index_dummy_pair);

      TEST_ASSERT(param_0!=Teuchos::null);
      TEST_ASSERT(param_1!=Teuchos::null);
      TEST_ASSERT(param_2!=Teuchos::null);

      RCP<const SpmdVector> spmd_0 = rcp_dynamic_cast<const SpmdVector>(param_0);
      RCP<const SpmdVector> spmd_1 = rcp_dynamic_cast<const SpmdVector>(param_1);
      RCP<const SpmdVector> spmd_2 = rcp_dynamic_cast<const SpmdVector>(param_2);
      
      TEST_ASSERT(spmd_0!=Teuchos::null);
      TEST_ASSERT(spmd_1!=Teuchos::null);
      TEST_ASSERT(spmd_2!=Teuchos::null);

      {
        Teuchos::ArrayRCP<const double> data;
        spmd_0->getLocalData(Teuchos::ptrFromRef(data));

        TEST_EQUALITY(data.size(),1);
        TEST_EQUALITY(data[0],1.0);
      }

      {
        Teuchos::ArrayRCP<const double> data;
        spmd_1->getLocalData(Teuchos::ptrFromRef(data));

        TEST_EQUALITY(data.size(),1);
        TEST_EQUALITY(data[0],3.0);
      }

      {
        Teuchos::ArrayRCP<const double> data;
        spmd_2->getLocalData(Teuchos::ptrFromRef(data));

        TEST_EQUALITY(data.size(),2);
        TEST_EQUALITY(data[0],4.0);
        TEST_EQUALITY(data[1],5.0);
      }
    }
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(thyra_model_evaluator, scalar_parameters_dfdp)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    PHX::KokkosDeviceSession session;

    bool parameter_on = true;
    AssemblyPieces ap;
  
    buildAssemblyPieces(parameter_on,false,ap);

    {
      typedef Thyra::ModelEvaluatorBase MEB;
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef Thyra::VectorBase<double> VectorType;
      typedef Thyra::LinearOpBase<double> OperatorType;
      typedef panzer::ModelEvaluator<double> PME;

      bool build_transient_support = false;
      RCP<PME> me = Teuchos::rcp(new PME(ap.lof,Teuchos::null,ap.gd,build_transient_support,0.0));
      me->addParameter("SOURCE_TEMPERATURE");
      me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                     *ap.eqset_factory,
                     *ap.bc_factory,
                     ap.cm_factory,
                     ap.cm_factory,
                      ap.closure_models,
                      ap.user_data,false,"");

      RCP<Thyra::VectorBase<double> > x = Thyra::createMember(me->get_x_space());
      Thyra::put_scalar(1.0,x.ptr());

      RCP<Thyra::VectorBase<double> > p = Thyra::createMember(me->get_p_space(0));
      Thyra::put_scalar(1.0,p.ptr());

      RCP<Thyra::VectorBase<double> > f1 = Thyra::createMember(*me->get_f_space());
      RCP<Thyra::VectorBase<double> > f2 = Thyra::createMember(*me->get_f_space());
      RCP<Thyra::VectorBase<double> > f3 = Thyra::createMember(*me->get_f_space());
      RCP<Thyra::VectorBase<double> > f4 = Thyra::createMember(*me->get_f_space());

      InArgs inArgs = me->createInArgs();
      inArgs.set_x(x);
      inArgs.set_p(0,p);

      OutArgs outArgs = me->createOutArgs();

      out << "evalModel(f1)" << std::endl;
      Thyra::put_scalar(1.0,p.ptr());
      outArgs.set_f(f1);
      outArgs.set_DfDp(0,MEB::Derivative<double>());
      me->evalModel(inArgs,outArgs);
      
      out << "evalModel(f2)" << std::endl;
      Thyra::put_scalar(1.0,p.ptr());
      outArgs.set_f(f2);
      outArgs.set_DfDp(0,MEB::Derivative<double>());
      me->evalModel(inArgs,outArgs);
      
      out << "evalModel(f3)" << std::endl;
      Thyra::put_scalar(20.0,p.ptr());
      outArgs.set_f(f3);
      outArgs.set_DfDp(0,MEB::Derivative<double>());
      me->evalModel(inArgs,outArgs);
      
      out << "evalModel(f4)" << std::endl;
      Thyra::put_scalar(1.0,p.ptr());
      outArgs.set_f(f4);
      outArgs.set_DfDp(0,MEB::Derivative<double>());
      me->evalModel(inArgs,outArgs);

      double tol = 10.0 * Teuchos::ScalarTraits<double>::eps();

      // f1 == f2
      TEST_EQUALITY_CONST(testEqualityOfVectorValues(*f1,*f2,tol), true);

      // f2 == f4
      TEST_EQUALITY_CONST(testEqualityOfVectorValues(*f2,*f4,tol), true);

      // f2 != f3
      TEST_EQUALITY_CONST(testEqualityOfVectorValues(*f2,*f3,tol), false);

      // TEST DfDp
      /////////////////////////////////////////////////////

      RCP<Thyra::VectorBase<double> > dfdp = Thyra::createMember(*me->get_f_space());
      Thyra::put_scalar(0.0,dfdp.ptr());
  
      Thyra::put_scalar(0.0,x.ptr());
      Thyra::put_scalar(0.0,f1.ptr());
      Thyra::put_scalar(20.0,p.ptr());
  
      out << "evalModel(dfdp)" << std::endl;
      outArgs.set_f(f1);
      outArgs.set_DfDp(0,MEB::Derivative<double>(dfdp,MEB::DERIV_MV_BY_COL));
      me->evalModel(inArgs,outArgs);
  
      Teuchos::ArrayRCP<const double> dfdp_data;
      Teuchos::ArrayRCP<const double> f1_data;
      dynamic_cast<const Thyra::SpmdVectorBase<double> &>(*dfdp).getLocalData(Teuchos::ptrFromRef(dfdp_data));
      dynamic_cast<const Thyra::SpmdVectorBase<double> &>(*f1).getLocalData(Teuchos::ptrFromRef(f1_data));
      for(int i=0;i<dfdp_data.size();i++) {
        if(dfdp_data[i]!=0.0)
        { TEST_FLOATING_EQUALITY(f1_data[i],20.0*dfdp_data[i],1e-10); }
        out << f1_data[i] << "    " << dfdp_data[i] << std::endl;
      }
    }
  }

  // Testing Ditributed Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator, distributed_parameters)
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
    typedef panzer::ModelEvaluator<double> PME;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    PHX::KokkosDeviceSession session;

    bool parameter_on = true;
    AssemblyPieces ap;

    buildAssemblyPieces(parameter_on,false,ap);

    int distributed_parameter_index = -1;
    Teuchos::RCP<panzer::LOCPair_GlobalEvaluationData> dataObject;
    RCP<PME> me;
    {
      bool build_transient_support = false;
      me = Teuchos::rcp(new PME(ap.lof,Teuchos::null,ap.gd,build_transient_support,0.0));
      me->addParameter("SOURCE_TEMPERATURE");
      
      // add a distributed parameter
      {
        // extract the vector space for the distributed parameter
        Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs 
            = Teuchos::rcp_dynamic_cast<const ThyraObjFactory<double> >(ap.lof,true)->getThyraDomainSpace();

        // setup an initial value for the parameter
        Teuchos::RCP<Thyra::VectorBase<double> > initial = Thyra::createMember(vs);
        Thyra::put_scalar(0.0,initial.ptr());

        // this object changes the parameter from a global object to a ghosted one
        dataObject = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(ap.lof,panzer::LinearObjContainer::X));

        // add the distributed parameter
        distributed_parameter_index = me->addDistributedParameter("Transient Predictor",
                                                                  vs,
                                                                  dataObject,
                                                                  initial); 
      }

      me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                     *ap.eqset_factory,
                     *ap.bc_factory,
                     ap.cm_factory,
                     ap.cm_factory,
                      ap.closure_models,
                      ap.user_data,false,"");
    }

    // build inputs
    ////////////////////////////////////////////////////////////////////////////////////

    // solution
    RCP<Thyra::VectorBase<double> > x = Thyra::createMember(me->get_x_space());
    Thyra::put_scalar(1.0,x.ptr());
    
    // locally replicated scalar parameter
    RCP<Thyra::VectorBase<double> > p = Thyra::createMember(me->get_p_space(0));
    Thyra::put_scalar(1.0,p.ptr());

    // distributed parameter
    RCP<Thyra::VectorBase<double> > distr_p = Thyra::createMember(me->get_p_space(distributed_parameter_index));
    Thyra::put_scalar(3.14,distr_p.ptr());

    // build outputs
    ////////////////////////////////////////////////////////////////////////////////////

    RCP<Thyra::VectorBase<double> > f = Thyra::createMember(me->get_f_space());

    // Test that the distributed parameter is updated correctly
    {
      InArgs in_args = me->createInArgs();
      OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.Np() == 2);

      in_args.set_x(x);
      in_args.set_p(0,p);
      in_args.set_p(distributed_parameter_index,distr_p);

      out_args.set_f(f);

      me->evalModel(in_args,out_args);
    
      // Export should have performed global to ghost, ghosted values should be 1
      // Create a gold standard to compare against
      RCP<const Thyra::VectorBase<double> > ghosted_distr_p 
          = Teuchos::rcp_dynamic_cast<ThyraObjContainer<double> >(dataObject->getGhostedLOC())->get_x_th();
      RCP<Thyra::VectorBase<double> > gold_standard = Thyra::createMember(ghosted_distr_p->range());
      Thyra::put_scalar(3.14,gold_standard.ptr());

      double tol = 10.0 * Teuchos::ScalarTraits<double>::eps();

      // note that all this tests is that the ghosted vector is correctly populated!
      TEST_EQUALITY_CONST(testEqualityOfVectorValues(*ghosted_distr_p,*gold_standard,tol,true), true);
    }

    
  }

  // Testing Parameter Support
  TEUCHOS_UNIT_TEST(thyra_model_evaluator, distro_parameters_dgdp)
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

    PHX::KokkosDeviceSession session;

    bool parameter_on = true;
    bool distr_param_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(parameter_on,distr_param_on,ap);

    int pIndex = -1;
    int rIndex = -1;

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    bool build_transient_support = true;
    RCP<PME> me 
        = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,Teuchos::null,ap.gd,build_transient_support,0.0));

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
      rIndex = me->addFlexibleResponse("DENSITY",blocks,builder); // integrate the density
    }

    // add distributed parameter
    {
      RCP<ThyraObjFactory<double> > th_param_lof = rcp_dynamic_cast<ThyraObjFactory<double> >(ap.param_lof);

      RCP<VectorType> param_density = Thyra::createMember(th_param_lof->getThyraDomainSpace());
      Thyra::assign(param_density.ptr(),3.7);
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

    // check that out args support DgDp
    {
      OutArgs outArgs = me->createOutArgs();
      TEST_ASSERT(outArgs.supports(MEB::OUT_ARG_DgDp,rIndex,pIndex).supports(MEB::DERIV_MV_GRADIENT_FORM));
    }

    // solution
    RCP<Thyra::VectorBase<double> > x = Thyra::createMember(me->get_x_space());
    Thyra::put_scalar(1.0,x.ptr());

    InArgs  in_args = me->createInArgs();
      
    RCP<VectorType> g = Thyra::createMember(*me->get_g_space(rIndex)); 
    RCP<VectorType> DgDp = Thyra::createMember(*me->get_p_space(pIndex)); 

    out << "DgDp = " << std::endl;

    {
      OutArgs out_args = me->createOutArgs();
      out_args.set_g(0,g);

      me->evalModel(in_args, out_args);

      out << "RESPONSE = " << std::endl;
      out << Teuchos::describe(*g,Teuchos::VERB_EXTREME) << std::endl;

      Teuchos::ArrayRCP<const double> g_data;
      rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(g)->getLocalData(Teuchos::ptrFromRef(g_data));
      TEST_FLOATING_EQUALITY(g_data[0],3.7,1e-12);
    }

    {
      OutArgs out_args = me->createOutArgs();
      out_args.set_g(0,g);
      out_args.set_DgDp(rIndex,pIndex,MEB::Derivative<double>(DgDp,MEB::DERIV_MV_GRADIENT_FORM));

      me->evalModel(in_args, out_args);

      out << "RESPONSE = " << std::endl;
      out << Teuchos::describe(*g,Teuchos::VERB_EXTREME) << std::endl;

      Teuchos::ArrayRCP<const double> g_data;
      rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(g)->getLocalData(Teuchos::ptrFromRef(g_data));
      TEST_FLOATING_EQUALITY(g_data[0],3.7,1e-12);
    }

    Teuchos::ArrayRCP<const double> DgDp_data;
    RCP<const SpmdVectorType> spmd_DgDp = rcp_dynamic_cast<const SpmdVectorType>(DgDp);
    rcp_dynamic_cast<const SpmdVectorType>(DgDp)->getLocalData(Teuchos::ptrFromRef(DgDp_data));
    double int_phi = 1.0/192.0;
    for(int i=0;i<DgDp_data.size();i++) {
      out << DgDp_data[i]  << " " << int_phi << std::endl;
      bool a = std::fabs(DgDp_data[i]-int_phi)/int_phi            <= 1e-14;
      bool b = std::fabs(DgDp_data[i]-2.0*int_phi)/(2.0*int_phi)  <= 1e-14;
      bool c = std::fabs(DgDp_data[i]-4.0*int_phi)/(4.0*int_phi)  <= 1e-14;

      TEST_ASSERT(a || b || c);
    }
  }

  // Testing that nominal values are correctly built and initialized
  //    specifically testing that adding distributed parameters doesn't wipe out
  //    previously set nominal values (like the inital condition)
  TEUCHOS_UNIT_TEST(model_evaluator, nominal_values)
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
    typedef panzer::ModelEvaluator<double> PME;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;

    PHX::KokkosDeviceSession session;

    double tol = 10.0 * Teuchos::ScalarTraits<double>::eps();

    bool parameter_on = true;
    AssemblyPieces ap;

    buildAssemblyPieces(parameter_on,false,ap);

    Teuchos::RCP<panzer::LOCPair_GlobalEvaluationData> dataObject;
    RCP<PME> me;

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;

    me = Teuchos::rcp(new PME(ap.lof,Teuchos::null,ap.gd,true,0.0));
    me->setupModel(ap.wkstContainer,ap.physicsBlocks,ap.bcs,
                   *ap.eqset_factory,
                   *ap.bc_factory,
                   ap.cm_factory,
                   ap.cm_factory,
                    ap.closure_models,
                    ap.user_data,false,"");

    // setup some initial conditions
    RCP<Thyra::VectorBase<double> > x0     = Thyra::createMember(me->get_x_space());
    RCP<Thyra::VectorBase<double> > x0_dot = Thyra::createMember(me->get_x_space());
    {
      InArgs nomValues = me->getNominalValues();

      // why is this const cast neccessary. What is "correct" approach to setting initial conditions
      Thyra::put_scalar(3.14,rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_x()).ptr());
      Thyra::put_scalar(6.28,rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_x_dot()).ptr());

      // setup the comparison vectors
      Thyra::put_scalar(3.14,x0.ptr());
      Thyra::put_scalar(6.28,x0_dot.ptr());
    }

    // check the initial conditions
    {
      InArgs nomValues = me->getNominalValues();

      RCP<const Thyra::VectorBase<double> > x0_n     = nomValues.get_x();
      RCP<const Thyra::VectorBase<double> > x0_dot_n = nomValues.get_x_dot();

      // check the error of x0
      TEST_ASSERT(testEqualityOfVectorValues(    *x0_n,     *x0, tol));
      TEST_ASSERT(testEqualityOfVectorValues(*x0_dot_n, *x0_dot, tol));
    }
 
    // add normal parameter
    /////////////////////////////////////////////////////////////
    int index_dummy = me->addParameter("SOURCE_TEMPERATURE");

    // check the initial conditions
    {
      InArgs nomValues = me->getNominalValues();

      RCP<const Thyra::VectorBase<double> > x0_n     = nomValues.get_x();
      RCP<const Thyra::VectorBase<double> > x0_dot_n = nomValues.get_x_dot();

      // check the error of x0
      TEST_ASSERT(testEqualityOfVectorValues(    *x0_n,     *x0, tol));
      TEST_ASSERT(testEqualityOfVectorValues(*x0_dot_n, *x0_dot, tol));
    }

    // set new parameter value
    RCP<Thyra::VectorBase<double> > p0     = Thyra::createMember(me->get_p_space(index_dummy));
    {
      InArgs nomValues = me->getNominalValues();

      Thyra::put_scalar(9.42,rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_p(index_dummy)).ptr());
      Thyra::put_scalar(9.42,p0.ptr());
    }

    // add first distributed parameter
    /////////////////////////////////////////////////////////////
    int index_distr_param1 = -1;
    {
      // the vector space for the distributed parameter
      Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs = me->get_x_space();
  
      // setup an initial value for the parameter
      Teuchos::RCP<Thyra::VectorBase<double> > initial = Thyra::createMember(vs);
      Thyra::put_scalar(0.0,initial.ptr());
  
      // this object changes the parameter from a global object to a ghosted one
      dataObject = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(ap.lof,panzer::LinearObjContainer::X));
  
      index_distr_param1 = me->addDistributedParameter("Distr Parameter 1",
                                                       vs,
                                                       dataObject,
                                                       initial); 
    }

    // check the initial conditions
    {
      InArgs nomValues = me->getNominalValues();

      RCP<const Thyra::VectorBase<double> > x0_n     = nomValues.get_x();
      RCP<const Thyra::VectorBase<double> > x0_dot_n = nomValues.get_x_dot();
      RCP<const Thyra::VectorBase<double> > p0_n     = nomValues.get_p(index_dummy);

      // check the error of x0
      TEST_ASSERT(testEqualityOfVectorValues(    *x0_n,     *x0, tol));
      TEST_ASSERT(testEqualityOfVectorValues(*x0_dot_n, *x0_dot, tol));
      TEST_ASSERT(testEqualityOfVectorValues(    *p0_n,     *p0, tol));
    }

    // set new parameter value
    RCP<Thyra::VectorBase<double> > p1     = Thyra::createMember(me->get_p_space(index_distr_param1));
    {
      typedef Thyra::VectorBase<double> Vector;

      InArgs nomValues = me->getNominalValues();

      Thyra::put_scalar(-3.14,rcp_const_cast<Vector>(nomValues.get_p(index_distr_param1)).ptr());
      Thyra::put_scalar(-3.14,p1.ptr());
    }

    // add second distributed parameter
    /////////////////////////////////////////////////////////////
    {
      // the vector space for the distributed parameter
      Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs = me->get_x_space();
  
      // setup an initial value for the parameter
      Teuchos::RCP<Thyra::VectorBase<double> > initial = Thyra::createMember(vs);
      Thyra::put_scalar(0.0,initial.ptr());
  
      // this object changes the parameter from a global object to a ghosted one
      dataObject = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(ap.lof,panzer::LinearObjContainer::X));
  
      // we are really just interested in what happens when we call this
      me->addDistributedParameter("Distr Parameter 2",
                                  vs,
                                  dataObject,
                                  initial); 
    }

    // check the initial conditions
    {
      InArgs nomValues = me->getNominalValues();

      RCP<const Thyra::VectorBase<double> > x0_n     = nomValues.get_x();
      RCP<const Thyra::VectorBase<double> > x0_dot_n = nomValues.get_x_dot();
      RCP<const Thyra::VectorBase<double> > p0_n     = nomValues.get_p(index_dummy);
      RCP<const Thyra::VectorBase<double> > p1_n     = nomValues.get_p(index_distr_param1);

      // check the error of x0
      TEST_ASSERT(testEqualityOfVectorValues(    *x0_n,     *x0, tol));
      TEST_ASSERT(testEqualityOfVectorValues(*x0_dot_n, *x0_dot, tol));
      TEST_ASSERT(testEqualityOfVectorValues(    *p0_n,     *p0, tol));
      TEST_ASSERT(testEqualityOfVectorValues(    *p1_n,     *p1, tol));
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
  
  void buildAssemblyPieces(bool parameter_on,bool distr_parameter_on,
                           AssemblyPieces & ap)
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
    ap.fmb->setupBCFieldManagers(bcs,ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,*ap.bc_factory,
                                 closure_models,*ap.lof,ap.user_data);
  }


}
