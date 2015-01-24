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
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#ifdef HAVE_STOKHOS
   #include "Panzer_SGEpetraLinearObjFactory.hpp"
 
   #include "Stokhos_HermiteBasis.hpp"
   #include "Stokhos_CompletePolynomialBasis.hpp"
   #include "Stokhos_QuadOrthogPolyExpansion.hpp"
   #include "Stokhos_TensorProductQuadrature.hpp"
   #include "Stokhos_SGModelEvaluator.hpp"
#endif

#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <cstdio> // for get char
#include <fstream>

namespace panzer {

  struct AssemblyPieces {
    RCP<panzer::FieldManagerBuilder> fmb;  
    RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    RCP<panzer::GlobalData> gd;
    RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager;
    Teuchos::ParameterList user_data;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::EquationSetFactory> eqset_factory;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    Teuchos::ParameterList closure_models;
  };

  RCP<Epetra_Vector> basic_ss_f;
  RCP<Epetra_Vector> basic_trans_f;
  RCP<Epetra_CrsMatrix> basic_ss_J;
  RCP<Epetra_CrsMatrix> basic_trans_J;

  #ifdef HAVE_STOKHOS
  Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> > sg_lof_null;
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > sg_exp_null;
  #endif

  struct RespFactoryFunc_Builder {
    MPI_Comm comm;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
    Teuchos::RCP<const panzer::UniqueGlobalIndexer<int,int> > globalIndexer;

    template <typename T>
    Teuchos::RCP<ResponseEvaluatorFactoryBase> build() const
    { return Teuchos::rcp(new ResponseEvaluatorFactory_Functional<T,int,int>(comm,1,true,"",linearObjFactory,globalIndexer)); }
  };

  // store steady-state me for testing parameters
  // RCP<panzer::ModelEvaluator_Epetra> ss_me;
  Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager_null;

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  bool testEqualityOfEpetraVectorValues(Epetra_Vector& a, Epetra_Vector& b, double tolerance, bool write_to_cout = false);

  void buildAssemblyPieces(bool parameter_on,
                           AssemblyPieces & ap
                           #ifdef HAVE_STOKHOS
                           , Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> > & sg_lof=sg_lof_null
                           , const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & sgExpansion=sg_exp_null
                           #endif
                           );

  TEUCHOS_UNIT_TEST(model_evaluator, basic)
  {
    using Teuchos::RCP;

    // panzer::pauseToAttach();

    bool parameter_on = true;
    AssemblyPieces ap;
  
    buildAssemblyPieces(parameter_on,ap);

    // Test a transient me
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = true;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,ap.ep_lof,p_names,ap.gd,build_transient_support));

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_W));

      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_sg));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_sg));
      
      
      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      RCP<Epetra_Vector> x_dot = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      x_dot->PutScalar(0.0);
      in_args.set_x(x);
      in_args.set_x_dot(x_dot);
      in_args.set_alpha(0.0);
      in_args.set_beta(1.0);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Operator> J_tmp = me->create_W();
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J_tmp);
      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);

      // compare or save residual vector for comparison with SG version
      if(basic_trans_f==Teuchos::null)
         basic_trans_f = out_args.get_f();
      else {
         double norm=0.0, diff=0.0;
         f->Norm2(&norm);
         f->Update(-1.0,*basic_trans_f,1.0);    
         f->Norm2(&diff);
         TEST_ASSERT(diff/norm <= 1e-15);
      }
    }

    // Test a steady-state me
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names(1);
      p_names[0] = Teuchos::rcp(new Teuchos::Array<std::string>(1));
      (*p_names[0])[0] = "SOURCE_TEMPERATURE";
      bool build_transient_support = false;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,ap.ep_lof,p_names,ap.gd,build_transient_support));
      
      // store to test parameter capabilities
      // ss_me = me;

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_W));

      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_sg));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_sg));
      
      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      in_args.set_x(x);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Operator> J_tmp = me->create_W();
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J_tmp);
      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);

      // compare or save residual vector for comparison with SG version
      if(basic_ss_f==Teuchos::null)
         basic_ss_f = out_args.get_f();
      else {
         double norm=0.0, diff=0.0;
         f->Norm2(&norm);
         f->Update(-1.0,*basic_ss_f,1.0);    
         f->Norm2(&diff);
         out << diff << " " << norm << std::endl;
         TEST_ASSERT(diff/norm <= 1e-15);
      }
    }

  }

  TEUCHOS_UNIT_TEST(model_evaluator, response)
  {
    bool parameter_on = true;
    AssemblyPieces ap;
  
    buildAssemblyPieces(parameter_on,ap);

    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = false;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,ap.ep_lof,p_names,ap.gd,build_transient_support));
 
      RespFactoryFunc_Builder builder;
      builder.comm = MPI_COMM_WORLD;
      builder.linearObjFactory = ap.ep_lof;
      builder.globalIndexer = ap.dofManager;

      std::vector<panzer::WorksetDescriptor> blocks;
      blocks.push_back(panzer::blockDescriptor("eblock-0_0"));
      blocks.push_back(panzer::blockDescriptor("eblock-1_0"));
      me->addResponse("TEMPERATURE",blocks,builder);
      me->buildResponses(ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,ap.closure_models,ap.user_data);

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();

      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->PutScalar(0.0);
      in_args.set_x(x);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Vector> g = Teuchos::rcp(new Epetra_Vector(*me->get_g_map(0)));
      RCP<Epetra_Vector> DgDx = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      out_args.set_f(f);
      out_args.set_g(0,g);
      out_args.set_DgDx(0,EpetraExt::ModelEvaluator::Derivative(Teuchos::rcp_static_cast<Epetra_MultiVector>(DgDx)));

      me->evalModel(in_args, out_args);

      out << "RESPONSE = " << std::endl;
      g->Print(out);

      out << "RESIDUAL | RESPONSE DERIVATIVE = " << std::endl;
      for(int i=0;i<DgDx->MyLength();i++)
        out << (*f)[i] << "     " << (*DgDx)[i] << std::endl;
    }
  }
  
/*
  TEUCHOS_UNIT_TEST(model_evaluator, basic_parameter)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb_off;  
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary_off; 
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof_off;
    Teuchos::RCP<panzer::GlobalData> gd_off;
    buildAssemblyPieces(false,fmb_off,rLibrary_off,gd_off,ep_lof_off);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb_on;  
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary_on; 
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof_on;
    Teuchos::RCP<panzer::GlobalData> gd_on;
    buildAssemblyPieces(true,fmb_on,rLibrary_on,gd_on,ep_lof_on);

    // transient test
    {
      // build base line residual and Jacobian
      /////////////////////////////////////////
      RCP<Epetra_Vector> f_off;
      RCP<Epetra_CrsMatrix> J_off;
      {
         std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names_off;
         bool build_transient_support = true;
         RCP<panzer::ModelEvaluator_Epetra> me_off 
              = rcp(new panzer::ModelEvaluator_Epetra(fmb_off,rLibrary_off,ep_lof_off,p_names_off,gd_off,true));
   
         EpetraExt::ModelEvaluator::InArgs in_args = me_off->createInArgs();
         EpetraExt::ModelEvaluator::OutArgs out_args = me_off->createOutArgs();
         
         RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*me_off->get_x_map()));
         RCP<Epetra_Vector> x_dot = rcp(new Epetra_Vector(*me_off->get_x_map()));
         x->Update(1.0, *(me_off->get_x_init()), 0.0);
         x_dot->PutScalar(0.0);
         in_args.set_x(x);
         in_args.set_x_dot(x_dot);
         in_args.set_alpha(0.0);
         in_args.set_beta(1.0);
   
         f_off = rcp(new Epetra_Vector(*me_off->get_f_map()));
         J_off = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(me_off->create_W());
         out_args.set_f(f_off);
         out_args.set_W(J_off);
   
         me_off->evalModel(in_args, out_args);
      }

      // build test of "on" but no parameters registered
      /////////////////////////////////////////////

      {
         std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names_on;
         bool build_transient_support = true;
         RCP<panzer::ModelEvaluator_Epetra> me_on 
              = rcp(new panzer::ModelEvaluator_Epetra(fmb_on,rLibrary_on,ep_lof_on,p_names_on,gd_on,true));
   
         EpetraExt::ModelEvaluator::InArgs in_args = me_on->createInArgs();
         EpetraExt::ModelEvaluator::OutArgs out_args = me_on->createOutArgs();
         
         RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*me_on->get_x_map()));
         RCP<Epetra_Vector> x_dot = rcp(new Epetra_Vector(*me_on->get_x_map()));
         x->Update(1.0, *(me_on->get_x_init()), 0.0);
         x_dot->PutScalar(0.0);
         in_args.set_x(x);
         in_args.set_x_dot(x_dot);
         in_args.set_alpha(0.0);
         in_args.set_beta(1.0);
   
         RCP<Epetra_Vector> f = rcp(new Epetra_Vector(*me_on->get_f_map()));
         RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(me_on->create_W());
         out_args.set_f(f);
         out_args.set_W(J);
   
         me_on->evalModel(in_args, out_args);

         double norm_f=0.0, norm_f_off=0.0;
         f->Norm2(&norm_f);
         f_off->Norm2(&norm_f_off);
         TEST_FLOATING_EQUALITY(norm_f,norm_f_off,1e-15);
      }

      // build test of "on" but with a parameter registered but not modified
      /////////////////////////////////////////////

      {
         std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names_on;
         std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names(1);
         p_names[0] = Teuchos::rcp(new Teuchos::Array<std::string>(1));
         (*p_names[0])[0] = "SOURCE_TEMPERATURE";
         bool build_transient_support = true;
         RCP<panzer::ModelEvaluator_Epetra> me_on 
              = rcp(new panzer::ModelEvaluator_Epetra(fmb_on,rLibrary_on,ep_lof_on,p_names_on,gd_on,true));
   
         EpetraExt::ModelEvaluator::InArgs in_args = me_on->createInArgs();
         EpetraExt::ModelEvaluator::OutArgs out_args = me_on->createOutArgs();
         
         RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*me_on->get_x_map()));
         RCP<Epetra_Vector> x_dot = rcp(new Epetra_Vector(*me_on->get_x_map()));
         x->Update(1.0, *(me_on->get_x_init()), 0.0);
         x_dot->PutScalar(0.0);
         in_args.set_x(x);
         in_args.set_x_dot(x_dot);
         in_args.set_alpha(0.0);
         in_args.set_beta(1.0);
   
         RCP<Epetra_Vector> f = rcp(new Epetra_Vector(*me_on->get_f_map()));
         RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(me_on->create_W());
         out_args.set_f(f);
         out_args.set_W(J);
   
         me_on->evalModel(in_args, out_args);

         double norm_f=0.0, norm_f_off=0.0;
         f->Norm2(&norm_f);
         f_off->Norm2(&norm_f_off);
         TEST_FLOATING_EQUALITY(norm_f,norm_f_off,1e-15);
      }
    }
  }
*/

  // Test for parameters and residual consistency
  /////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(model_evaluator, parameter_library)
  {
    // NOTE: this test must be run AFTER the basic test above so that
    // ss_me is created!
    // RCP<panzer::ModelEvaluator_Epetra> me = ss_me; This appears to cause seg faults for some reason!!!!

    RCP<panzer::ModelEvaluator_Epetra> me;
    {
      bool parameter_on = true;
      AssemblyPieces ap;
    
      buildAssemblyPieces(parameter_on,ap);
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names(1);
      p_names[0] = Teuchos::rcp(new Teuchos::Array<std::string>(1));
      (*p_names[0])[0] = "SOURCE_TEMPERATURE";
      bool build_transient_support = false;
      me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,ap.ep_lof,p_names,ap.gd,build_transient_support));
    }

    TEUCHOS_ASSERT(nonnull(me));

    EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
    // solution
    RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
    x->PutScalar(1.0);
    
    // parameter
    TEST_ASSERT(in_args.Np() == 1);
    RCP<Epetra_Vector> p = Teuchos::rcp(new Epetra_Vector(*me->get_p_map(0)));
    p->PutScalar(1.0);

    // create f
    RCP<Epetra_Vector> f1 = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
    RCP<Epetra_Vector> f2 = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
    RCP<Epetra_Vector> f3 = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
    RCP<Epetra_Vector> f4 = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));

    RCP<Epetra_Vector> dfdp = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));

    // set values and evaluate
    in_args.set_x(x);
    in_args.set_p(0,p);

    out << "evalModel(f1)" << std::endl;
    out_args.set_f(f1);
    out_args.set_DfDp(0,EpetraExt::ModelEvaluator::Derivative());
    me->evalModel(in_args,out_args);
    
    out << "evalModel(f2)" << std::endl;
    out_args.set_f(f2);
    out_args.set_DfDp(0,EpetraExt::ModelEvaluator::Derivative());
    me->evalModel(in_args,out_args);
    
    out << "evalModel(f3)" << std::endl;
    p->PutScalar(20.0);
    out_args.set_f(f3);
    out_args.set_DfDp(0,EpetraExt::ModelEvaluator::Derivative());
    me->evalModel(in_args,out_args);
    
    out << "evalModel(f4)" << std::endl;
    p->PutScalar(1.0);
    out_args.set_f(f4);
    out_args.set_DfDp(0,EpetraExt::ModelEvaluator::Derivative());
    me->evalModel(in_args,out_args);
    
    // f1 == f2 == f4, f3 is evaluated with p=20 instead of p=1

    double tol = 10.0 * Teuchos::ScalarTraits<double>::eps();

    // f1 == f2
    // ROGER: This demonstrates the broken functionality!!!
    // Uncomment line below when fixed!!!
    TEST_EQUALITY_CONST(testEqualityOfEpetraVectorValues(*f1,*f2,tol), true);

    // f2 == f4
    TEST_EQUALITY_CONST(testEqualityOfEpetraVectorValues(*f2,*f4,tol), true);

    // f2 != f3
    TEST_EQUALITY_CONST(testEqualityOfEpetraVectorValues(*f2,*f3,tol), false);
 

    // TEST DfDp
    /////////////////////////////////////////////////////

    x->PutScalar(0.0);
    f1->PutScalar(0.0);

    out << "evalModel(f2)" << std::endl;
    p->PutScalar(20.0);
    out_args.set_f(f1);
    out_args.set_DfDp(0,EpetraExt::ModelEvaluator::Derivative(dfdp,EpetraExt::ModelEvaluator::DERIV_MV_BY_COL));
    me->evalModel(in_args,out_args);

    for(int i=0;i<f1->MyLength();i++) {
      if((*dfdp)[i]!=0.0)
      { TEST_FLOATING_EQUALITY((*f1)[i],20.0*(*dfdp)[i],1e-10); }
      out << (*f1)[i] << "    " << (*dfdp)[i] << std::endl;
    }

  }

  // Testing Ditributed Parameter Support
  TEUCHOS_UNIT_TEST(model_evaluator, distributed_parameters)
  {

    RCP<panzer::ModelEvaluator_Epetra> me;
    int distributed_parameter_index = -1;
    Teuchos::RCP<Epetra_Vector> ghosted_distributed_parameter;
    {
      bool parameter_on = true;
      AssemblyPieces ap;
    
      buildAssemblyPieces(parameter_on,ap);
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names(1);
      p_names[0] = Teuchos::rcp(new Teuchos::Array<std::string>(1));
      (*p_names[0])[0] = "SOURCE_TEMPERATURE";
      bool build_transient_support = false;
      me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,ap.ep_lof,p_names,ap.gd,build_transient_support));

      // add a distributed parameter
      ghosted_distributed_parameter = Teuchos::rcp(new Epetra_Vector(*ap.ep_lof->getGhostedMap()));

      distributed_parameter_index = me->addDistributedParameter("Transient Predictor",
								ap.ep_lof->getMap(),
								ap.ep_lof->getGhostedImport(),
								ghosted_distributed_parameter);
    }

    TEUCHOS_ASSERT(nonnull(me));

    // solution
    RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
    x->PutScalar(1.0);
    
    // locally replicated scalar parameter
    RCP<Epetra_Vector> p = Teuchos::rcp(new Epetra_Vector(*me->get_p_map(0)));
    p->PutScalar(1.0);

    RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));

    // Test that the distributed parameter is updated correctly
    {
      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.Np() == 2);

      in_args.set_x(x);
      in_args.set_p(0,p);

      // Set ghosted distributed parameter to 0
      ghosted_distributed_parameter->PutScalar(0.0);

      // Set global distributed parameter to 1 in in_args
      RCP<Epetra_Vector> global_distributed_parameter = 
	Teuchos::rcp(new Epetra_Vector(*me->get_p_map(distributed_parameter_index)));
      global_distributed_parameter->PutScalar(1.0);
      in_args.set_p(distributed_parameter_index,global_distributed_parameter);

      out_args.set_f(f);
      me->evalModel(in_args,out_args);
    
      // Export should have performed global to ghost, ghosted values should be 1
      // Create a gold standard to compare against
      RCP<Epetra_Vector> gold_standard = Teuchos::rcp(new Epetra_Vector(ghosted_distributed_parameter->Map()));
      gold_standard->PutScalar(1.0);

      double tol = 10.0 * Teuchos::ScalarTraits<double>::eps();

      TEST_EQUALITY_CONST(testEqualityOfEpetraVectorValues(*ghosted_distributed_parameter,*gold_standard,tol), true);
    }

  }

#ifdef HAVE_STOKHOS
  Teuchos::RCP<Stokhos::SGModelEvaluator> buildStochModel(const Teuchos::RCP<const Epetra_Comm> & Comm,
                                                          const Teuchos::RCP<EpetraExt::ModelEvaluator> & model,
                                                          Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion,
                                                          bool fullExpansion);
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > buildExpansion(int numDim,int order,bool fullExpansion);

  TEUCHOS_UNIT_TEST(model_evaluator, sg)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    bool parameter_on = true;
    AssemblyPieces ap;
    Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> > sg_lof;

    bool fullExpansion = true;
    RCP<Stokhos::OrthogPolyExpansion<int,double> > sgExpansion = buildExpansion(2,4,fullExpansion);
  
    buildAssemblyPieces(parameter_on,ap,sg_lof,sgExpansion);
  
    // Test a transient me, with basic values (no SG)
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = true;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,sg_lof,p_names,ap.gd,build_transient_support));

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_W));

      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_sg));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_sg)); // not currently supported
      
      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      RCP<Epetra_Vector> x_dot = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      x_dot->PutScalar(0.0);
      in_args.set_x(x);
      in_args.set_x_dot(x_dot);
      in_args.set_alpha(0.0);
      in_args.set_beta(1.0);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Operator> J_tmp = me->create_W();
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J_tmp);
      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);

      // compare or save residual vector for comparison with SG version
      // basic_trans_f = out_args.get_f();
      if(basic_trans_f==Teuchos::null)
      {   TEST_ASSERT(false); }
      else {
         double norm=0.0, diff=0.0;
         f->Norm2(&norm);
         f->Update(-1.0,*basic_trans_f,1.0);    
         f->Norm2(&diff);
         TEST_ASSERT(diff/norm <= 1e-15);
      }
    }

    // Test a steady-state me, basic values (no SG)
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = false;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,sg_lof,p_names,ap.gd,build_transient_support));

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_W));

      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_sg));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_sg));
      
      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      in_args.set_x(x);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Operator> J_tmp = me->create_W();
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J_tmp);
      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);

      // compare or save residual vector for comparison with SG version
      if(basic_ss_f==Teuchos::null)
      {   TEST_ASSERT(false); }
      else {
         double norm=0.0, diff=0.0;
         f->Norm2(&norm);
         f->Update(-1.0,*basic_ss_f,1.0);    
         f->Norm2(&diff);
         out << diff << " " << norm << std::endl;
         TEST_ASSERT(diff/norm <= 1e-15);
      }
    }

    // Test a steady-state me, SG values
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = false;
      RCP<panzer::ModelEvaluator_Epetra> pan_me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(ap.fmb,ap.rLibrary,sg_lof,p_names,ap.gd,build_transient_support));
      RCP<EpetraExt::ModelEvaluator> me = buildStochModel(Comm,pan_me,sgExpansion,fullExpansion);

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      RCP<const Epetra_Map> x_map = me->get_x_map();
      RCP<const Epetra_Map> f_map = me->get_f_map();

      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(me->create_W());

      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      in_args.set_x(x);

      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);
    }

  }

  // Hepler functions
  ////////////////////////////////////////////////////////////

  Teuchos::RCP<Stokhos::SGModelEvaluator> buildStochModel(const Teuchos::RCP<const Epetra_Comm> & Comm,
                                                          const Teuchos::RCP<EpetraExt::ModelEvaluator> & model,
                                                          Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion,
                                                          bool fullExpansion)
  {
     Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = expansion->getTripleProduct();
     Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = expansion->getBasis();
  
     // build Parallel data
     Teuchos::ParameterList parallelParams; // empty!
     Teuchos::RCP<Stokhos::ParallelData> sgParallelData 
        = Teuchos::rcp(new Stokhos::ParallelData(basis,Cijk,Comm,parallelParams));
  
     // build parameter list
     Teuchos::RCP<Teuchos::ParameterList> sgParams = Teuchos::rcp(new Teuchos::ParameterList);
     if(!fullExpansion) {
        sgParams->set("Parameter Expansion Type","Linear");
        sgParams->set("Jacobian Expansion Type","Linear");
     }
     Teuchos::ParameterList & sgOpParams = sgParams->sublist("SG Operator");
     sgOpParams.set("Operator Method","Fully Assembled");
  
     Teuchos::ParameterList & sgPrecParams = sgParams->sublist("SG Preconditioner");
     sgPrecParams.set("Preconditioner Method","Mean-based");
     sgPrecParams.set("Mean Preconditioner Type","ML");
     Teuchos::ParameterList & precParams = sgPrecParams.sublist("Mean Preconditioner Parameters");
     precParams.set("default values","SA");
  
     // build model evaluator
     Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model 
        = Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
                                                     expansion, sgParallelData, 
                                                     sgParams));
  
     return sg_model;
  }

  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > buildExpansion(int numDim,int order,bool fullExpansion)
  { 
     Teuchos::Array<Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(numDim);
     for(int i=0;i<numDim;i++)
        bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(order));
     Teuchos::RCP<Stokhos::ProductBasis<int,double> > basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    
     // build Cijk and "expansion"
     Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;
     if(!fullExpansion)
       Cijk = basis->computeLinearTripleProductTensor();
     else
       Cijk = basis->computeTripleProductTensor();
     Teuchos::RCP<Stokhos::Quadrature<int,double> > quadrature = Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    
     return Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis,Cijk,quadrature));
  }

#endif

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
  
  /** Compares Epetra_Vector values and returns true if difference is
      less than tolerance
  */
  bool testEqualityOfEpetraVectorValues(Epetra_Vector& a, Epetra_Vector& b, double tolerance, bool write_to_cout)
  {  
    bool is_equal = true;

    TEUCHOS_ASSERT(a.Map().NumMyElements() == b.Map().NumMyElements());

    for (int i = 0; i < a.Map().NumMyElements(); ++i) {
      
      std::string output = "    equal!: ";
      
      if (std::fabs(a[i] - b[i]) > tolerance) {
	is_equal = false;
	output = "NOT equal!: ";
      }
      
      if (write_to_cout)
	std::cout << output << a[i] << " - " << b[i] << " = " << (a[i] - b[i]) << std::endl;
    }

    int globalSuccess = -1;
    Teuchos::RCP<const Teuchos::Comm<Teuchos::Ordinal> > comm = Teuchos::DefaultComm<Teuchos::Ordinal>::getComm();
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, is_equal ? 0 : 1, Teuchos::outArg(globalSuccess) );
    return (globalSuccess==0);
  }

  void buildAssemblyPieces(bool parameter_on,
                           AssemblyPieces & ap
                           #ifdef HAVE_STOKHOS
                           , Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> > & sg_lof
                           , const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & sgExpansion
                           #endif
                           )
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
    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
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

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk_classic::STKConnManager<int>(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    ap.dofManager = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),ap.physicsBlocks,conn_manager);

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    ap.ep_lof = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm.getConst(),ap.dofManager));
    linObjFactory = ap.ep_lof;
    #ifdef HAVE_STOKHOS
    if(sgExpansion!=Teuchos::null) {
       sg_lof = Teuchos::rcp(new panzer::SGEpetraLinearObjFactory<panzer::Traits,int>(ap.ep_lof,sgExpansion,Teuchos::null));
       linObjFactory = sg_lof;
    }
    #endif

    ap.rLibrary = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,ap.dofManager,linObjFactory)); 

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    // Add in the application specific closure model factory
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    ap.cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    if(parameter_on)
       closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<std::string>("Type","Parameter");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("DENSITY").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("HEAT_CAPACITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    #ifdef HAVE_STOKHOS
    if(sgExpansion!=Teuchos::null)
       closure_models.sublist("ion solid").sublist("ION_DENSITY").set("Expansion",sgExpansion);
    else
       closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    #else
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    #endif
    closure_models.sublist("ion solid").sublist("ION_HEAT_CAPACITY").set<double>("Value",1.0);
    ap.closure_models = closure_models;

    ap.user_data = Teuchos::ParameterList("User Data");

    ap.fmb->setWorksetContainer(wkstContainer);
    ap.fmb->setupVolumeFieldManagers(ap.physicsBlocks,ap.cm_factory,closure_models,*linObjFactory,ap.user_data);
    ap.fmb->setupBCFieldManagers(bcs,ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,bc_factory,closure_models,*linObjFactory,ap.user_data);
  }


}
