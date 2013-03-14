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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

#include "Kokkos_DefaultNode.hpp"

#include "Panzer_ConfigDefs.hpp"
#include "Panzer_STK_ModelEvaluatorFactory_Epetra.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_PauseToAttach.hpp"

#ifdef Panzer_BUILD_PAPI_SUPPORT
#include "Panzer_PAPI_Counter.hpp"
#endif

#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include "user_app_NOXObserverFactory_Epetra.hpp"
#include "user_app_RythmosObserverFactory_Epetra.hpp"
#include "user_app_ResponseAggregatorFactory.hpp"
#include "user_app_ResponseEvaluatorFactory_HOFlux.hpp"

#include <Ioss_SerializeIO.h>

#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
  int status = 0;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  if (mpiSession.getNProc() > 1) {
    out->setShowProcRank(true);
    out->setOutputToRootOnly(0);
  }

#ifdef Panzer_BUILD_PAPI_SUPPORT
    panzer::PAPICounter papi_counter("Panzer: Total Execution", mpiSession.getRank(), MPI_COMM_WORLD);
    papi_counter.start();
#endif

  try {
    
    Teuchos::RCP<Teuchos::Time> total_time = 
      Teuchos::TimeMonitor::getNewTimer("User App: Total Time");
    
    Teuchos::TimeMonitor timer(*total_time); 
    
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    
    // Parse the command line arguments
    std::string input_file_name = "user_app.xml";
    int exodus_io_num_procs = 0;
    bool pauseToAttachOn = false;
    bool fluxCalculation = false;
    {
      Teuchos::CommandLineProcessor clp;
      
      clp.setOption("i", &input_file_name, "User_App input xml filename");
      clp.setOption("exodus-io-num-procs", &exodus_io_num_procs, "Number of processes that can access the file system at the same time to read their portion of a sliced exodus file in parallel.  Defaults to 0 - implies all processes for the run can access the file system at the same time.");
      clp.setOption("pause-to-attach","disable-pause-to-attach", &pauseToAttachOn, "Call pause to attach, default is off.");
      clp.setOption("flux-calc","disable-flux-calc", &fluxCalculation, "Enable the flux calculation.");
      
      Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = 
	clp.parse(argc,argv,&std::cerr);
      
      TEUCHOS_TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL, 
			 std::runtime_error, "Failed to parse command line!");
    }

    if(pauseToAttachOn)
       panzer::pauseToAttach();

    TEUCHOS_ASSERT(exodus_io_num_procs <= comm->getSize());
    if (exodus_io_num_procs != 0)
      Ioss::SerializeIO::setGroupFactor(exodus_io_num_procs);

    // Parse the input file and broadcast to other processes
    Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("User_App Parameters"));
    Teuchos::updateParametersFromXmlFileAndBroadcast(input_file_name, input_params.ptr(), *comm);
    
    *out << *input_params << std::endl;
    
    // Add in the application specific equation set factory
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);

    // Add in the application specific closure model factory
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;  
    cm_factory.buildObjects(cm_builder);

    // Add in the application specific bc factory
    user_app::BCFactory bc_factory; 

    // Create the global data
    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();

    // A GlobalStatistics closure model requires the comm to be set in the user data.
    input_params->sublist("User Data").set("Comm", comm);

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
       = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>());

    // Add in the application specific observer factories
    {
      // Rythmos
      {
	Teuchos::RCP<const panzer_stk::RythmosObserverFactory> rof = 
	  Teuchos::rcp(new user_app::RythmosObserverFactory_Epetra(stkIOResponseLibrary));
	input_params->sublist("Solver Factories").set("Rythmos Observer Factory", rof);
      }

      // NOX
      {
	Teuchos::RCP<user_app::NOXObserverFactory_Epetra> nof = 
	  Teuchos::rcp(new user_app::NOXObserverFactory_Epetra(stkIOResponseLibrary));

	Teuchos::RCP<Teuchos::ParameterList> observers_to_build = 
	  Teuchos::parameterList(input_params->sublist("Solver Factories").sublist("NOX Observers"));
	  
	input_params->sublist("Solver Factories").remove("NOX Observers");

	nof->setParameterList(observers_to_build);

	input_params->sublist("Solver Factories").set<Teuchos::RCP<const panzer_stk::NOXObserverFactory> >("NOX Observer Factory", nof);
      }
    }

    Teuchos::RCP<Thyra::ModelEvaluator<double> > physics;
    Teuchos::RCP<Thyra::ModelEvaluator<double> > solver;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    {
      panzer_stk::ModelEvaluatorFactory_Epetra<double> me_factory;
      Teuchos::RCP<user_app::MyResponseAggregatorFactory<panzer::Traits> > ra_factory = 
	Teuchos::rcp(new user_app::MyResponseAggregatorFactory<panzer::Traits>);
      me_factory.setParameterList(input_params);
      me_factory.buildObjects(comm,global_data,eqset_factory,bc_factory,cm_factory,ra_factory.ptr()); 

      physics = me_factory.getPhysicsModelEvaluator();
      solver = me_factory.getResponseOnlyModelEvaluator();
      rLibrary = me_factory.getResponseLibrary();
      physicsBlocks = me_factory.getPhysicsBlocks();
      linObjFactory = me_factory.getLinearObjFactory();
    }  

    // setup outputs to mesh on the stkIOResponseLibrary
    ////////////////////////////////////////////////////////////////

    stkIOResponseLibrary->initialize(*rLibrary);
    
    // 1. Register correct aggregator and reserve response - this was done in the appropriate observer object

    // 2. Build volume field managers
    {
      Teuchos::ParameterList user_data(input_params->sublist("User Data"));
      user_data.set<int>("Workset Size",input_params->sublist("Assembly").get<unsigned long>("Workset Size"));
    
      stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                        cm_factory,
                                        input_params->sublist("Closure Models"),
                                        user_data);
    }
    
    // setup outputs to mesh on the fluxResponseLibrary
    ////////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > fluxResponseLibrary 
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>);
    
    if(fluxCalculation) {
      fluxResponseLibrary->initialize(*rLibrary);
  
      // build high-order flux response
      {
        user_app::HOFluxResponse_Builder builder;
        builder.comm = MPI_COMM_WORLD;
        builder.cubatureDegree = 2;
  
        std::vector<panzer::WorksetDescriptor> sidesets;
        sidesets.push_back(panzer::sidesetVolumeDescriptor("eblock-0_0","left"));
  
        fluxResponseLibrary->addResponse("HO-Flux",sidesets,builder);
      }
  
      {
        Teuchos::ParameterList user_data(input_params->sublist("User Data"));
        user_data.set<int>("Workset Size",input_params->sublist("Assembly").get<unsigned long>("Workset Size"));
      
        fluxResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                          *eqset_factory,
                                          cm_factory,
                                          input_params->sublist("Closure Models"),
                                          user_data);
      }
  
      {
        Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Residual> > resp
            = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Residual> >(fluxResponseLibrary->getResponse<panzer::Traits::Residual>("HO-Flux"),true);
  
        // allocate a vector
        Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*resp->getMap()));
        resp->setVector(vec);
      }
    }
    
    ////////////////////////////////////////////////////////////////

    // solve the system
    {
      
      // Set inputs
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = solver->createInArgs();
      const Thyra::ModelEvaluatorBase::InArgs<double> inArgsNominal = solver->getNominalValues();

      // Set outputs
      Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = solver->createOutArgs();

      // Solution vector is returned as extra respons vector
      Teuchos::RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*physics->get_x_space());
      for(std::size_t i=0;i<rLibrary->getLabeledResponseCount();i++)
         outArgs.set_g(i,Teuchos::null);
      outArgs.set_g(rLibrary->getLabeledResponseCount(),gx);

      // Now, solve the problem and return the responses
      solver->evalModel(inArgs, outArgs);

      // number of responses, minus the solution vector
      TEUCHOS_ASSERT(rLibrary->getLabeledResponseCount()==Teuchos::as<std::size_t>(outArgs.Ng()-1));
      
      // get responses if there are any
      //////////////////////////////////////////////
      if(rLibrary->getLabeledResponseCount()>0) {

         std::vector<std::string> labels;
         rLibrary->getVolumeResponseLabels(labels);

         Thyra::ModelEvaluatorBase::InArgs<double> respInArgs = physics->createInArgs();
         Thyra::ModelEvaluatorBase::OutArgs<double> respOutArgs = physics->createOutArgs();

         TEUCHOS_ASSERT(rLibrary->getLabeledResponseCount()==Teuchos::as<std::size_t>(respOutArgs.Ng()));
   
         respInArgs.set_x(gx);
   
         // set up response out args
         for(int i=0;i<respOutArgs.Ng();i++) {
	   Teuchos::RCP<Thyra::VectorBase<double> > response = Thyra::createMember(*physics->get_g_space(i));
           respOutArgs.set_g(i,response);
         }
   
         rLibrary->printVolumeContainers(*out);
   
         // Now, solve the problem and return the responses
         physics->evalModel(respInArgs, respOutArgs);
   
         // loop over out args for printing
         for(int i=0;i<respOutArgs.Ng();i++) {
	   Teuchos::RCP<Thyra::VectorBase<double> > response = respOutArgs.get_g(i);

            TEUCHOS_ASSERT(response!=Teuchos::null); // should not be null!

            *out << "Response Value \"" << labels[i] << "\": " << *response << std::endl;
         }
      }

      if(fluxCalculation) {
        // initialize the assembly container
        panzer::AssemblyEngineInArgs ae_inargs;
        ae_inargs.container_ = linObjFactory->buildLinearObjContainer();
        ae_inargs.ghostedContainer_ = linObjFactory->buildGhostedLinearObjContainer();
        ae_inargs.alpha = 0.0;
        ae_inargs.beta = 1.0;
        ae_inargs.evaluate_transient_terms = false;

        // initialize the ghosted container
        linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

        const Teuchos::RCP<panzer::EpetraLinearObjContainer> epGlobalContainer
           = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs.container_,true);
        epGlobalContainer->set_x_th(gx);

        // evaluate current on contacts
        fluxResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
        fluxResponseLibrary->evaluate<panzer::Traits::Residual>(ae_inargs);

        // output current values
        *out << "\nFlux values: \n";
        {
          std::string currentRespName = "HO-Flux";

          Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > resp
              = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(fluxResponseLibrary->getResponse<panzer::Traits::Residual>(currentRespName),true);

          *out << "   " << currentRespName << " = " << resp->value << std::endl;
        }
      }
    }

    // Clear out a circular dependency on rcp's
    stkIOResponseLibrary->free();
    fluxResponseLibrary->free();

  }
  catch (std::exception& e) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << e.what() << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }
  catch (std::string& msg) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << msg << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }
  catch (...) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << "Caught UNKOWN exception" << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }
  
  Teuchos::TimeMonitor::summarize(*out,false,true,false);

#ifdef Panzer_BUILD_PAPI_SUPPORT
    papi_counter.stop();
    papi_counter.report(std::cout);
#endif

  if (status == 0)
    *out << "panzer::MainDriver run completed." << std::endl;

  return status;
}
