// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

#include "Panzer_NodeType.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#ifdef Panzer_BUILD_PAPI_SUPPORT
#include "Panzer_PAPI_Counter.hpp"
#endif

#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include "user_app_NOXObserverFactory.hpp"
#ifdef PANZER_HAVE_TEMPUS
#include "user_app_TempusObserverFactory.hpp"
#endif
#include "user_app_ResponseEvaluatorFactory_HOFlux.hpp"

#include "Panzer_ResponseEvaluatorFactory_Probe.hpp"

#include <Ioss_SerializeIO.h>

#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
  int status = 0;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  Kokkos::initialize(argc,argv);

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
     const auto stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Panzer Main Driver"));
     Teuchos::TimeMonitor::setStackedTimer(stackedTimer);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // Parse the command line arguments
    std::string input_file_name = "user_app.xml";
    int exodus_io_num_procs = 0;
    bool pauseToAttachOn = false;
    bool fluxCalculation = false;
    bool pointCalculation = false;
    bool printTimers = false;
    bool printInputPL = false;
    {
      Teuchos::CommandLineProcessor clp;

      clp.setOption("i", &input_file_name, "User_App input xml filename");
      clp.setOption("exodus-io-num-procs", &exodus_io_num_procs, "Number of processes that can access the file system at the same time to read their portion of a sliced exodus file in parallel.  Defaults to 0 - implies all processes for the run can access the file system at the same time.");
      clp.setOption("pause-to-attach","disable-pause-to-attach", &pauseToAttachOn, "Call pause to attach, default is off.");
      clp.setOption("flux-calc","disable-flux-calc", &fluxCalculation, "Enable the flux calculation.");
      clp.setOption("point-calc","disable-point-calc", &pointCalculation, "Enable the probe evaluator unit test.");
      clp.setOption("time","no-time", &printTimers, "Print the timing information.");
      clp.setOption("pl","no-pl", &printTimers, "Print the input ParameterList at the start of the run.");

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

    if (printInputPL)
      *out << *input_params << std::endl;

    Teuchos::ParameterList solver_factories = input_params->sublist("Solver Factories");
    input_params->remove("Solver Factories");

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

    // extract and then remove the volume responses
    Teuchos::ParameterList responses = input_params->sublist("Responses");
    input_params->remove("Responses");

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
       = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>());

    Teuchos::RCP<Thyra::ModelEvaluator<double> > physics;
    Teuchos::RCP<Thyra::ModelEvaluator<double> > solver;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    std::map<int,std::string> responseIndexToName;
    {
      panzer_stk::ModelEvaluatorFactory<double> me_factory;

      me_factory.setParameterList(input_params);
      me_factory.buildObjects(comm,global_data,eqset_factory,bc_factory,cm_factory);

      // add a volume response functional for each field
      for(Teuchos::ParameterList::ConstIterator itr=responses.begin();itr!=responses.end();++itr) {
        const std::string name = responses.name(itr);
        TEUCHOS_ASSERT(responses.entry(itr).isList());
        Teuchos::ParameterList & lst = Teuchos::getValue<Teuchos::ParameterList>(responses.entry(itr));

        // parameterize the builder
        panzer::FunctionalResponse_Builder<int,int> builder;
        builder.comm = MPI_COMM_WORLD; // good enough
        builder.cubatureDegree = 2;
        builder.requiresCellIntegral = lst.isType<bool>("Requires Cell Integral") ? lst.get<bool>("Requires Cell Integral"): false;
        builder.quadPointField = lst.get<std::string>("Field Name");

        // add the response
        std::vector<std::string> eblocks;
        panzer::StringTokenizer(eblocks,lst.get<std::string>("Element Blocks"),",",true);

        std::vector<panzer::WorksetDescriptor> wkst_descs;
        for(std::size_t i=0;i<eblocks.size();i++)
          wkst_descs.push_back(panzer::blockDescriptor(eblocks[i]));

        int respIndex = me_factory.addResponse(name,wkst_descs,builder);
        responseIndexToName[respIndex] = name;
      }

      // enusre all the responses are built
      me_factory.buildResponses(cm_factory);

      physicsBlocks = me_factory.getPhysicsBlocks();
      physics = me_factory.getPhysicsModelEvaluator();
      rLibrary = me_factory.getResponseLibrary();
      linObjFactory = me_factory.getLinearObjFactory();

      // Add in the application specific observer factories
      {
	      // NOX
        Teuchos::RCP<user_app::NOXObserverFactory> nof;
        {
                nof = Teuchos::rcp(new user_app::NOXObserverFactory(stkIOResponseLibrary));

          Teuchos::RCP<Teuchos::ParameterList> observers_to_build =
            Teuchos::parameterList(solver_factories.sublist("NOX Observers"));

          nof->setParameterList(observers_to_build);
        }

#ifdef PANZER_HAVE_TEMPUS
        Teuchos::RCP<const panzer_stk::TempusObserverFactory> tof;
	      {
          tof = Teuchos::rcp(new user_app::TempusObserverFactory(stkIOResponseLibrary,rLibrary->getWorksetContainer()));
	      }
        solver = me_factory.buildResponseOnlyModelEvaluator(physics,global_data,Teuchos::null,nof.ptr(),tof.ptr());
#else
        solver = me_factory.buildResponseOnlyModelEvaluator(physics,global_data,nof.ptr());
#endif
      }
    }

    // setup outputs to mesh on the stkIOResponseLibrary
    ////////////////////////////////////////////////////////////////

    stkIOResponseLibrary->initialize(*rLibrary);

    // 1. Register correct aggregator and reserve response - this was done in the appropriate observer object

    // 2. Build volume field managers
    {
      Teuchos::ParameterList user_data(input_params->sublist("User Data"));
      user_data.set<int>("Workset Size",input_params->sublist("Assembly").get<int>("Workset Size"));

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
        user_data.set<int>("Workset Size",input_params->sublist("Assembly").get<int>("Workset Size"));

        fluxResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                                     *eqset_factory,
                                                     cm_factory,
                                                     input_params->sublist("Closure Models"),
                                                     user_data);
      }

      {
        Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Residual> > resp
            = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Residual> >(fluxResponseLibrary->getResponse<panzer::Traits::Residual>("HO-Flux"),true);

        const auto vec = Thyra::createMember(*resp->getVectorSpace(),"HO-Flux Vector");
        resp->setVector(vec);
      }
    }

    // setup outputs for the point calculation
    ////////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > pointResponseLibrary
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>);

    if(pointCalculation) {
      pointResponseLibrary->initialize(*rLibrary);

      {
        panzer::ProbeResponse_Builder<panzer::LocalOrdinal,panzer::GlobalOrdinal> builder;
        builder.comm = MPI_COMM_WORLD;
        builder.point = Teuchos::Array<double>{0.5,0.5}; // Bottom
        builder.cubatureDegree = 2;
        builder.fieldName = "TEMPERATURE";
        builder.applyDirichletToDerivative = false;

        std::vector<panzer::WorksetDescriptor> descriptors;
        descriptors.push_back(panzer::WorksetDescriptor("eblock-0_0"));

        pointResponseLibrary->addResponse("Value In Middle",descriptors,builder);
      }

      {
        Teuchos::ParameterList user_data(input_params->sublist("User Data"));
        user_data.set<int>("Workset Size",input_params->sublist("Assembly").get<int>("Workset Size"));

        pointResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                                      *eqset_factory,
                                                      cm_factory,
                                                      input_params->sublist("Closure Models"),
                                                      user_data);
      }

      {
        Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Residual> > resp
            = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Residual> >(pointResponseLibrary->getResponse<panzer::Traits::Residual>("Value In Middle"),true);

        const auto vec = Thyra::createMember(*resp->getVectorSpace(),"Value In Middle Response Thyra Vector");
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

      // Solution vector is returned as extra response vector
      Teuchos::RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*physics->get_x_space());
      for(int i=0;i<outArgs.Ng()-1;i++)
         outArgs.set_g(i,Teuchos::null);
      outArgs.set_g(outArgs.Ng()-1,gx);

      // Now, solve the problem and return the responses
      solver->evalModel(inArgs, outArgs);

      // get responses if there are any
      //////////////////////////////////////////////
      if(physics->Ng()>0) {

         Thyra::ModelEvaluatorBase::InArgs<double> respInArgs = physics->createInArgs();
         Thyra::ModelEvaluatorBase::OutArgs<double> respOutArgs = physics->createOutArgs();

         TEUCHOS_ASSERT(physics->Ng()==respOutArgs.Ng());

         respInArgs.set_x(gx);

         // set up response out args
         for(int i=0;i<respOutArgs.Ng();i++) {
	   Teuchos::RCP<Thyra::VectorBase<double> > response = Thyra::createMember(*physics->get_g_space(i));
           respOutArgs.set_g(i,response);
         }

         // Now, solve the problem and return the responses
         physics->evalModel(respInArgs, respOutArgs);

         // loop over out args for printing
         for(int i=0;i<respOutArgs.Ng();i++) {
	   Teuchos::RCP<Thyra::VectorBase<double> > response = respOutArgs.get_g(i);

            TEUCHOS_ASSERT(response!=Teuchos::null); // should not be null!

            *out << "Response Value \"" << responseIndexToName[i] << "\": " << Thyra::get_ele(*response,0) << std::endl;
         }
      }

      if(fluxCalculation) {
        stackedTimer->start("Flux Response Calculation");
        // initialize the assembly container
        panzer::AssemblyEngineInArgs ae_inargs;
        ae_inargs.container_ = linObjFactory->buildLinearObjContainer();
        ae_inargs.ghostedContainer_ = linObjFactory->buildGhostedLinearObjContainer();
        ae_inargs.alpha = 0.0;
        ae_inargs.beta = 1.0;
        ae_inargs.evaluate_transient_terms = false;

        // initialize the ghosted container
        linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

        const Teuchos::RCP<panzer::ThyraObjContainer<double>> thGlobalContainer
          = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double>>(ae_inargs.container_,true);
        thGlobalContainer->set_x_th(gx);

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
        stackedTimer->stop("Flux Response Calculation");
      }

      if(pointCalculation) {
        stackedTimer->start("Point Value Response Calculation");
        // initialize the assembly container
        panzer::AssemblyEngineInArgs ae_inargs;
        ae_inargs.container_ = linObjFactory->buildLinearObjContainer();
        ae_inargs.ghostedContainer_ = linObjFactory->buildGhostedLinearObjContainer();
        ae_inargs.alpha = 0.0;
        ae_inargs.beta = 1.0;
        ae_inargs.evaluate_transient_terms = false;

        // initialize the ghosted container
        linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

        const Teuchos::RCP<panzer::ThyraObjContainer<double>> thGlobalContainer
          = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double>>(ae_inargs.container_,true);
        thGlobalContainer->set_x_th(gx);

        // evaluate current on contacts
        pointResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
        pointResponseLibrary->evaluate<panzer::Traits::Residual>(ae_inargs);

        // output current values
        *out << "\nPoint Values: \n";
        {
          std::string currentRespName = "Value In Middle";

          Teuchos::RCP<panzer::Response_Probe<panzer::Traits::Residual> > resp
              = Teuchos::rcp_dynamic_cast<panzer::Response_Probe<panzer::Traits::Residual> >(pointResponseLibrary->getResponse<panzer::Traits::Residual>(currentRespName),true);

          // Linear problem with analytic solution
          const double gold_value = 0.5;
          const double tol = 1.0e-8;
          *out << "   " << currentRespName << " = " << resp->value << ", error = " << fabs(resp->value - gold_value) << ", tol = " << tol << std::endl;
          TEUCHOS_ASSERT(fabs(resp->value - gold_value) < tol);
        }

        stackedTimer->stop("Point Value Response Calculation");
      }
    }

    stackedTimer->stopBaseTimer();
    if (printTimers) {
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = true;
      options.output_minmax = true;
      options.output_histogram = true;
      options.num_histogram = 5;
      stackedTimer->report(std::cout, Teuchos::DefaultComm<int>::getComm(), options);
    }

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

  // Teuchos::TimeMonitor::summarize(*out,false,true,false);

#ifdef Panzer_BUILD_PAPI_SUPPORT
    papi_counter.stop();
    papi_counter.report(std::cout);
#endif

  if (status == 0)
    *out << "panzer::MainDriver run completed." << std::endl;

  return status;
}
