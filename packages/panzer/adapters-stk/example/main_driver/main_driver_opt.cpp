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
#include "Teuchos_YamlParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

#include "Panzer_NodeType.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "NOX_Utils.H"
#include "NOX_Observer_Print.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorObserverComposite.hpp"
#include "Tempus_IntegratorObserverBasic.hpp"
#include "Tempus_StepperFactory.hpp"

#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include "user_app_NOXObserverFactory.hpp"
#ifdef PANZER_HAVE_TEMPUS
#include "user_app_TempusObserverFactory.hpp"
#endif
#include "user_app_ResponseEvaluatorFactory_HOFlux.hpp"

#include "Panzer_ROLTempusReducedObjective.hpp"

#include "Panzer_ResponseEvaluatorFactory_Probe.hpp"

#include <Ioss_SerializeIO.h>

#include <string>
#include <iostream>

// *************************************************
// applicaiton functions defined in this file
// *************************************************
namespace user_app {
  void addResponsesToModelEvaluatorFactory(const Teuchos::ParameterList& response_sublist,
                                           panzer_stk::ModelEvaluatorFactory<double>& me_factory);
  void computeAndPrintResponses(const Teuchos::RCP<Thyra::ModelEvaluator<double>>& me,
                                const auto& x,
                                const auto& x_dot,
                                const auto& t,
                                const auto& global_data,
                                std::ostream& out);
}

// *************************************************
// *************************************************
//  Main Driver for a transient optimization problem using ROL and
//  Tempus
// *************************************************
// *************************************************

int main(int argc, char *argv[])
{
  int status = 0;

  std::ostream blackhole{nullptr};
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  Kokkos::initialize(argc,argv);

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  if (mpiSession.getNProc() > 1) {
    out->setShowProcRank(true);
    out->setOutputToRootOnly(0);
  }

  try {
    const auto stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Panzer Main Driver: Transient Optimization"));
    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // Parse the command line arguments
    std::string input_file_name = "energy-transient-tempus-opt-blocked.xml";
    int exodus_io_num_procs = 0;
    bool printTimers = true;
    bool printInputPL = false;
    bool overrideNoxOutput = true;
    {
      Teuchos::CommandLineProcessor clp;

      clp.setOption("i", &input_file_name, "Input file name, supports xml or yaml");
      clp.setOption("exodus-io-num-procs", &exodus_io_num_procs, "Number of processes that can access the file system at the same time to read their portion of a sliced exodus file in parallel.  Defaults to 0 - implies all processes for the run can access the file system at the same time.");
      clp.setOption("time","no-time", &printTimers, "Print the timing information.");
      clp.setOption("pl","no-pl", &printInputPL, "Print the input ParameterList at the start of the run.");
      clp.setOption("override-nox-output","default-nox-output", &overrideNoxOutput, "Override default nox printing with new print observer.");

      auto parse_return = clp.parse(argc,argv,&std::cerr);

      TEUCHOS_TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL,
                                 std::runtime_error, "Failed to parse command line!");
    }

    TEUCHOS_ASSERT(exodus_io_num_procs <= comm->getSize());
    if (exodus_io_num_procs != 0)
      Ioss::SerializeIO::setGroupFactor(exodus_io_num_procs);

    // Parse the input file and broadcast to other processes
    Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("User_App Parameters"));
    {
      const std::string check_yaml = ".yaml";
      const auto search_yaml = input_file_name.find(check_yaml);
      if (search_yaml != std::string::npos) {
        Teuchos::updateParametersFromYamlFileAndBroadcast(input_file_name, input_params.ptr(), *comm);
      } else {
        Teuchos::updateParametersFromXmlFileAndBroadcast(input_file_name, input_params.ptr(), *comm);
      }
    }

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
    input_params->sublist("User Data").set("Comm", comm); // Required for GlobalStatistics

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
       = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>());

    Teuchos::RCP<Thyra::ModelEvaluator<double> > physics;
    Teuchos::RCP<Thyra::ModelEvaluator<double> > solver;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    Teuchos::RCP<panzer::GlobalIndexer> globalIndexer;
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    {
      Teuchos::ParameterList responses = input_params->sublist("Responses");
      input_params->remove("Responses");

      panzer_stk::ModelEvaluatorFactory<double> me_factory;
      me_factory.setParameterList(input_params);
      auto strat_list = Teuchos::parameterList();
      *strat_list = input_params->sublist("Solution Control",true).sublist("Tempus",true).sublist("My Example Stepper",true).sublist("My Example Solver",true).sublist("NOX",true).sublist("Direction",true).sublist("Newton",true).sublist("Stratimikos Linear Solver",true).sublist("Stratimikos",true);
      me_factory.setStratimikosList(strat_list);
      me_factory.buildObjects(comm,global_data,eqset_factory,bc_factory,cm_factory);

      user_app::addResponsesToModelEvaluatorFactory(responses,me_factory);
      me_factory.buildResponses(cm_factory);

      physicsBlocks = me_factory.getPhysicsBlocks();
      physics = me_factory.getPhysicsModelEvaluator();
      rLibrary = me_factory.getResponseLibrary();
      linObjFactory = me_factory.getLinearObjFactory();
      globalIndexer = me_factory.getGlobalIndexer();
      mesh = me_factory.getMesh();
    }

    // setup outputs to mesh on the stkIOResponseLibrary
    {
      stkIOResponseLibrary->initialize(*rLibrary);

      Teuchos::ParameterList user_data(input_params->sublist("User Data"));
      user_data.set<int>("Workset Size",input_params->sublist("Assembly").get<int>("Workset Size"));

      std::vector<std::string> eBlocks;
      mesh->getElementBlockNames(eBlocks);

      panzer_stk::RespFactorySolnWriter_Builder builder;
      builder.mesh = mesh;
      stkIOResponseLibrary->addResponse("Main Field Output",eBlocks,builder);

      stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                                    cm_factory,
                                                    input_params->sublist("Closure Models"),
                                                    user_data);
    }

    // Tempus stepper
    using namespace Teuchos;
    using namespace Tempus;

    RCP<ParameterList> tempus_pl = parameterList(std::string("Tempus"));
    *tempus_pl = input_params->sublist("Solution Control").sublist("Tempus");

    if (overrideNoxOutput) {
      auto nox_utils = Teuchos::rcp(new NOX::Utils(NOX::Utils::Error,comm->getRank(),0,3));
      auto print_nonlinear = Teuchos::rcp(new NOX::ObserverPrint(nox_utils,10));
      tempus_pl->sublist("My Example Stepper",true).sublist("My Example Solver",true).sublist("NOX",true).sublist("Solver Options",true).set<RCP<NOX::Observer>>("Observer",print_nonlinear);
    }

    // Tempus is overriding NOX parameters with its own defaults. Need
    // to fix this. It also does not validate because of arbitrary
    // solver name. Should validate but use
    // disableRecursiveValidation() to avoid errors with arbitrary
    // names.

    const bool doInitialization = false;
    auto integrator = createIntegratorBasic<double>(tempus_pl, physics, doInitialization);

    RCP<ParameterList> noxList = parameterList("Correct NOX Params");
    *noxList = tempus_pl->sublist("My Example Stepper",true).sublist("My Example Solver",true).sublist("NOX",true);
    // noxList->print(std::cout);
    integrator->getStepper()->getSolver()->setParameterList(noxList);
    integrator->initialize();

    // Setting observers on tempus breaks the screen output of the
    // time steps! It replaces IntegratorObserverBasic which handles
    // IO. Is this at least documented?
    {
      RCP<Tempus::IntegratorObserverComposite<double>> tempus_observers = rcp(new Tempus::IntegratorObserverComposite<double>);
      RCP<const panzer_stk::TempusObserverFactory> tof = Teuchos::rcp(new user_app::TempusObserverFactory(stkIOResponseLibrary,rLibrary->getWorksetContainer()));
      tempus_observers->addObserver(Teuchos::rcp(new Tempus::IntegratorObserverBasic<double>));
      tempus_observers->addObserver(tof->buildTempusObserver(mesh,globalIndexer,linObjFactory));
      integrator->setObserver(tempus_observers);
    }

    RCP<Thyra::VectorBase<double>> x0 = physics->getNominalValues().get_x()->clone_v();
    integrator->initializeSolutionHistory(0.0, x0);

    bool integratorStatus = integrator->advanceTime();
    TEUCHOS_ASSERT(integratorStatus);

    auto x = integrator->getSolutionHistory()->getCurrentState()->getX();
    auto x_dot = integrator->getSolutionHistory()->getCurrentState()->getXDot();
    auto t = integrator->getSolutionHistory()->getCurrentState()->getTime();

    user_app::computeAndPrintResponses(physics,x,x_dot,t,global_data,*out);

    {
      Teuchos::RCP<Teuchos::ParameterList> rol_params = Teuchos::rcp(new Teuchos::ParameterList("ROL Parameters"));
      {
        const std::string rol_input_file = "rol_params.xml";
        Teuchos::updateParametersFromXmlFileAndBroadcast(rol_input_file, rol_params.ptr(), *comm);

        RCP<ParameterList> pl = sublist(rol_params, "Tempus", true);
        RCP<ParameterList> objective_params = sublist(rol_params, "Objective", true);
        ROL::Ptr<ROL::TempusReducedObjective<double> > objective =
          ROL::makePtr<ROL::TempusReducedObjective<double> >(physics, pl, objective_params);

        // Create target -- do forward integration with perturbed parameter values
        // RCP<ROL::Vector<double> > zs = objective->create_design_vector();
        // zs->setScalar(1.5);
        // RCP<ROL::Vector<double> > r = objective->create_response_vector();
        // objective->run_tempus(*r, *zs);
        // objective->set_target(r);

      }
    }

    stackedTimer->stopBaseTimer();
    if (printTimers) {
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = true;
      options.output_minmax = false;
      options.output_histogram = false;
      options.num_histogram = 5;
      std::string timing_file_name = input_file_name+"_timing.log";
      std::fstream timing_file{timing_file_name,std::ios::out | std::ios::trunc};
      stackedTimer->report(timing_file, Teuchos::DefaultComm<int>::getComm(), options);
    }

  }
  catch (std::exception& e) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << e.what() << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }

  if (status == 0)
    *out << "panzer::MainDriver run completed." << std::endl;

  return status;
}

// *************************************************
// *************************************************

void user_app::addResponsesToModelEvaluatorFactory(const Teuchos::ParameterList& responses,
                                                   panzer_stk::ModelEvaluatorFactory<double>& me_factory)
{
  for(Teuchos::ParameterList::ConstIterator itr=responses.begin();itr!=responses.end();++itr) {
    const std::string name = responses.name(itr);
    TEUCHOS_ASSERT(responses.entry(itr).isList());
    Teuchos::ParameterList & lst = Teuchos::getValue<Teuchos::ParameterList>(responses.entry(itr));

    if (lst.get<std::string>("Type") == "Functional") {

      panzer::FunctionalResponse_Builder<int,int> builder;
      builder.comm = MPI_COMM_WORLD;
      builder.cubatureDegree = 2;
      builder.requiresCellIntegral = true;
      builder.quadPointField = lst.get<std::string>("Field Name");

      std::vector<std::string> eblocks;
      panzer::StringTokenizer(eblocks,lst.get<std::string>("Element Blocks"),",",true);

      std::vector<panzer::WorksetDescriptor> wkst_descs;
      for(std::size_t i=0;i<eblocks.size();i++)
        wkst_descs.push_back(panzer::blockDescriptor(eblocks[i]));

      me_factory.addResponse(name,wkst_descs,builder);
    }
    else if (lst.get<std::string>("Type") == "Point Value") {
      panzer::ProbeResponse_Builder<panzer::LocalOrdinal,panzer::GlobalOrdinal> builder;
      builder.comm = MPI_COMM_WORLD;
      builder.point = Teuchos::Array<double>{0.5,0.5};
      builder.cubatureDegree = 2;
      builder.fieldName = lst.get<std::string>("Field Name");
      builder.applyDirichletToDerivative = false;

      std::vector<std::string> eblocks;
      panzer::StringTokenizer(eblocks,lst.get<std::string>("Element Blocks"),",",true);

      std::vector<panzer::WorksetDescriptor> descriptors;
      for(std::size_t i=0;i<eblocks.size();i++)
        descriptors.push_back(panzer::blockDescriptor(eblocks[i]));

      me_factory.addResponse("Value In Middle",descriptors,builder);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Error: Response type of \"" << lst.get<std::string>("Type") << "\" is not supported!");
    }
  }
}

void user_app::computeAndPrintResponses(const Teuchos::RCP<Thyra::ModelEvaluator<double>>& physics,
                                        const auto& x,
                                        const auto& x_dot,
                                        const auto& t,
                                        const auto& global_data,
                                        std::ostream& os)
{
  if(physics->Ng()>0) {

    os << "Ng = " << physics->Ng() << std::endl;
    os << "Np = " << physics->Np() << std::endl;
    Thyra::ModelEvaluatorBase::InArgs<double> respInArgs = physics->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<double> respOutArgs = physics->createOutArgs();

    TEUCHOS_ASSERT(physics->Ng()==respOutArgs.Ng());

    respInArgs.set_x(x);
    respInArgs.set_x_dot(x_dot);

    for(int g=0;g<respOutArgs.Ng();++g) {
      Teuchos::RCP<Thyra::VectorBase<double> > response = Thyra::createMember(*physics->get_g_space(g));
      respOutArgs.set_g(g,response);

      for (int p = 0; p < physics->Np(); ++p) {
        bool derivative_supported = respOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,g,p).supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        os << "DgDp(" << g << "," << p << ") supports = " << derivative_supported << std::endl;
        if (derivative_supported) {
          auto dgdp = physics->create_DgDp_op(g,p);
          respOutArgs.set_DgDp(g,p,Thyra::ModelEvaluatorBase::Derivative<double>(dgdp));
        }
      }
    }

    physics->evalModel(respInArgs, respOutArgs);

    {
      auto parameter_library = global_data->pl;
      TEUCHOS_ASSERT(nonnull(parameter_library));
      for (auto i=parameter_library->begin();i != parameter_library->end(); ++i) {
        os << "Sacado::ParameterLibrary: " << i->first << std::endl;
      }
    }

    {
      auto nominalValues = physics->getNominalValues();
      for (int i=0; i < respInArgs.Np();++i) {
        auto p = nominalValues.get_p(i);
        auto p_names = physics->get_p_names(i);
        os << "ModelEvaluator::NominalValues Parameter Value: \"" << (*p_names)[0] << "\" = " << Thyra::get_ele(*p,0) << std::endl;
      }
    }

    for (int i=0;i<respOutArgs.Ng();i++) {
      Teuchos::RCP<Thyra::VectorBase<double> > response = respOutArgs.get_g(i);
      TEUCHOS_ASSERT(response!=Teuchos::null);
      os << "Response Value: \"" << physics->get_g_names(i)[0] << "\" = " << Thyra::get_ele(*response,0) << std::endl;
      for (int j=0; j < respOutArgs.Np(); ++j) {
        // os << "  dg(" << i << ")/dp(" << j << ") supports(GRAD_FORM) = " << respOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,j).supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) << std::endl;
        // os << "  dg(" << i << ")/dp(" << j << ") supports(JAC_FORM)  = " << respOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,j).supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) << std::endl;
      }
    }
  }
}
