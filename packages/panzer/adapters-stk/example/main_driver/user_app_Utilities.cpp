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
#include "user_app_TransientObjective.hpp"

#include "Panzer_ResponseEvaluatorFactory_Probe.hpp"

#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_RandomVector.hpp"

#include <Ioss_SerializeIO.h>

#include <string>
#include <iostream>
#include <tuple>

#include "user_app_Utilities.hpp"

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

std::tuple< Teuchos::RCP<Thyra::ModelEvaluator<double>>,
            Teuchos::RCP<panzer::GlobalData>,
            Teuchos::RCP<panzer_stk::STK_Interface>,
            Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>,
            Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>,
            Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>>,
            Teuchos::RCP<panzer::GlobalIndexer>
            >
user_app::buildModelEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& input_params,
                              const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
{
  using namespace Teuchos;
  using namespace Tempus;

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

    Teuchos::RCP<Thyra::ModelEvaluator<double>> physics;
    Teuchos::RCP<Thyra::ModelEvaluator<double>> solver;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> rLibrary;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock>> physicsBlocks;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>> linObjFactory;
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

    return {physics,global_data,mesh,rLibrary,stkIOResponseLibrary,linObjFactory,globalIndexer};
}

Teuchos::RCP<Tempus::IntegratorBasic<double>>
user_app::buildTimeIntegrator(const Teuchos::RCP<Teuchos::ParameterList>& input_params,
                              const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              Teuchos::RCP<Thyra::ModelEvaluator<double>> physics,
                              Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                              Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> rLibrary,
                              Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> stkIOResponseLibrary,
                              Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>> linObjFactory,
                              Teuchos::RCP<panzer::GlobalIndexer> globalIndexer,
                              const bool overrideNoxOutput)
{
  using namespace Teuchos;
  using namespace Tempus;

  RCP<ParameterList> tempus_pl = parameterList(input_params->sublist("Solution Control").sublist("Tempus"));

  if (overrideNoxOutput) {
    auto nox_utils = Teuchos::rcp(new NOX::Utils(NOX::Utils::Error,comm->getRank(),0,3));
    auto print_nonlinear = Teuchos::rcp(new NOX::ObserverPrint(nox_utils,10));
    tempus_pl->sublist("My Example Stepper",true)
      .sublist("My Example Solver",true)
      .sublist("NOX",true)
      .sublist("Solver Options",true)
      .set<RCP<NOX::Observer>>("Observer",print_nonlinear);
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
    RCP<const panzer_stk::TempusObserverFactory> tof =
      Teuchos::rcp(new user_app::TempusObserverFactory(stkIOResponseLibrary,rLibrary->getWorksetContainer()));
    tempus_observers->addObserver(Teuchos::rcp(new Tempus::IntegratorObserverBasic<double>));
    tempus_observers->addObserver(tof->buildTempusObserver(mesh,globalIndexer,linObjFactory));
    integrator->setObserver(tempus_observers);
  }

  RCP<Thyra::VectorBase<double>> x0 = physics->getNominalValues().get_x()->clone_v();
  integrator->initializeSolutionHistory(0.0, x0);

  return integrator;
}

std::tuple<int,int> user_app::findParameterIndex(const std::string& p_name,const Thyra::ModelEvaluator<double>& me)
{
  bool found = false;
  int index = -1;
  int sub_index = -1;
  for (int p = 0; p < me.Np(); ++p) {
    auto p_names = me.get_p_names(p);
    for (Teuchos::Ordinal p_sub=0; p_sub < p_names->size(); ++p_sub) {
      if ((*p_names)[p_sub] == p_name) {
        found = true;
        index = p;
        sub_index = p_sub;
        break;
      }
    }
    if (found) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!found,std::runtime_error,
    "ERROR: the parameter \"" << p_name << "\" is not a valid parameter name in the model evaluator!");

  return {index,sub_index};
}

std::tuple<int,int> user_app::findResponseIndex(const std::string& g_name,const Thyra::ModelEvaluator<double>& me)
{
  bool found = false;
  int index = -1;
  int sub_index = -1;
  for (int g = 0; g < me.Ng(); ++g) {
    auto g_names = me.get_g_names(g);
    for (Teuchos::Ordinal g_sub=0; g_sub < g_names.size(); ++g_sub) {
      std::cout << "g(" << g << "," << g_sub << ")=" << g_names[g_sub] << std::endl;
      if (g_names[g_sub] == g_name) {
        found = true;
        index = g;
        sub_index = g_sub;
        break;
      }
    }
    if (found) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!found,std::runtime_error,
    "ERROR: the response \"" << g_name << "\" is not a valid response name in the model evaluator!");
  
  return {index,sub_index};
}
