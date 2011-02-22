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

#include "Kokkos_DefaultNode.hpp"

#include "Panzer_ConfigDefs.hpp"
#include "Panzer_STK_ModelEvaluatorFactory_Epetra.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"

#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include "user_app_NOXObserverFactory_Epetra.hpp"
#include "user_app_RythmosObserverFactory_Epetra.hpp"

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

  try {
    
    Teuchos::RCP<Teuchos::Time> total_time = 
      Teuchos::TimeMonitor::getNewTimer("Helios: Total Time");
    
    Teuchos::TimeMonitor timer(*total_time); 
    
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    
    // Parse the command line arguments
    std::string input_file_name = "user_app.xml";
    {
      Teuchos::CommandLineProcessor clp;
      
      clp.setOption("i", &input_file_name, "User_App input xml filename");
      
      Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = 
	clp.parse(argc,argv,&std::cerr);
      
      TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL, 
			 std::runtime_error, "Failed to parse command line!");
    }

    // Parse the input file and broadcast to other processes
    Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("User_App Parameters"));
    Teuchos::updateParametersFromXmlFileAndBroadcast(input_file_name, input_params.get(), *comm);
    
    *out << *input_params << std::endl;
    
    // Add in the application specific equation set factory
    Teuchos::RCP<const panzer::EquationSetFactory> eqset_factory = 
      Teuchos::rcp(new user_app::MyFactory);
    input_params->sublist("Assembly").set("Equation Set Factory", eqset_factory);

    // Add in the application specific closure model factory
    Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > cm_factory = 
      Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    (Teuchos::rcp_const_cast<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >(cm_factory))->buildObjects(cm_builder);
    input_params->sublist("Assembly").set("Closure Model Factory", cm_factory);

    // Add in the application specific bc factory
    Teuchos::RCP<const panzer::BCStrategyFactory> bc_factory = 
      Teuchos::rcp(new user_app::BCFactory);
    input_params->sublist("Assembly").set("BC Factory", bc_factory);

    // Add in the application specific observer factories
    {
      Teuchos::RCP<const panzer_stk::RythmosObserverFactory_Epetra> rof = 
	Teuchos::rcp(new user_app::RythmosObserverFactory_Epetra);
      input_params->sublist("Solver Factories").set("Rythmos Observer Factory", rof);
      Teuchos::RCP<const panzer_stk::NOXObserverFactory_Epetra> nof = 
	Teuchos::rcp(new user_app::NOXObserverFactory_Epetra);
      input_params->sublist("Solver Factories").set("NOX Observer Factory", nof);
    }

    Teuchos::RCP<Thyra::ModelEvaluator<double> > physics;
    Teuchos::RCP<Thyra::ModelEvaluator<double> > solver;
    {
      panzer_stk::ModelEvaluatorFactory_Epetra<double> me_factory;
      me_factory.setParameterList(input_params);
      me_factory.buildObjects(comm); 
      physics = me_factory.getPhysicsModelEvaluator();
      solver = me_factory.getResponseOnlyModelEvaluator();
    }
    
    // solve the system
    {
      

      // Set inputs
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = solver->createInArgs();
      const Thyra::ModelEvaluatorBase::InArgs<double> inArgsNominal = solver->getNominalValues();

      // Set outputs
      Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = solver->createOutArgs();
      int num_g = outArgs.Ng();
      assert (num_g == 1);  // for now only solution is returned

      // Solution vector is returned as extra respons vector
      RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*physics->get_x_space());
      outArgs.set_g(0,gx);

      // Now, solve the problem and return the responses
      solver->evalModel(inArgs, outArgs);
      
      //std::cout << *gx << std::endl;
      
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
  
  Teuchos::TimeMonitor::summarize(*out,false,true,false);

  if (status == 0)
    *out << "panzer::MainDriver run completed." << std::endl;

  return status;
}
