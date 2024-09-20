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

#include "user_app_Utilities.hpp"

#include <string>
#include <iostream>
#include <tuple>

// *************************************************
//  Main Driver for a transient optimization problem using ROL and
//  Tempus
// *************************************************

int main(int argc, char *argv[])
{
  using namespace Teuchos;
  using namespace Tempus;

  int status = 0;

  std::ostream os_dev_null{nullptr};
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &os_dev_null);
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

    // Pull off objective and rol lists so the panzer validation does not fail
    RCP<ParameterList> objective_params = parameterList(input_params->sublist("Objective",true));
    RCP<ParameterList> rol_params = parameterList(input_params->sublist("ROL", true));
    input_params->remove("Objective");
    input_params->remove("ROL");

    {
      RCP<ParameterList> tempus_params = parameterList(input_params->sublist("Solution Control",true).sublist("Tempus",true));
      auto objective = ROL::makePtr<ROL::TransientReducedObjective<double>>(input_params,comm,objective_params,out);

      // Create target -- do forward integration with perturbed parameter values
      RCP<ROL::Vector<double>> p = objective->create_design_vector();
      p->setScalar(2.0);
      RCP<ROL::Vector<double>> r = objective->create_response_vector();
      objective->run_tempus(*r, *p);
      objective->set_target(r);

      p->setScalar(1.5);
      ROL::OptimizationProblem<double> problem(objective, p);
      ROL::OptimizationSolver<double> solver(problem, *rol_params);
      ROL::Ptr<std::ostream> outStream = ROL::makePtrFromRef(std::cout);
      if (comm->getRank() == 0)
        solver.solve(*outStream);
      else
        solver.solve(os_dev_null);

      {
        const ROL::ThyraVector<double>& thyra_p = dyn_cast<const ROL::ThyraVector<double> >(*p);
        ROL::ThyraVector<double>& thyra_r = dyn_cast<ROL::ThyraVector<double> >(*r);
        *out << "Final Values: p = " << Thyra::get_ele(*(thyra_p.getVector()),0)
             << ", g = " << Thyra::get_ele(*(thyra_r.getVector()),0)
             << std::endl;
      }
    }

    // user_app::computeAndPrintResponses(physics,x,x_dot,t,global_data,*out);

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
