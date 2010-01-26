#include "Piro_Thyra_PerformAnalysis.hpp"
#include <iostream>
#include <string>

#ifdef Piro_ENABLE_TriKota
#include "TriKota_Driver.hpp"
#include "TriKota_ThyraDirectApplicInterface.hpp"
#endif

#ifdef Piro_ENABLE_MOOCHO
#include "MoochoPack_MoochoThyraSolver.hpp"
#endif

#ifdef Piro_ENABLE_OptiPack
#include "OptiPack_NonlinearCG.hpp"
#include "GlobiPack_BrentsLineSearch.hpp"
#endif

using std::cout; using std::endl; using std::string;
using Teuchos::RCP; using Teuchos::rcp; using Teuchos::ParameterList;
using Teuchos::null; using Teuchos::outArg;

int
Piro::Thyra::PerformAnalysis(
    ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& analysisParams)
{
  int status;

  string analysis = analysisParams.get<string>("Analysis Package");

  cout << "\n\nPiro::Thyra::PerformAnalysis() requests: " << analysis << endl;

#ifdef Piro_ENABLE_TriKota
  if (analysis=="Dakota") {
    cout << "Piro PerformAnalysis: Dakota Analysis Being Performed " << endl;

    status = Piro::Thyra::PerformDakotaAnalysis(piroModel,
                         analysisParams.sublist("Dakota"));

  } else
#endif
#ifdef Piro_ENABLE_MOOCHO
  if (analysis == "MOOCHO") {
    cout << "Piro PerformAnalysis: MOOCHO Optimization Being Performed " << endl;
    status = Piro::Thyra::PerformMoochoAnalysis(piroModel,
                          analysisParams.sublist("MOOCHO"));

  } else
#endif
#ifdef Piro_ENABLE_OptiPack
  if (analysis == "OptiPack") {
    cout << "Piro PerformAnalysis: Optipack Optimization Being Performed " << endl;
    status = Piro::Thyra::PerformOptiPackAnalysis(piroModel,
                    analysisParams.sublist("OptiPack"),
                    analysisParams.sublist("GlobiPack"));

  } else
#endif
  {
    if (analysis == "Dakota" || analysis == "OptiPack" || analysis == "MOOCHO")
      cout << "ERROR: Trilinos/Piro was not configured to include \n " 
           << "       analysis type: " << analysis << endl;
    else
      cout << "ERROR: Piro: Unknown analysis type: " << analysis << "\n"
           << "       Valid analysis types are: Dakota, MOOCHO, OptiPack\n" << endl;
    status = 0; // Should not fail tests
  }

  return status;
}

int
Piro::Thyra::PerformMoochoAnalysis(
    ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& moochoParams)
{
#ifdef Piro_ENABLE_MOOCHO
  MoochoPack::MoochoThyraSolver solver;

  // Set the model and parameter list
  solver.setModel(rcp(&piroModel,false));
  solver.setParameterList(rcp(&moochoParams,false));

  // Solve the NLP
  cout << "DGM_Moocho: Calling solver.solve()" << endl;
  const  MoochoPack::MoochoSolver::ESolutionStatus
    solution_status = solver.solve();

  // Extract the final solution
  RCP<const ::Thyra::VectorBase<double> > p =
     solver.getFinalPoint().get_p(0);

  cout << "\nAfter MOOCHO optimization, parameters are: " 
       << "\n\tp = " << Teuchos::describe(*p, Teuchos::VERB_EXTREME ) << endl;

  return (int) solution_status;
#else
 cout << "ERROR: Trilinos/Piro was not configured to include MOOCHO analysis."
      << endl;
 return 0;  // should not fail tests
#endif
}

int
Piro::Thyra::PerformDakotaAnalysis(
    ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& dakotaParams)
{
#ifdef Piro_ENABLE_TriKota
  string dakotaIn = dakotaParams.get("Input File","dakota.in");
  string dakotaOut= dakotaParams.get("Output File","dakota.out");
  string dakotaErr= dakotaParams.get("Error File","dakota.err");
  string dakotaRes= dakotaParams.get("Restart File","dakota_restart.out");

  TriKota::Driver dakota(dakotaIn.c_str(), dakotaOut.c_str(),
                         dakotaErr.c_str(), dakotaRes.c_str());

  RCP<TriKota::ThyraDirectApplicInterface> trikota_interface =
    rcp(new TriKota::ThyraDirectApplicInterface
         (dakota.getProblemDescDB(), rcp(&piroModel,false)), false);

  dakota.run(trikota_interface.get());


  Dakota::RealVector finalValues =
    dakota.getFinalSolution().all_continuous_variables();

  cout << "\nAfter Dakota analysis, final parameters are: " 
       << "\n\tp = " << finalValues << endl;

  return 0;
#else
 cout << "ERROR: Trilinos/Piro was not configured to include Dakota analysis."
      << "\nYou must enable TriKota." << endl;
 return 0;  // should not fail tests
#endif
}

int
Piro::Thyra::PerformOptiPackAnalysis(
    ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& optipackParams,
    Teuchos::ParameterList& globipackParams)
{
#ifdef Piro_ENABLE_OptiPack
  // First, Linesearch stuff
  const RCP<GlobiPack::BrentsLineSearch<double> >
    linesearch = GlobiPack::brentsLineSearch<double>();
  const RCP<ParameterList> lsPL = rcp(&globipackParams, false);
  linesearch->setParameterList(lsPL);

  // Temporary Debug
  cout << "\nCurrent LineSearch parameters" << endl;
  lsPL->print(cout);

  // Second, Optimization stuff

  const RCP< ::Thyra::VectorBase<double> > p =
    ::Thyra::createMember(piroModel.get_p_space(0));

  RCP<const ::Thyra::VectorBase<double> > p_init = piroModel.getNominalValues().get_p(0);

  ::Thyra::copy(*p_init, p.ptr());

  const RCP<OptiPack::NonlinearCG<double> > cgSolver =
    OptiPack::nonlinearCG<double>(rcp(&piroModel,false), 0, 0, linesearch);

  const RCP<ParameterList> pl = rcp(&optipackParams,false);
  cgSolver->setParameterList(pl);
  
  // Temporary Debug Info
  cout << "\nCurrent nonlinearCG parameter list" << endl;
  pl->print(cout);

  cout << "\nStarting optimization parameter values  "
       << Teuchos::describe(*p, Teuchos::VERB_EXTREME ) << endl;

  // Solve the prob
  double g_opt;  // optimal value of the response
  int numIters;  // number of iteration taken
  const OptiPack::NonlinearCGUtils::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
                       null, null, null, outArg(numIters) );

  cout << "\nAfter OptiPack optimization of " << numIters 
       << "  iterations \n\tg = " << g_opt << "\n\tp = " << 
         Teuchos::describe(*p, Teuchos::VERB_EXTREME ) << endl;

  return (int) solveResult;
#else
 cout << "ERROR: Trilinos/Piro was not configured to include OptiPack analysis."
      << endl;
 return 0;  // should not fail tests
#endif
}
