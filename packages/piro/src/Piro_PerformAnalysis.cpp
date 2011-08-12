// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Piro_PerformAnalysis.hpp"
#include "Teuchos_FancyOStream.hpp"
#include <iostream>
#include <string>
#include "Thyra_DetachedVectorView.hpp"

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
Piro::PerformAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& analysisParams,
    RCP< Thyra::VectorBase<double> >& p)
{

  analysisParams.validateParameters(*Piro::getValidPiroAnalysisParameters(),0);

  int status;
  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  string analysis = analysisParams.get<string>("Analysis Package");
  *out << "\n\nPiro::PerformAnalysis() requests: " << analysis << endl;

#ifdef Piro_ENABLE_TriKota
  if (analysis=="Dakota") {
    *out << "Piro PerformAnalysis: Dakota Analysis Being Performed " << endl;

    status = Piro::PerformDakotaAnalysis(piroModel,
                         analysisParams.sublist("Dakota"), p);

  } else
#endif
#ifdef Piro_ENABLE_MOOCHO
  if (analysis == "MOOCHO") {
    *out << "Piro PerformAnalysis: MOOCHO Optimization Being Performed " << endl;
    status = Piro::PerformMoochoAnalysis(piroModel,
                          analysisParams.sublist("MOOCHO"), p);

  } else
#endif
#ifdef Piro_ENABLE_OptiPack
  if (analysis == "OptiPack") {
    *out << "Piro PerformAnalysis: Optipack Optimization Being Performed " << endl;
    status = Piro::PerformOptiPackAnalysis(piroModel,
                    analysisParams.sublist("OptiPack"),
                    analysisParams.sublist("GlobiPack"), p);

  } else
#endif
  {
    if (analysis == "Dakota" || analysis == "OptiPack" || analysis == "MOOCHO")
      *out << "ERROR: Trilinos/Piro was not configured to include \n " 
           << "       analysis type: " << analysis << endl;
    else
      *out << "ERROR: Piro: Unknown analysis type: " << analysis << "\n"
           << "       Valid analysis types are: Dakota, MOOCHO, OptiPack\n" << endl;
    status = 0; // Should not fail tests
  }

  // Output status and paramters
  if (status==0)  *out << "\nPiro Analysis Finished successfully." << endl;
  else  *out << "\nPiro Analysis failed with status: " << status << endl; 

  if ( analysisParams.get("Output Final Parameters", true) )
    if (p != Teuchos::null) {
       *out << "\tFinal parameters are: " << "\n\tp = ";
       *out << Teuchos::describe(*p, Teuchos::VERB_EXTREME ) << endl;
    }

  return status;
}

int
Piro::PerformMoochoAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& moochoParams,
    RCP< Thyra::VectorBase<double> >& p)
{
#ifdef Piro_ENABLE_MOOCHO
  MoochoPack::MoochoThyraSolver solver;

  // Set the model and parameter list
  solver.setModel(rcp(&piroModel,false));
  solver.setParameterList(rcp(&moochoParams,false));

  // Solve the NLP
  const  MoochoPack::MoochoSolver::ESolutionStatus
    solution_status = solver.solve();

  // Extract the final solution
  p = Thyra::createMember(piroModel.get_p_space(0));
  RCP<const Thyra::VectorBase<double> > p_final = solver.getFinalPoint().get_p(0);
  Thyra::copy(*p_final, p.ptr());

  return (int) solution_status;
#else
 RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
 *out << "ERROR: Trilinos/Piro was not configured to include MOOCHO analysis."
      << endl;
 return 0;  // should not fail tests
#endif
}

int
Piro::PerformDakotaAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& dakotaParams,
    RCP< Thyra::VectorBase<double> >& p)
{
#ifdef Piro_ENABLE_TriKota
  string dakotaIn = dakotaParams.get("Input File","dakota.in");
  string dakotaOut= dakotaParams.get("Output File","dakota.out");
  string dakotaErr= dakotaParams.get("Error File","dakota.err");
  string dakotaRes= dakotaParams.get("Restart File","dakota_restart.out");
  int p_index = dakotaParams.get("Parameter Vector Index", 0);
  int g_index = dakotaParams.get("Response Vector Index", 0);

  TriKota::Driver dakota(dakotaIn.c_str(), dakotaOut.c_str(),
                         dakotaErr.c_str(), dakotaRes.c_str());

  RCP<TriKota::ThyraDirectApplicInterface> trikota_interface =
    rcp(new TriKota::ThyraDirectApplicInterface
         (dakota.getProblemDescDB(), rcp(&piroModel,false), p_index, g_index), 
	false);

  dakota.run(trikota_interface.get());

  Dakota::RealVector finalValues;
  if (dakota.rankZero()) 
    finalValues = dakota.getFinalSolution().all_continuous_variables();

  // Copy Dakota parameters into Thyra
  p = Thyra::createMember(piroModel.get_p_space(p_index));
  {
      Thyra::DetachedVectorView<double> global_p(p);
      for (int i = 0; i < finalValues.length(); ++i) 
        global_p[i] = finalValues[i];
  }

  return 0;
#else
 RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
 *out << "ERROR: Trilinos/Piro was not configured to include Dakota analysis."
      << "\nYou must enable TriKota." << endl;
 return 0;  // should not fail tests
#endif
}

int
Piro::PerformOptiPackAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& optipackParams,
    Teuchos::ParameterList& globipackParams,
    RCP< Thyra::VectorBase<double> >& p)
{
   RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
#ifdef Piro_ENABLE_OptiPack
  // First, Linesearch stuff
  const RCP<GlobiPack::BrentsLineSearch<double> >
    linesearch = GlobiPack::brentsLineSearch<double>();
  const RCP<ParameterList> lsPL = rcp(&globipackParams, false);
  linesearch->setParameterList(lsPL);

  // Temporary Debug
  *out << "\nCurrent LineSearch parameters" << endl;
  lsPL->print(*out);

  // Second, Optimization stuff

  p = Thyra::createMember(piroModel.get_p_space(0));

  RCP<const Thyra::VectorBase<double> > p_init = piroModel.getNominalValues().get_p(0);

  Thyra::copy(*p_init, p.ptr());

  const RCP<OptiPack::NonlinearCG<double> > cgSolver =
    OptiPack::nonlinearCG<double>(rcp(&piroModel,false), 0, 0, linesearch);

  const RCP<ParameterList> pl = rcp(&optipackParams,false);
  cgSolver->setParameterList(pl);
  
  // Temporary Debug Info
  *out << "\nCurrent nonlinearCG parameter list" << endl;
  pl->print(*out);

  // Solve the prob
  double g_opt;  // optimal value of the response
  int numIters;  // number of iteration taken
  const OptiPack::NonlinearCGUtils::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
                       null, null, null, outArg(numIters) );

  return (int) solveResult;
#else
 *out << "ERROR: Trilinos/Piro was not configured to include OptiPack analysis."
      << endl;
 return 0;  // should not fail tests
#endif
}


RCP<const Teuchos::ParameterList>
Piro::getValidPiroAnalysisParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("Valid Piro Analysis Params"));;

  validPL->set<std::string>("Analysis Package", "","Must be: MOOCHO, Dakota, or OptiPack.");
  validPL->set<bool>("Output Final Parameters", false, "");
  validPL->sublist("MOOCHO",   false, "");
  validPL->sublist("OptiPack", false, "");
  validPL->sublist("GlobiPack", false, "");
  validPL->sublist("Dakota",   false, "");

  return validPL;
}


