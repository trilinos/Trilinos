// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_PerformAnalysis.hpp"

#include "Piro_PerformSolve.hpp"

#include "Teuchos_FancyOStream.hpp"
#include <iostream>
#include <string>
#include "Thyra_DetachedVectorView.hpp"

#ifdef HAVE_PIRO_TRIKOTA
#include "TriKota_Driver.hpp"
#include "TriKota_ThyraDirectApplicInterface.hpp"
#endif

#ifdef HAVE_PIRO_MOOCHO
#include "MoochoPack_MoochoThyraSolver.hpp"
#endif

#ifdef HAVE_PIRO_OPTIPACK
#include "OptiPack_NonlinearCG.hpp"
#include "GlobiPack_BrentsLineSearch.hpp"
#endif

#ifdef HAVE_PIRO_ROL
#include "ROL_ThyraVector.hpp"
#include "ROL_Thyra_BoundConstraint.hpp"
#include "ROL_ThyraME_Objective.hpp"
#include "ROL_ThyraProductME_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#endif


using std::cout; using std::endl; using std::string;
using Teuchos::RCP; using Teuchos::rcp; using Teuchos::ParameterList;
using Teuchos::null; using Teuchos::outArg;

int
Piro::PerformAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& analysisParams,
    RCP< Thyra::VectorBase<double> >& result)
{

  analysisParams.validateParameters(*Piro::getValidPiroAnalysisParameters(),0);

  int status;
  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  string analysis = analysisParams.get<string>("Analysis Package");
  *out << "\n\nPiro::PerformAnalysis() requests: " << analysis << endl;

  if (analysis=="Solve") {
    *out << "Piro PerformAnalysis: Model Solve Being Performed " << endl;
    Piro::PerformSolveBase(piroModel, analysisParams.sublist("Solve"), result);
    status = 0; // Succeeds or throws
  }
#ifdef HAVE_PIRO_TRIKOTA
  else if (analysis=="Dakota") {
    *out << "Piro PerformAnalysis: Dakota Analysis Being Performed " << endl;

    status = Piro::PerformDakotaAnalysis(piroModel,
                         analysisParams.sublist("Dakota"), result);

  }
#endif
#ifdef HAVE_PIRO_MOOCHO
  else if (analysis == "MOOCHO") {
    *out << "Piro PerformAnalysis: MOOCHO Optimization Being Performed " << endl;
    status = Piro::PerformMoochoAnalysis(piroModel,
                          analysisParams.sublist("MOOCHO"), result);

  }
#endif
#ifdef HAVE_PIRO_OPTIPACK
  else if (analysis == "OptiPack") {
    *out << "Piro PerformAnalysis: Optipack Optimization Being Performed " << endl;
    status = Piro::PerformOptiPackAnalysis(piroModel,
                    analysisParams.sublist("OptiPack"),
                    analysisParams.sublist("GlobiPack"), result);

  }
#endif
#ifdef HAVE_PIRO_ROL
  else if (analysis == "ROL") {
    *out << "Piro PerformAnalysis: ROL Optimization Being Performed " << endl;
    status = Piro::PerformROLAnalysis(piroModel,
                          analysisParams, result);

  }
#endif
  else {
    if (analysis == "Dakota" || analysis == "OptiPack" || analysis == "MOOCHO" || analysis == "ROL")
      *out << "ERROR: Trilinos/Piro was not configured to include \n "
           << "       analysis type: " << analysis << endl;
    else
      *out << "ERROR: Piro: Unknown analysis type: " << analysis << "\n"
           << "       Valid analysis types are: Solve, Dakota, MOOCHO, OptiPack, ROL\n" << endl;
    status = 0; // Should not fail tests
  }

  // Output status and paramters
  if (status==0)  *out << "\nPiro Analysis Finished successfully." << endl;
  else  *out << "\nPiro Analysis failed with status: " << status << endl;

  if ( analysisParams.get("Output Final Parameters", true) )
    if (result != Teuchos::null) {
       *out << "\tFinal parameters are: " << "\n\tp = ";
       *out << Teuchos::describe(*result, Teuchos::VERB_EXTREME ) << endl;
    }

  return status;
}

int
Piro::PerformMoochoAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& moochoParams,
    RCP< Thyra::VectorBase<double> >& p)
{
#ifdef HAVE_PIRO_MOOCHO
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
#ifdef HAVE_PIRO_TRIKOTA
  dakotaParams.validateParameters(*Piro::getValidPiroAnalysisDakotaParameters(),0);
  using std::string;

  string dakotaIn  = dakotaParams.get("Input File","dakota.in");
  string dakotaOut = dakotaParams.get("Output File","dakota.out");
  string dakotaErr = dakotaParams.get("Error File","dakota.err");
  string dakotaRes = dakotaParams.get("Restart File","dakota_restart.out");
  string dakotaRestartIn;
  if (dakotaParams.isParameter("Restart File To Read"))
    dakotaRestartIn = dakotaParams.get<string>("Restart File To Read");

  int dakotaRestartEvals= dakotaParams.get("Restart Evals To Read", 0);

  int p_index = dakotaParams.get("Parameter Vector Index", 0);
  int g_index = dakotaParams.get("Response Vector Index", 0);

  TriKota::Driver dakota(dakotaIn, dakotaOut, dakotaErr, dakotaRes,
                         dakotaRestartIn, dakotaRestartEvals);

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
#ifdef HAVE_PIRO_OPTIPACK
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

int
Piro::PerformROLAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& analysisParams,
    RCP< Thyra::VectorBase<double> >& p)
{
  auto rolParams = analysisParams.sublist("ROL");
#ifdef HAVE_PIRO_ROL
  using std::string;

  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  int g_index = rolParams.get<int>("Response Vector Index", 0);

  int num_parameters = rolParams.get<int>("Number of Parameters", 1);
  std::vector<int> p_indices(num_parameters);
  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    p_indices[i] = rolParams.get<int>(ss.str(), i);
  }

  Teuchos::Array<Teuchos::RCP<Thyra::VectorSpaceBase<double> const>> p_spaces(num_parameters);
  Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<double>>> p_vecs(num_parameters);
  for (auto i = 0; i < num_parameters; ++i) {
    p_spaces[i] = piroModel.get_p_space(p_indices[i]);
    p_vecs[i] = Thyra::createMember(p_spaces[i]);
  }
  Teuchos::RCP<Thyra::DefaultProductVectorSpace<double> const> p_space = Thyra::productVectorSpace<double>(p_spaces);
  Teuchos::RCP<Thyra::DefaultProductVector<double>> p_prod = Thyra::defaultProductVector<double>(p_space, p_vecs());
  p = p_prod;

//  p = Thyra::createMember(piroModel.get_p_space(p_index));

  for (auto i = 0; i < num_parameters; ++i) {
    RCP<const Thyra::VectorBase<double> > p_init = piroModel.getNominalValues().get_p(p_indices[i]);
    Thyra::copy(*p_init, p_prod->getNonconstVectorBlock(i).ptr());
  }

  ROL::ThyraVector<double> rol_p(p_prod);


  ROL::ThyraProductME_Objective<double> obj(piroModel, g_index, p_indices, Teuchos::rcp(&analysisParams.sublist("Optimization Status"),false));

  bool print = rolParams.get<bool>("Print Output", false);

  int seed = rolParams.get<int>("Seed For Thyra Randomize", 42);

  //! set initial guess (or use the one provided by the Model Evaluator)
  std::string init_guess_type = rolParams.get<string>("Parameter Initial Guess Type", "From Model Evaluator");
  if(init_guess_type == "Uniform Vector")
    rol_p.putScalar(rolParams.get<double>("Uniform Parameter Guess", 1.0));
  else if(init_guess_type == "Random Vector") {
    Teuchos::Array<double> minmax(2); minmax[0] = -1; minmax[1] = 1;
    minmax = rolParams.get<Teuchos::Array<double> >("Min And Max Of Random Parameter Guess", minmax);
    ::Thyra::randomize<double>( minmax[0], minmax[1],  rol_p.getVector().ptr());
  }
  else if(init_guess_type != "From Model Evaluator") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
              std::endl << "Error in Piro::PerformROLAnalysis:  " <<
              "Parameter Initial Guess Type \"" << init_guess_type << "\" is not Known.\nValid options are: \"Parameter Scalar Guess\", \"Uniform Vector\" and \"Random Vector\""<<std::endl);
  }

  //! test thyra implementation of ROL vector
  if(rolParams.get<bool>("Test Vector", false)) {
    Teuchos::RCP<Thyra::VectorBase<double> > rand_vec_x = p->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > rand_vec_y = p->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > rand_vec_z = p->clone_v();
    ::Thyra::seed_randomize<double>( seed );

    int num_tests = rolParams.get<int>("Number Of Vector Tests", 1);

    for(int i=0; i< num_tests; i++) {

      *out << "\nROL performing vector test " << i+1 << " of " << num_tests  << std::endl;

      ::Thyra::randomize<double>( -1.0, 1.0,  rand_vec_x.ptr());
      ::Thyra::randomize<double>( -1.0, 1.0,  rand_vec_y.ptr());
      ::Thyra::randomize<double>( -1.0, 1.0,  rand_vec_z.ptr());

      ROL::ThyraVector<double> rol_x(rand_vec_x);
      ROL::ThyraVector<double> rol_y(rand_vec_y);
      ROL::ThyraVector<double> rol_z(rand_vec_z);

      rol_x.checkVector(rol_y, rol_z,print, *out);
    }
  }

  //! check correctness of Gradient prvided by Model Evaluator
  if(rolParams.get<bool>("Check Gradient", false)) {
    Teuchos::RCP<Thyra::VectorBase<double> > rand_vec = p->clone_v();
    ::Thyra::seed_randomize<double>( seed );

    int num_checks = rolParams.get<int>("Number Of Gradient Checks", 1);
    double norm_p = rol_p.norm();

    for(int i=0; i< num_checks; i++) {

      *out << "\nROL performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

      ::Thyra::randomize<double>( -1.0, 1.0,  rand_vec.ptr());

      ROL::ThyraVector<double> rol_direction(rand_vec);

      double norm_d = rol_direction.norm();
      if(norm_d*norm_p > 0.0)
        rol_direction.scale(norm_p/norm_d);

      obj.checkGradient(rol_p, rol_direction, print, *out);
    }
  }

  // Define Step
  Teuchos::RCP<ROL::LineSearchStep<double> > stepLS = Teuchos::rcp(new ROL::LineSearchStep<double>(rolParams.sublist("ROL Options")));
  Teuchos::RCP<ROL::TrustRegionStep<double> > stepTR = Teuchos::rcp(new ROL::TrustRegionStep<double>(rolParams.sublist("ROL Options")));

  *out << "\nROL options:" << std::endl;
  rolParams.sublist("ROL Options").print(*out);
  *out << std::endl;


  // Define Status Test
  double gtol  = rolParams.get("Gradient Tolerance", 1e-5);  // norm of gradient tolerance
  double stol  = rolParams.get("Step Tolerance", 1e-5);  // norm of step tolerance
  int   maxit = rolParams.get("Max Iterations", 100);    // maximum number of iterations
  Teuchos::RCP<ROL::StatusTest<double> > status =
    Teuchos::rcp(new ROL::StatusTest<double>(gtol, stol, maxit));

  // Define Algorithm
  ROL::Algorithm<double> algoLS(stepLS,status,print);
  ROL::Algorithm<double> algoTR(stepTR,status,print);

  // Run Algorithm
  std::vector<std::string> output;
  if(rolParams.get<bool>("Bound Constrained", false)) {
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<double>>> p_lo_vecs(num_parameters);
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<double>>> p_up_vecs(num_parameters);
    double eps_bound = rolParams.get<double>("epsilon bound", 1e-6);
    for (auto i = 0; i < num_parameters; ++i) {
      p_lo_vecs[i] = piroModel.getLowerBounds().get_p(p_indices[i]);
      p_up_vecs[i] = piroModel.getUpperBounds().get_p(p_indices[i]);
      TEUCHOS_TEST_FOR_EXCEPTION((p_lo_vecs[i] == Teuchos::null)  || (p_up_vecs[i] == Teuchos::null), Teuchos::Exceptions::InvalidParameter,
          std::endl << "Error in Piro::PerformROLAnalysis:  " <<
          "Lower and/or Upper bounds pointers are null, cannot perform bound constrained optimization"<<std::endl);
    }
    Teuchos::RCP<const Thyra::VectorBase<double>> p_lo = Thyra::defaultProductVector<double>(p_space, p_lo_vecs());
    Teuchos::RCP<const Thyra::VectorBase<double>> p_up = Thyra::defaultProductVector<double>(p_space, p_up_vecs());

    ROL::Thyra_BoundConstraint<double>  boundConstraint(p_lo->clone_v(), p_up->clone_v(), eps_bound);
    if(rolParams.get<std::string>("Step Method", "Line Search") == "Line Search") {
      *out << "\nUsing Line Search Algorithm" << std::endl;
      output  = algoLS.run(rol_p, obj, boundConstraint, print, *out);
    }
    else {
      *out << "\nUsing Trust Region Algorithm" << std::endl;
      output  = algoTR.run(rol_p, obj, boundConstraint, print, *out);
    }
  }
  else
    if(rolParams.get<std::string>("Step Method", "Line Search") == "Line Search") {
      *out << "\nUsing Line Search Algorithm" << std::endl;
      output  = algoLS.run(rol_p, obj, print, *out);
    }
    else {
      *out << "\nUsing Trust Region Algorithm" << std::endl;
      output  = algoTR.run(rol_p, obj, print, *out);
    }


  for ( unsigned i = 0; i < output.size(); i++ ) {
    *out << output[i];
  }

  return 0;
#else
 RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
 *out << "ERROR: Trilinos/Piro was not configured to include ROL analysis."
      << "\nYou must enable ROL." << endl;
 return 0;  // should not fail tests
#endif
}


RCP<const Teuchos::ParameterList>
Piro::getValidPiroAnalysisParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("Valid Piro Analysis Params"));;

  validPL->set<std::string>("Analysis Package", "","Must be: Solve, MOOCHO, Dakota, or OptiPack.");
  validPL->set<bool>("Output Final Parameters", false, "");
  validPL->sublist("Solve",     false, "");
  validPL->sublist("MOOCHO",    false, "");
  validPL->sublist("OptiPack",  false, "");
  validPL->sublist("GlobiPack", false, "");
  validPL->sublist("Dakota",    false, "");
  validPL->sublist("ROL",       false, "");
  validPL->set<int>("Write Interval", 1, "Iterval between writes to mesh");

  return validPL;
}


RCP<const Teuchos::ParameterList>
Piro::getValidPiroAnalysisDakotaParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("Valid Piro Analysis Dakota Params"));;

  validPL->set<std::string>("Input File", "","Defaults to dakota.in");
  validPL->set<std::string>("Output File", "","Defaults to dakota.out");
  validPL->set<std::string>("Error File", "","Defaults to dakota.err");
  validPL->set<std::string>("Restart File", "","Defaults to dakota_restart.out");
  validPL->set<std::string>("Restart File To Read", "","Defaults to NULL (no restart file read)");
  validPL->set<int>("Restart Evals To Read", 0,
                    "Number of evaluations to read from restart. Defaults to 0 (all)");
  validPL->set<int>("Parameter Vector Index", 0,"");
  validPL->set<int>("Response Vector Index", 0,"");

  return validPL;
}
