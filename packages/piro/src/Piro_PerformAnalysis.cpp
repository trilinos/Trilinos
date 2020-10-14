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

#include "Piro_SteadyStateSolver.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#endif

#ifdef HAVE_PIRO_ROL
#include "ROL_ThyraVector.hpp"
#include "ROL_ScaledThyraVector.hpp"
#include "ROL_Thyra_BoundConstraint.hpp"
#include "ROL_ThyraME_Objective.hpp"
#include "ROL_ThyraProductME_Objective.hpp"
#include "ROL_ThyraProductME_Objective_SimOpt.hpp"
#include "ROL_ThyraProductME_Constraint_SimOpt.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"
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

#ifdef HAVE_PIRO_ROL
  else if (analysis == "ROL") {
    *out << "Piro PerformAnalysis: ROL Optimization Being Performed " << endl;
    status = Piro::PerformROLAnalysis(piroModel,
                          analysisParams, result);

  }
#endif
  else {
    if (analysis == "Dakota" || 
        analysis == "ROL")
      *out << "ERROR: Trilinos/Piro was not configured to include \n "
           << "       analysis type: " << analysis << endl;
    else
      *out << "ERROR: Piro: Unknown analysis type: " << analysis << "\n"
           << "       Valid analysis types are: Solve, Dakota, ROL\n" << endl;
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
  (void)piroModel;
  (void)dakotaParams;
  (void)p;
 RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
 *out << "ERROR: Trilinos/Piro was not configured to include Dakota analysis."
      << "\nYou must enable TriKota." << endl;
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

  int verbose = rolParams.get<int>("Verbosity Level", 3);
  Teuchos::EVerbosityLevel verbosityLevel;
  switch(verbose) {
    case 1: verbosityLevel= Teuchos::VERB_LOW; break;
    case 2: verbosityLevel= Teuchos::VERB_MEDIUM; break;
    case 3: verbosityLevel= Teuchos::VERB_HIGH; break;
    case 4: verbosityLevel= Teuchos::VERB_EXTREME; break;
    default: verbosityLevel= Teuchos::VERB_NONE;
  }

  if(rolParams.isParameter("Use Old Reduced Space Interface") && rolParams.get<bool>("Use Old Reduced Space Interface")) {

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

    for (auto i = 0; i < num_parameters; ++i) {
      RCP<const Thyra::VectorBase<double> > p_init = piroModel.getNominalValues().get_p(p_indices[i]);
      Thyra::copy(*p_init, p_prod->getNonconstVectorBlock(i).ptr());
    }

    ROL::ThyraVector<double> rol_p(p_prod);


    ROL::ThyraProductME_Objective<double> obj(piroModel, g_index, p_indices, Teuchos::rcp(&analysisParams.sublist("Optimization Status"),false),verbosityLevel);


    bool print = rolParams.get<bool>("Print Output", false);

    int seed = rolParams.get<int>("Seed For Thyra Randomize", 42);

    //! set initial guess (or use the one provided by the Model Evaluator)
    std::string init_guess_type = rolParams.get<string>("Parameter Initial Guess Type", "From Model Evaluator");
    if(init_guess_type == "Uniform Vector")
      rol_p.putScalar(rolParams.get<double>("Uniform Parameter Guess", 1.0));
    else if(init_guess_type == "Random Vector") {
      Teuchos::Array<double> minmax(2); minmax[0] = -1; minmax[1] = 1;
      minmax = rolParams.get<Teuchos::Array<double> >("Min And Max Of Random Parameter Guess", minmax);
      ::Thyra::randomize<double>( minmax[0], minmax[1], rol_p.getVector().ptr());
    }
    else if(init_guess_type != "From Model Evaluator") {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                std::endl << "Error in Piro::PerformROLAnalysis: " <<
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

        *out << "\nROL performing vector test " << i+1 << " of " << num_tests << std::endl;

        ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_x.ptr());
        ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_y.ptr());
        ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_z.ptr());

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

        ::Thyra::randomize<double>( -1.0, 1.0, rand_vec.ptr());

        ROL::ThyraVector<double> rol_p_direction(rand_vec);

        double norm_d = rol_p_direction.norm();
        if(norm_d*norm_p > 0.0)
          rol_p_direction.scale(norm_p/norm_d);

        obj.checkGradient(rol_p, rol_p_direction, print, *out);
      }
    }

    // Define Step
    Teuchos::RCP<ROL::LineSearchStep<double> > stepLS = Teuchos::rcp(new ROL::LineSearchStep<double>(rolParams.sublist("ROL Options")));
    Teuchos::RCP<ROL::TrustRegionStep<double> > stepTR = Teuchos::rcp(new ROL::TrustRegionStep<double>(rolParams.sublist("ROL Options")));

    *out << "\nROL options:" << std::endl;
    rolParams.sublist("ROL Options").print(*out);
    *out << std::endl;


    // Define Status Test
    double gtol  = rolParams.get("Gradient Tolerance", 1e-5); // norm of gradient tolerance
    double stol  = rolParams.get("Step Tolerance", 1e-5);     // norm of step tolerance
    int    maxit = rolParams.get("Max Iterations", 100);      // maximum number of iterations
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
      //double eps_bound = rolParams.get<double>("epsilon bound", 1e-6);
      for (auto i = 0; i < num_parameters; ++i) {
        p_lo_vecs[i] = piroModel.getLowerBounds().get_p(p_indices[i]);
        p_up_vecs[i] = piroModel.getUpperBounds().get_p(p_indices[i]);
        TEUCHOS_TEST_FOR_EXCEPTION((p_lo_vecs[i] == Teuchos::null) || (p_up_vecs[i] == Teuchos::null), Teuchos::Exceptions::InvalidParameter,
            std::endl << "Error in Piro::PerformROLAnalysis: " <<
            "Lower and/or Upper bounds pointers are null, cannot perform bound constrained optimization"<<std::endl);
      }
      Teuchos::RCP<const Thyra::VectorBase<double>> p_lo = Thyra::defaultProductVector<double>(p_space, p_lo_vecs());
      Teuchos::RCP<const Thyra::VectorBase<double>> p_up = Thyra::defaultProductVector<double>(p_space, p_up_vecs());


      ROL::ThyraVector<double> plo(p_lo->clone_v());
      ROL::ThyraVector<double> pup(p_up->clone_v());
      Teuchos::RCP<ROL::BoundConstraint<double> > boundConstraint =
          rcp( new ROL::Bounds<double>(ROL::makePtrFromRef(plo), ROL::makePtrFromRef(pup)));


      //ROL::Thyra_BoundConstraint<double> boundConstraint(p_lo->clone_v(), p_up->clone_v(), eps_bound);
      if(rolParams.get<std::string>("Step Method", "Line Search") == "Line Search") {
        *out << "\nUsing Line Search Algorithm" << std::endl;
        output = algoLS.run(rol_p, obj, *boundConstraint, print, *out);
      }
      else {
        *out << "\nUsing Trust Region Algorithm" << std::endl;
        output = algoTR.run(rol_p, obj, *boundConstraint, print, *out);
      }
    }
    else
      if(rolParams.get<std::string>("Step Method", "Line Search") == "Line Search") {
        *out << "\nUsing Line Search Algorithm" << std::endl;
        output = algoLS.run(rol_p, obj, print, *out);
      }
      else {
        *out << "\nUsing Trust Region Algorithm" << std::endl;
        output = algoTR.run(rol_p, obj, print, *out);
      }


    for ( unsigned i = 0; i < output.size(); i++ ) {
      *out << output[i];
    }

    return 0;

  } else { //This is not supported for Piro Epetra solver
    using std::string;
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double>> model;
    Teuchos::RCP<Piro::SteadyStateSolver<double>> piroSSSolver;
#ifdef HAVE_PIRO_NOX
    auto piroNOXSolver = Teuchos::rcp_dynamic_cast<Piro::NOXSolver<double>>(Teuchos::rcpFromRef(piroModel));
    if(Teuchos::nonnull(piroNOXSolver)) {
      piroSSSolver = Teuchos::rcp_dynamic_cast<Piro::SteadyStateSolver<double>>(piroNOXSolver);
      model = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getSubModel());
   } else
#endif
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Error in Piro::PerformROLAnalysis: " <<
          "only Piro::NOXSolver is currently supported for piroModel\n"
          "Set \"Use Old Reduced Space Interface\" to true in input file if using Piro::Epetra::NOXSolver"<<std::endl);
    }

    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    int g_index = rolParams.get<int>("Response Vector Index", 0);

    int num_parameters = rolParams.get<int>("Number of Parameters", 1);
    std::vector<int> p_indices(num_parameters);
    std::vector<std::string> p_names;

    for(int i=0; i<num_parameters; ++i) {
      std::ostringstream ss; ss << "Parameter Vector Index " << i;
      p_indices[i] = rolParams.get<int>(ss.str(), i);
      const auto names_array = *piroSSSolver->getModel().get_p_names(p_indices[i]);
      for (int k=0; k<names_array.size(); k++) {
        p_names.push_back(names_array[k]);
      }
    }

    auto opt_paramList = Teuchos::rcp(&analysisParams.sublist("Optimization Status"),false);
    opt_paramList->set("Parameter Names", Teuchos::rcpFromRef(p_names));

    Teuchos::Array<Teuchos::RCP<Thyra::VectorSpaceBase<double> const>> p_spaces(num_parameters);
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<double>>> p_vecs(num_parameters);
    for (auto i = 0; i < num_parameters; ++i) {
     p_spaces[i] = model->get_p_space(p_indices[i]);
     p_vecs[i] = Thyra::createMember(p_spaces[i]);
    }
    Teuchos::RCP<Thyra::DefaultProductVectorSpace<double> const> p_space = Thyra::productVectorSpace<double>(p_spaces);
    Teuchos::RCP<Thyra::DefaultProductVector<double>> p_prod = Thyra::defaultProductVector<double>(p_space, p_vecs());
    p = p_prod;

    //  p = Thyra::createMember(piroModel.get_p_space(p_index));

    for (auto i = 0; i < num_parameters; ++i) {
     RCP<const Thyra::VectorBase<double> > p_init = model->getNominalValues().get_p(p_indices[i]);
     Thyra::copy(*p_init, p_prod->getNonconstVectorBlock(i).ptr());
    }

    ROL::ThyraVector<double> rol_p(p_prod);
    //Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space;
    Teuchos::RCP<Thyra::VectorSpaceBase<double> const> x_space = model->get_x_space();

    Teuchos::RCP<Thyra::VectorBase<double>> x = Thyra::createMember(x_space);
    Thyra::copy(*model->getNominalValues().get_x(), x.ptr());

    ROL::ThyraVector<double> rol_x(x);
    Teuchos::RCP<Thyra::VectorBase<double>> lambda_vec = Thyra::createMember(x_space);
    ROL::ThyraVector<double> rol_lambda(lambda_vec);

    ThyraProductME_Objective_SimOpt<double> obj(*model, g_index, p_indices, Teuchos::rcp(&analysisParams.sublist("Optimization Status"),false),verbosityLevel);
    ThyraProductME_Constraint_SimOpt<double> constr(*model, g_index, p_indices, Teuchos::rcp(&analysisParams.sublist("Optimization Status"),false),verbosityLevel);

    constr.setSolveParameters(rolParams.sublist("ROL Options"));

    if(rolParams.isParameter("Use NOX Solver") && rolParams.get<bool>("Use NOX Solver"))
      constr.setExternalSolver(Teuchos::rcpFromRef(piroModel));
    constr.setNumResponses(piroSSSolver->num_g());


    ROL::Ptr<ROL::Objective_SimOpt<double> > obj_ptr = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint_SimOpt<double> > constr_ptr = ROL::makePtrFromRef(constr);

    ROL::Ptr<ROL::Vector<double> > rol_p_ptr = ROL::makePtrFromRef(rol_p);
    ROL::Ptr<ROL::Vector<double> > rol_x_ptr = ROL::makePtrFromRef(rol_x);
    ROL::Ptr<ROL::Vector<double> > rol_lambda_ptr = ROL::makePtrFromRef(rol_lambda);
    ROL::Reduced_Objective_SimOpt<double> reduced_obj(obj_ptr,constr_ptr,rol_x_ptr,rol_p_ptr,rol_lambda_ptr);

    bool print = rolParams.get<bool>("Print Output", false);

    int seed = rolParams.get<int>("Seed For Thyra Randomize", 42);

    //! set initial guess (or use the one provided by the Model Evaluator)
    std::string init_guess_type = rolParams.get<string>("Parameter Initial Guess Type", "From Model Evaluator");
    if(init_guess_type == "Uniform Vector")
     rol_p.putScalar(rolParams.get<double>("Uniform Parameter Guess", 1.0));
    else if(init_guess_type == "Random Vector") {
     Teuchos::Array<double> minmax(2); minmax[0] = -1; minmax[1] = 1;
     minmax = rolParams.get<Teuchos::Array<double> >("Min And Max Of Random Parameter Guess", minmax);
     ::Thyra::randomize<double>( minmax[0], minmax[1], rol_p.getVector().ptr());
    }
    else if(init_guess_type != "From Model Evaluator") {
     TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
               std::endl << "Error in Piro::PerformROLAnalysis: " <<
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

       *out << "\nROL performing vector test " << i+1 << " of " << num_tests << std::endl;

       ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_x.ptr());
       ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_y.ptr());
       ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_z.ptr());

       ROL::ThyraVector<double> rol_x(rand_vec_x);
       ROL::ThyraVector<double> rol_y(rand_vec_y); 
       ROL::ThyraVector<double> rol_z(rand_vec_z);

       rol_x.checkVector(rol_y, rol_z,print, *out);
     }
    }





    //! check correctness of Gradient prvided by Model Evaluator
    if(rolParams.get<bool>("Check Gradient", false)) {
     Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec1 = p->clone_v();
     Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec1 = x->clone_v();
     Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec2 = p->clone_v();
     Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec2 = x->clone_v();

     ::Thyra::seed_randomize<double>( seed );

     auto rol_x_zero = rol_x.clone(); rol_x_zero->zero();
     auto rol_p_zero = rol_p.clone(); rol_p_zero->zero();

     int num_checks = rolParams.get<int>("Number Of Gradient Checks", 1);
     double norm_p = rol_p.norm();
     double norm_x = rol_x.norm();

     ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x),ROL::makePtrFromRef(rol_p));

     for(int i=0; i< num_checks; i++) {

       *out << "\nROL performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

       // compute direction 1
       ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec1.ptr());
       ::Thyra::randomize<double>( -1.0, 1.0, x_rand_vec1.ptr());

       ROL::ThyraVector<double> rol_p_direction1(p_rand_vec1);
       ROL::ThyraVector<double> rol_x_direction1(x_rand_vec1);

       double norm_d = rol_p_direction1.norm();
       if(norm_d*norm_p > 0.0)
         rol_p_direction1.scale(norm_p/norm_d);
       norm_d = rol_x_direction1.norm();
       if(norm_d*norm_x > 0.0)
         rol_x_direction1.scale(norm_x/norm_d);

       ROL::Vector_SimOpt<double> sopt_vec_direction1(ROL::makePtrFromRef(rol_x_direction1),ROL::makePtrFromRef(rol_p_direction1));
       ROL::Vector_SimOpt<double> sopt_vec_direction1_x(ROL::makePtrFromRef(rol_x_direction1),rol_p_zero);
       ROL::Vector_SimOpt<double> sopt_vec_direction1_p(rol_x_zero,ROL::makePtrFromRef(rol_p_direction1));

       // compute direction 2
       ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec2.ptr());
       ::Thyra::randomize<double>( -1.0, 1.0, x_rand_vec2.ptr());

       ROL::ThyraVector<double> rol_p_direction2(p_rand_vec2);
       ROL::ThyraVector<double> rol_x_direction2(x_rand_vec2);

       norm_d = rol_p_direction2.norm();
       if(norm_d*norm_p > 0.0)
         rol_p_direction2.scale(norm_p/norm_d);
       norm_d = rol_x_direction2.norm();
       if(norm_d*norm_x > 0.0)
         rol_x_direction2.scale(norm_x/norm_d);

       ROL::Vector_SimOpt<double> sopt_vec_direction2(ROL::makePtrFromRef(rol_x_direction2),ROL::makePtrFromRef(rol_p_direction2));
       ROL::Vector_SimOpt<double> sopt_vec_direction2_x(ROL::makePtrFromRef(rol_x_direction2),rol_p_zero);
       ROL::Vector_SimOpt<double> sopt_vec_direction2_p(rol_x_zero,ROL::makePtrFromRef(rol_p_direction2));


       int num_steps = 10;
       int order = 2;

       if(rolParams.get<bool>("Expensive Derivative Checks", false)) {
         *out << "Checking Reduced Gradient Accuracy" << std::endl;
         reduced_obj.checkGradient(rol_p, rol_p, rol_p_direction1, print, *out);
       }
       // Check derivatives.

       *out << "Checking Accuracy of Objective Gradient " << std::endl;
       obj.checkGradient(sopt_vec,sopt_vec_direction1,true,*out,num_steps,order);
       *out << "Checking Accuracy of Objective Gradient in x direction" << std::endl;
       obj.checkGradient(sopt_vec,sopt_vec_direction1_x,true,*out,num_steps,order);
       *out << "Checking Accuracy of Objective Gradient in p direction" << std::endl;
       obj.checkGradient(sopt_vec,sopt_vec_direction1_p,true,*out,num_steps,order);


       *out << "Checking Accuracy of Constraint Gradient " << std::endl;
       constr.checkApplyJacobian(sopt_vec,sopt_vec_direction1,rol_x_direction1, true,*out,num_steps,order);
       *out << "Checking Accuracy of Constraint Gradient in x direction (Jacobian) " << std::endl;
       constr.checkApplyJacobian(sopt_vec,sopt_vec_direction1_x,rol_x_direction1,true,*out,num_steps,order);
       *out << "Checking Accuracy of Constraint Gradient in p direction" << std::endl;
       constr.checkApplyJacobian(sopt_vec,sopt_vec_direction1_p,rol_x_direction1,true,*out,num_steps,order);

       if(rolParams.get<bool>("Expensive Derivative Checks", false))
         constr.checkApplyAdjointJacobian(sopt_vec,rol_x_direction1,rol_x_direction1,sopt_vec,true,*out,num_steps);

       *out << "Checking Consistency of Constraint Gradient and its adjoint" << std::endl;
       constr.checkAdjointConsistencyJacobian(rol_x_direction1, sopt_vec_direction2, sopt_vec,true,*out);

       obj.update(rol_x,rol_p);
       constr.update(rol_x,rol_p);
       *out << "Checking Symmetry of objective Hessian" << std::endl;
       obj.checkHessSym(sopt_vec,sopt_vec_direction1, sopt_vec_direction2, true,*out);

       *out << "Checking Symmetry of objective Hessian (H_xx = H_xx^T)" << std::endl;
       obj.checkHessSym(sopt_vec,sopt_vec_direction1_x, sopt_vec_direction2_x, true,*out);
       *out << "Checking Symmetry of objective Hessian (H_xp = H_px^T)" << std::endl;
       obj.checkHessSym(sopt_vec,sopt_vec_direction1_x, sopt_vec_direction2_p, true,*out);
       *out << "Checking Symmetry of objective Hessian (H_pp = H_pp^T)" << std::endl;
       obj.checkHessSym(sopt_vec,sopt_vec_direction1_p, sopt_vec_direction2_p, true,*out);

       *out << "Checking Accuracy of objective Hessian" << std::endl;
       obj.checkHessVec(sopt_vec,sopt_vec_direction1,true,*out,num_steps,order);

       *out << "Checking Accuracy of constraint Hessian" << std::endl;
       constr.checkApplyAdjointHessian(sopt_vec, rol_x_direction1, sopt_vec_direction2, sopt_vec_direction2, true,*out,num_steps,order);

     }
    }

    bool useFullSpace = rolParams.get("Full Space",false);

    *out << "\nROL options:" << std::endl;
    rolParams.sublist("ROL Options").print(*out);
    *out << std::endl;


    ROL::Ptr<ROL::StatusTest<double>> status = ROL::makePtr<ROL::StatusTest<double>>(rolParams.sublist("ROL Options"));
    ROL::Ptr<ROL::Step<double>> step;
    if(rolParams.get<std::string>("Step Method", "Line Search") == "Line Search")
      step = ROL::makePtr<ROL::LineSearchStep<double>>(rolParams.sublist("ROL Options"));
    else
      step = ROL::makePtr<ROL::TrustRegionStep<double>>(rolParams.sublist("ROL Options"));
    ROL::Ptr<ROL::Algorithm<double> > algo;
    algo = ROL::makePtr<ROL::Algorithm<double>>(step, status,false);

    //this is for testing the PrimalScaledThyraVector. At the moment the scaling is set to 1, so it is not changing the dot product
    Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_p = p->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_x = x->clone_v();
    ::Thyra::put_scalar<double>( 1.0, scaling_vector_p.ptr());
    ::Thyra::put_scalar<double>( 1.0, scaling_vector_x.ptr());
    //::Thyra::randomize<double>( 0.5, 2.0, scaling_vector_p.ptr());
    //::Thyra::randomize<double>( 0.5, 2.0, scaling_vector_x.ptr());
    ROL::PrimalScaledThyraVector<double> rol_x_primal(x, scaling_vector_x);
    ROL::PrimalScaledThyraVector<double> rol_p_primal(p, scaling_vector_p);

    // Run Algorithm
    std::vector<std::string> output;
    Teuchos::RCP<ROL::BoundConstraint<double> > boundConstraint;
    bool boundConstrained = rolParams.get<bool>("Bound Constrained", false);

    if(boundConstrained) {
     Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<double>>> p_lo_vecs(num_parameters);
     Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<double>>> p_up_vecs(num_parameters);
     //double eps_bound = rolParams.get<double>("epsilon bound", 1e-6);
     for (auto i = 0; i < num_parameters; ++i) {
       p_lo_vecs[i] = piroModel.getLowerBounds().get_p(p_indices[i]);
       p_up_vecs[i] = piroModel.getUpperBounds().get_p(p_indices[i]);
       TEUCHOS_TEST_FOR_EXCEPTION((p_lo_vecs[i] == Teuchos::null) || (p_up_vecs[i] == Teuchos::null), Teuchos::Exceptions::InvalidParameter,
           std::endl << "Error in Piro::PerformROLAnalysis: " <<
           "Lower and/or Upper bounds pointers are null, cannot perform bound constrained optimization"<<std::endl);
     }
     Teuchos::RCP<Thyra::VectorBase<double>> p_lo = Thyra::defaultProductVector<double>(p_space, p_lo_vecs());
     Teuchos::RCP<Thyra::VectorBase<double>> p_up = Thyra::defaultProductVector<double>(p_space, p_up_vecs());

     //ROL::Thyra_BoundConstraint<double> boundConstraint(p_lo->clone_v(), p_up->clone_v(), eps_bound);
     boundConstraint = rcp( new ROL::Bounds<double>(ROL::makePtr<ROL::ThyraVector<double> >(p_lo), ROL::makePtr<ROL::ThyraVector<double> >(p_up)));
    }

     if ( useFullSpace ) {
       //ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x),ROL::makePtrFromRef(rol_p));
       ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x_primal),ROL::makePtrFromRef(rol_p_primal));
       auto r_ptr = rol_x.clone();
       double tol = 1e-5;
       constr.solve(*r_ptr,rol_x,rol_p,tol);
       if(boundConstrained) {
         ROL::BoundConstraint<double> u_bnd(rol_x);
         ROL::Ptr<ROL::BoundConstraint<double> > bnd = ROL::makePtr<ROL::BoundConstraint_SimOpt<double> >(ROL::makePtrFromRef(u_bnd),boundConstraint);
         ROL::OptimizationProblem<double> prob(ROL::makePtrFromRef(obj), ROL::makePtrFromRef(sopt_vec), bnd, ROL::makePtrFromRef(constr), r_ptr);
         ROL::OptimizationSolver<double> optSolver(prob, rolParams.sublist("ROL Options"));
         optSolver.solve(*out);
       } else {
         ROL::OptimizationProblem<double> prob(ROL::makePtrFromRef(obj), ROL::makePtrFromRef(sopt_vec), ROL::makePtrFromRef(constr), r_ptr);
         ROL::OptimizationSolver<double> optSolver(prob, rolParams.sublist("ROL Options"));
         optSolver.solve(*out);
       }
     } else {
       if(boundConstrained)
         output = algo->run(rol_p_primal, reduced_obj, *boundConstraint, print, *out);
       else
         output = algo->run(rol_p_primal, reduced_obj, print, *out);
     }

    for ( unsigned i = 0; i < output.size(); i++ ) {
     *out << output[i];
    }

    return 0;
  }
#else
  (void)piroModel;
  (void)p;
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

  validPL->set<std::string>("Analysis Package", "","Must be: Solve, ROL or Dakota.");
  validPL->set<bool>("Output Final Parameters", false, "");
  validPL->sublist("Solve",     false, "");
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
