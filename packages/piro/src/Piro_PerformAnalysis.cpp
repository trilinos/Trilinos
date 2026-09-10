// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_PerformAnalysis.hpp"

#include "Piro_PerformSolve.hpp"

#include "Teuchos_FancyOStream.hpp"
#include <iostream>
#include <string>
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"

#include "Piro_SteadyStateSolver.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#endif

#ifdef HAVE_PIRO_ROL
#include "ROL_ThyraVector.hpp"
#include "ROL_ScaledThyraVector.hpp"
#include "ROL_Thyra_BoundConstraint.hpp"
#include "Piro_ThyraProductME_Objective_SimOpt.hpp"
#include "Piro_ThyraProductME_Constraint_SimOpt.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TypeB_AlgorithmFactory.hpp"
#include "ROL_TypeU_AlgorithmFactory.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Solver.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Piro_CustomLBFGSSecant.hpp"
#include "ROL_LinearOpScaledThyraVector.hpp"

#ifdef HAVE_PIRO_HDSALIB
#include "Piro_HDSA_MD_ROL_Data_Interface.hpp"
#include "Piro_HDSA_MD_ROL_Elliptic_u_Prior_Interface.hpp"
#include "Piro_HDSA_MD_ROL_Elliptic_z_Prior_Interface.hpp"
#include "HDSA_MD_ROL_Opt_Prob_Interface.hpp"
#include "HDSA_MD_Posterior_Sampling.hpp"
#include "HDSA_MD_Prior_Sampling.hpp"
#include "HDSA_MD_Hessian_Analysis.hpp"
#include "HDSA_MD_Update.hpp"
#include "HDSA_MD_OED.hpp"
#endif

#endif

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TransientSolver.hpp"
#include "Piro_TempusSolver.hpp"
#ifdef HAVE_PIRO_ROL
#include "Piro_ThyraProductME_ROL_DynamicObjective.hpp"
#include "Piro_ThyraProductME_Tempus_FinalObjective.hpp"
#include "Piro_ThyraProductME_ROL_DynamicConstraint.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_ReducedDynamicStationaryControlsObjective.hpp"
#endif
#endif

#ifdef HAVE_PIRO_TEKO
#include "Teko_InverseLibrary.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_SolveInverseFactory.hpp"
#endif

using std::cout; using std::endl; using std::string;
using Teuchos::RCP; using Teuchos::rcp; using Teuchos::ParameterList;
using Teuchos::null; using Teuchos::outArg;

template<class CharT, class Traits>
const std::stringstream& 
Piro::RolOutputBuffer<CharT,Traits>::getStringStream() const {
   return ss; 
}

template<class CharT, class Traits>
int
Piro::RolOutputBuffer<CharT,Traits>::overflow(int c) {
  if (c != Traits::eof())          ss << static_cast<CharT>(c);
  if (putchar(c) == Traits::eof()) return Traits::eof();
  return c;
}


Teuchos::RCP<Thyra::ProductVectorBase<double> > Piro::createProductVector(const Teuchos::RCP<Thyra::VectorBase<double> >& vec) {
  Teuchos::RCP<Thyra::ProductVectorBase<double> > prd_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double>>(vec);
  if(Teuchos::is_null(prd_vec)) {
      Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> spaces(1);
      spaces[0] = vec->space();
      Teuchos::RCP<const Thyra::DefaultProductVectorSpace<double>> prd_space = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(spaces()));
      Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<double>>> vecs(1);
      vecs[0] = vec;
      prd_vec = Thyra::defaultProductVector<double>(prd_space, vecs());
    }
  return prd_vec;
}

int
Piro::PerformAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& piroParams,
    RCP< Thyra::VectorBase<double> >& result,
    RCP< Piro::ROL_ObserverBase<double> > observer,
     const std::vector<Teuchos::RCP< Thyra::VectorBase<double> > >& x_diff_at_samples,
     const std::vector<Teuchos::RCP< Thyra::VectorBase<double> > >& p_samples)
{
  auto analysisParams = piroParams.sublist("Analysis");
  analysisParams.validateParameters(*Piro::getValidPiroAnalysisParameters(),0);

  int analysisVerbosity = analysisParams.get<int>("Output Level",2);
  RCP<std::ostream> out;
  if(analysisVerbosity > 0)
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  else // no output
    out = Teuchos::rcp(new Teuchos::oblackholestream());

  string analysis = analysisParams.get<string>("Analysis Package");

  int status;
  if (analysis=="Solve") {
    *out << "Piro::PerformAnalysis: Model Solve Being Performed " << endl;
    Piro::PerformSolveBase(piroModel, analysisParams.sublist("Solve"), result);
    status = 0; // Succeeds or throws
  }

#ifdef HAVE_PIRO_ROL
  else if (analysis == "ROL") {
    *out << "Piro::PerformAnalysis: ROL Optimization Being Performed " << endl;
    status = Piro::PerformROLAnalysis(piroModel,
                          piroParams, result, observer);

  } else if (analysis == "HDSA") {
    *out << "Piro::PerformAnalysis: HDSA Post Optimality Analysis Being Performed " << endl;
    status = Piro::PerformHDSAAnalysis(piroModel,
                          piroParams, result, observer, x_diff_at_samples, p_samples);

  }
#endif
  else {
    if ((analysis == "ROL") || (analysis == "HDSA"))
      *out << "ERROR: Trilinos/Piro was not configured to include \n "
           << "       analysis type: " << analysis << endl;
    else
      *out << "ERROR: Piro: Unknown analysis type: " << analysis << "\n"
           << "       Valid analysis types are: Solve and ROL\n" << endl;
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
Piro::PerformROLSteadyAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& piroParams,
    RCP< Thyra::VectorBase<double> >& p,
    RCP< Piro::ROL_ObserverBase<double> > observer)
{
  auto analysisParams = piroParams.sublist("Analysis");
  int analysisVerbosity = analysisParams.get<int>("Output Level",2);

  RCP<std::ostream> out;
  if(analysisVerbosity > 0)
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  else // no output
    out = Teuchos::rcp(new Teuchos::oblackholestream());


#ifdef HAVE_PIRO_ROL

  using std::string;
  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double>> model, adjointModel;
  Teuchos::RCP<Piro::SteadyStateSolver<double>> piroSSSolver;

  auto rolParams = analysisParams.sublist("ROL");  
  int num_parameters = rolParams.get<int>("Number Of Parameters", 1);
  
#ifdef HAVE_PIRO_NOX
  auto piroNOXSolver = Teuchos::rcp_dynamic_cast<Piro::NOXSolver<double>>(Teuchos::rcpFromRef(piroModel));
  if(Teuchos::nonnull(piroNOXSolver)) {
    piroSSSolver = Teuchos::rcp_dynamic_cast<Piro::SteadyStateSolver<double>>(piroNOXSolver);

    std::vector<int> p_indices(num_parameters);

    for(int i=0; i<num_parameters; ++i) {
      std::ostringstream ss; ss << "Parameter Vector Index " << i;
      p_indices[i] = rolParams.get<int>(ss.str(), i);
    }


    Teuchos::RCP<const Thyra::ProductVectorBase<double> > prodvec_p 
      = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double>>(piroNOXSolver->getSubModel()->getNominalValues().get_p(0));

    if ( prodvec_p.is_null()) {
      model = Teuchos::rcp(new Piro::ProductModelEvaluator<double>(
        Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getSubModel()),
        p_indices));

      if (!piroNOXSolver->getAdjointSubModel().is_null()) {
        adjointModel = Teuchos::rcp(new Piro::ProductModelEvaluator<double>(
          Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getAdjointSubModel()),
          p_indices));
      }
    }
    else {
      model = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getSubModel());
      adjointModel = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getAdjointSubModel());
    }
  } else
#endif
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
        std::endl << "Piro::PerformROLSteadyAnalysis, ERROR: " <<
        "only Piro::NOXSolver is currently supported for piroModel"<<std::endl);
  }

  rolParams.validateParameters(*Piro::getValidPiroAnalysisROLParameters(num_parameters),0);

  int g_index = rolParams.get<int>("Response Vector Index", 0);  
  auto p_names = Teuchos::rcp(new std::vector<std::string>());

  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    const auto names_array = *piroSSSolver->getModel().get_p_names(0);
    for (int k=0; k<names_array.size(); k++) {
      p_names->push_back(names_array[k]);
    }
  }

  //set names of parameters in the "Optimization Status" sublist
  piroParams.sublist("Optimization Status").set("Parameter Names", p_names);

  if(rolParams.isParameter("Objective Recovery Value"))
    piroParams.sublist("Optimization Status").set("Objective Recovery Value", rolParams.get<double>("Objective Recovery Value"));

  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space = model->get_p_space(0);
  p = model->getNominalValues().get_p(0)->clone_v();

  ROL::Ptr<ROL::Vector<double> > rol_p_ptr = ROL::makePtr<ROL::ThyraVector<double>>(p);
  //Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space;
  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> x_space = model->get_x_space();

  Teuchos::RCP<Thyra::VectorBase<double>> x = Thyra::createMember(x_space);
  Thyra::copy(*model->getNominalValues().get_x(), x.ptr());

  ROL::Ptr<ROL::Vector<double> > rol_x_ptr = ROL::makePtr<ROL::ThyraVector<double>>(x);
  Teuchos::RCP<Thyra::VectorBase<double>> lambda_vec = Thyra::createMember(x_space);
  ROL::Ptr<ROL::Vector<double> > rol_lambda_ptr = ROL::makePtr<ROL::ThyraVector<double>>(lambda_vec);

  Teuchos::EVerbosityLevel analysisVerbosityLevel;
  switch(analysisVerbosity) {
    case 1: analysisVerbosityLevel= Teuchos::VERB_LOW; break;
    case 2: analysisVerbosityLevel= Teuchos::VERB_MEDIUM; break;
    case 3: analysisVerbosityLevel= Teuchos::VERB_HIGH; break;
    case 4: analysisVerbosityLevel= Teuchos::VERB_EXTREME; break;
    default: analysisVerbosityLevel= Teuchos::VERB_NONE;
  }  
  ROL::Ptr<ROL::Objective_SimOpt<double> > obj_ptr = ROL::makePtr<Piro::ThyraProductME_Objective_SimOpt<double>>(model, g_index, piroParams, analysisVerbosityLevel, observer);
  ROL::Ptr<ROL::Constraint_SimOpt<double> > constr_ptr = ROL::makePtr<Piro::ThyraProductME_Constraint_SimOpt<double>>(model, adjointModel, piroParams, analysisVerbosityLevel, observer);

  constr_ptr->setSolveParameters(rolParams.sublist("ROL Options"));

  {
    auto thyra_constr_ptr = ROL::dynamicPtrCast<Piro::ThyraProductME_Constraint_SimOpt<double>>(constr_ptr);
    if(rolParams.isParameter("Use NOX Solver") && rolParams.get<bool>("Use NOX Solver"))
      thyra_constr_ptr->setExternalSolver(Teuchos::rcpFromRef(piroModel));
    thyra_constr_ptr->setNumResponses(piroSSSolver->num_g());
  }

  ROL::Reduced_Objective_SimOpt<double> reduced_obj(obj_ptr,constr_ptr,rol_x_ptr,rol_p_ptr,rol_lambda_ptr);

  int seed = rolParams.get<int>("Seed For Thyra Randomize", 42);

  //! set initial guess (or use the one provided by the Model Evaluator)
  std::string init_guess_type = rolParams.get<string>("Parameter Initial Guess Type", "From Model Evaluator");
  if(init_guess_type == "Uniform Vector")
    rol_p_ptr->setScalar(rolParams.get<double>("Uniform Parameter Guess", 1.0));
  else if(init_guess_type == "Random Vector") {
    Teuchos::Array<double> minmax(2); minmax[0] = -1; minmax[1] = 1;
    minmax = rolParams.get<Teuchos::Array<double> >("Min And Max Of Random Parameter Guess", minmax);
    ::Thyra::randomize<double>( minmax[0], minmax[1], ROL::dynamicPtrCast<ROL::ThyraVector<double>>(rol_p_ptr)->getVector().ptr());
  }
  else if(init_guess_type != "From Model Evaluator") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
              std::endl << "Piro::PerformROLSteadyAnalysis, ERROR: " <<
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

      *out << "\nPiro::PerformROLSteadyAnalysis: Performing vector test " << i+1 << " of " << num_tests << std::endl;

      ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_x.ptr());
      ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_y.ptr());
      ::Thyra::randomize<double>( -1.0, 1.0, rand_vec_z.ptr());

      ROL::ThyraVector<double> rol_vec_x(rand_vec_x);
      ROL::ThyraVector<double> rol_vec_y(rand_vec_y);
      ROL::ThyraVector<double> rol_vec_z(rand_vec_z);

      rol_vec_x.checkVector(rol_vec_y, rol_vec_z, true, *out);
    }
  }

  //! check correctness of Gradient prvided by Model Evaluator
  if(rolParams.get<bool>("Check Derivatives", false)) {
    Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec1 = p->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec1 = x->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec2 = p->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec2 = x->clone_v();

    ::Thyra::seed_randomize<double>( seed );

    auto rol_x_zero = rol_x_ptr->clone(); rol_x_zero->zero();
    auto rol_p_zero = rol_p_ptr->clone(); rol_p_zero->zero();

    int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
    double norm_p = rol_p_ptr->norm();
    double norm_x = rol_x_ptr->norm();

    ROL::Vector_SimOpt<double> sopt_vec(rol_x_ptr,rol_p_ptr);

    for(int i=0; i< num_checks; i++) {

      *out << "\nPiro::PerformROLSteadyAnalysis: Performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

      // compute direction 1
      ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec1.ptr());
      ::Thyra::randomize<double>( -1.0, 1.0, x_rand_vec1.ptr());

      ROL::Ptr<ROL::Vector<double> > rol_p_direction1 = ROL::makePtr<ROL::ThyraVector<double>>(p_rand_vec1);
      ROL::Ptr<ROL::Vector<double> > rol_x_direction1 = ROL::makePtr<ROL::ThyraVector<double>>(x_rand_vec1);

      double norm_d = rol_p_direction1->norm();
      if(norm_d*norm_p > 0.0)
        rol_p_direction1->scale(norm_p/norm_d);
      norm_d = rol_x_direction1->norm();
      if(norm_d*norm_x > 0.0)
        rol_x_direction1->scale(norm_x/norm_d);

      ROL::Vector_SimOpt<double> sopt_vec_direction1(rol_x_direction1, rol_p_direction1);
      ROL::Vector_SimOpt<double> sopt_vec_direction1_x(rol_x_direction1, rol_p_zero);
      ROL::Vector_SimOpt<double> sopt_vec_direction1_p(rol_x_zero, rol_p_direction1);

      // compute direction 2
      ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec2.ptr());
      ::Thyra::randomize<double>( -1.0, 1.0, x_rand_vec2.ptr());

      ROL::Ptr<ROL::Vector<double> > rol_p_direction2 = ROL::makePtr<ROL::ThyraVector<double>>(p_rand_vec2);
      ROL::Ptr<ROL::Vector<double> > rol_x_direction2 = ROL::makePtr<ROL::ThyraVector<double>>(x_rand_vec2);

      norm_d = rol_p_direction2->norm();
      if(norm_d*norm_p > 0.0)
        rol_p_direction2->scale(norm_p/norm_d);
      norm_d = rol_x_direction2->norm();
      if(norm_d*norm_x > 0.0)
        rol_x_direction2->scale(norm_x/norm_d);

      ROL::Vector_SimOpt<double> sopt_vec_direction2(rol_x_direction2, rol_p_direction2);
      ROL::Vector_SimOpt<double> sopt_vec_direction2_x(rol_x_direction2,rol_p_zero);
      ROL::Vector_SimOpt<double> sopt_vec_direction2_p(rol_x_zero, rol_p_direction2);


      int num_steps = 10;
      int order = 2;

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Reduced Derivative Checks", false)) {
        *out << "Piro::PerformROLSteadyAnalysis: Checking Reduced Gradient Accuracy" << std::endl;
        reduced_obj.checkGradient(*rol_p_ptr, *rol_p_direction1, true, *out);
      }
      // Check derivatives.

      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Objective Gradient " << std::endl;
      obj_ptr->checkGradient(sopt_vec,sopt_vec_direction1,true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Objective Gradient in x direction" << std::endl;
      obj_ptr->checkGradient(sopt_vec,sopt_vec_direction1_x,true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Objective Gradient in p direction" << std::endl;
      obj_ptr->checkGradient(sopt_vec,sopt_vec_direction1_p,true,*out,num_steps,order);


      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Constraint Gradient " << std::endl;
      constr_ptr->checkApplyJacobian(sopt_vec,sopt_vec_direction1,*rol_x_direction1, true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Constraint Gradient in x direction (Jacobian) " << std::endl;
      constr_ptr->checkApplyJacobian(sopt_vec,sopt_vec_direction1_x,*rol_x_direction1,true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Constraint Gradient in p direction" << std::endl;
      constr_ptr->checkApplyJacobian(sopt_vec,sopt_vec_direction1_p,*rol_x_direction1,true,*out,num_steps,order);

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Expensive Derivative Checks", false))
        constr_ptr->checkApplyAdjointJacobian(sopt_vec,*rol_x_direction1,*rol_x_direction1,sopt_vec,true,*out,num_steps);

      *out << "Piro::PerformROLSteadyAnalysis: Checking Consistency of Constraint Gradient and its adjoint" << std::endl;
      constr_ptr->checkAdjointConsistencyJacobian(*rol_x_direction1, sopt_vec_direction2, sopt_vec,true,*out);

      obj_ptr->update(*rol_x_ptr,*rol_p_ptr,ROL::UpdateType::Temp);
      constr_ptr->update(*rol_x_ptr,*rol_p_ptr,ROL::UpdateType::Temp);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian" << std::endl;
      obj_ptr->checkHessSym(sopt_vec,sopt_vec_direction1, sopt_vec_direction2, true,*out);

      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian (H_xx = H_xx^T)" << std::endl;
      obj_ptr->checkHessSym(sopt_vec,sopt_vec_direction1_x, sopt_vec_direction2_x, true,*out);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian (H_xp = H_px^T)" << std::endl;
      obj_ptr->checkHessSym(sopt_vec,sopt_vec_direction1_x, sopt_vec_direction2_p, true,*out);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian (H_pp = H_pp^T)" << std::endl;
      obj_ptr->checkHessSym(sopt_vec,sopt_vec_direction1_p, sopt_vec_direction2_p, true,*out);

      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of objective Hessian" << std::endl;
      obj_ptr->checkHessVec(sopt_vec,sopt_vec_direction1,true,*out,num_steps,order);

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Reduced Derivative Checks", false)) {
        *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of reduced objective Hessian" << std::endl;
        reduced_obj.update(*rol_p_ptr,ROL::UpdateType::Temp);
        auto hsymCheck = reduced_obj.checkHessSym(*rol_p_ptr, *rol_p_direction1, *rol_p_direction2, false,*out);
        *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of reduced objective Hessian - output:" << std::endl;
        *out << std::right
                << std::setw(20) << "<w, H(x)v>"
                << std::setw(20) << "<v, H(x)w>"
                << std::setw(20) << "abs error"
                << "\n";
        *out << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << hsymCheck[0]
                << std::setw(20) << hsymCheck[1]
                << std::setw(20) << hsymCheck[2]
                << "\n";
        *out << "Piro::PerformROLAnalysis: Checking Accuracy of reduced objective Hessian" << std::endl;
        reduced_obj.checkHessVec(*rol_p_ptr, *rol_p_direction1,true,*out,num_steps,order);
      }

      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of constraint Hessian" << std::endl;
      constr_ptr->checkApplyAdjointHessian(sopt_vec, *rol_x_direction1, sopt_vec_direction2, sopt_vec_direction2, true,*out,num_steps,order);

    }
  }

  bool useFullSpace = rolParams.get("Full Space",false);

  if(analysisVerbosity >= 3) {
    *out << "\nPiro PerformAnalysis: ROL options:" << std::endl;
    rolParams.sublist("ROL Options").print(*out);
    *out << std::endl;
  }

  bool useCustomDotProduct = false;
  bool lumpHessianMatrix = false;
  int response_index_dotProd = -1;
  if(rolParams.isSublist("Matrix Based Dot Product")) {
    const Teuchos::ParameterList& matrixDotProductList = rolParams.sublist("Matrix Based Dot Product");
    auto matrixType = matrixDotProductList.get<std::string>("Matrix Type");
    if(matrixType == "Hessian Of Response") {
      useCustomDotProduct = true;
      response_index_dotProd = matrixDotProductList.sublist("Matrix Types").sublist("Hessian Of Response").get<int>("Response Index");
      lumpHessianMatrix = matrixDotProductList.sublist("Matrix Types").sublist("Hessian Of Response").get<bool>("Lump Matrix");
    }
    else if (matrixType == "Identity")
      useCustomDotProduct = false;
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLSteadyAnalysis, ERROR: " <<
          "Matrix Type not recognized. Available options are: \n" <<
          "\"Identity\" and \"Hessian Of Response\""<<std::endl);
    }
  }

  bool useCustomSecant = false;
  int secantMaxStorage = -1;
  double secantScaling(1.0);
  int response_index_secant = -1;
  if(rolParams.isSublist("Custom Secant")) {
    Teuchos::ParameterList customSecantList = rolParams.sublist("Custom Secant");
    secantMaxStorage = customSecantList.get<int>("Maximum Storage");
    secantScaling = customSecantList.get<double>("Scaling",1.0);
    useCustomSecant = true;
    auto type = customSecantList.get<std::string>("Type", "Limited-Memory BFGS");

    TEUCHOS_TEST_FOR_EXCEPTION(type != "Limited-Memory BFGS", Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLAnalysis, ERROR: " <<
          "Type of Custom Secant not recognized. Available options are: \n" <<
          "\"Limited-Memory BFGS\""<<std::endl);

    auto initializationType = customSecantList.get<std::string>("Initialization Type");

    if(initializationType == "Hessian Of Response") {
      response_index_secant = customSecantList.sublist("Initialization Types").sublist("Hessian Of Response").get<int>("Response Index");
    }
    else if(initializationType == "Identity") {
      response_index_secant = -1;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLAnalysis, ERROR: " <<
          "Approximate Hessian not recognized. Available options are: \n" <<
          "\"Identity\",\"Hessian Of Response\""<<std::endl);
    }
  }


  
  Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_p = Teuchos::null;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > H_dotP(Teuchos::null), invH_dotP(Teuchos::null), H_sec(Teuchos::null), invH_sec(Teuchos::null);

  #ifdef HAVE_PIRO_TEKO
  {
    if(analysisVerbosity > 2)
      *out << "\nPiro::PerformROLSteadyAnalysis: Start the computation of H_pp" << std::endl;

    Teuchos::RCP<const Piro::ProductModelEvaluator<double>> model_PME = getProductModelEvaluator(model);

    Teko::BlockedLinearOp bH_dotP, bH_sec;

    if (useCustomDotProduct && !model_PME.is_null()) {
      bH_dotP = Teko::createBlockedOp();
      model_PME->block_diagonal_hessian_22(bH_dotP, *rol_x_ptr, *rol_p_ptr, response_index_dotProd);
    }
    if(useCustomSecant && (response_index_secant != -1 ) && !model_PME.is_null()) {

      if (response_index_dotProd == response_index_secant)
        bH_sec = bH_dotP;
      else {
        bH_sec = Teko::createBlockedOp();
        model_PME->block_diagonal_hessian_22(bH_sec, *rol_x_ptr, *rol_p_ptr, response_index_secant);
      }
    }
    
    if(analysisVerbosity > 2)
      *out << "Piro::PerformROLSteadyAnalysis: End of the computation of H_pp" << std::endl;

    if (useCustomDotProduct) {
      if(lumpHessianMatrix) {
        auto ones_vector_p = p->clone_v();
        ::Thyra::put_scalar<double>( 1.0, ones_vector_p.ptr());
        auto ones_vector_p_prod = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(ones_vector_p);

        scaling_vector_p = p->clone_v();
        auto scaling_vector_p_prod = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(scaling_vector_p);
        Teko::applyOp(bH_dotP, ones_vector_p_prod, scaling_vector_p_prod);
      } else {
        int numBlocks = bH_dotP->productRange()->numBlocks();
        std::vector<Teko::LinearOp> diag(numBlocks);
        for (int i=0; i<numBlocks; ++i) {
          auto linOp = Teuchos::rcp_dynamic_cast<Thyra::LinearOpWithSolveBase<double>>(
          Teuchos::rcp_const_cast<Thyra::LinearOpBase<double>>(Teko::getBlock(i, i, bH_dotP)));
          diag[i] = Thyra::nonconstInverse(linOp);
        }
        H_dotP = Teko::toLinearOp(bH_dotP);
        invH_dotP = Teko::createBlockUpperTriInverseOp(bH_dotP, diag);
      }
    }

    if(useCustomSecant) {
      if(response_index_secant == -1 ) {// identity initialization 
        invH_sec = H_sec = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(p_space));
      } else if ((response_index_dotProd == response_index_secant) && Teuchos::nonnull(H_dotP) && Teuchos::nonnull(invH_dotP)) {
        H_sec = H_dotP;
        invH_sec = invH_dotP;
      } else {
        int numBlocks = bH_sec->productRange()->numBlocks();
        std::vector<Teko::LinearOp> diag(numBlocks);
        for (int i=0; i<numBlocks; ++i) {
          auto linOp = Teuchos::rcp_dynamic_cast<Thyra::LinearOpWithSolveBase<double>>(
          Teuchos::rcp_const_cast<Thyra::LinearOpBase<double>>(Teko::getBlock(i, i, bH_sec)));
          diag[i] = Thyra::nonconstInverse(linOp);
        }
        H_sec = Teko::toLinearOp(bH_sec);
        invH_sec = Teko::createBlockUpperTriInverseOp(bH_sec, diag);
      }
    }
  }

#else
  (void)response_index_dotProd;
  (void)response_index_secant;
  TEUCHOS_TEST_FOR_EXCEPTION(useCustomDotProduct||useCustomSecant, Teuchos::Exceptions::InvalidParameter,
      std::endl << "Piro::PerformROLSteadyAnalysis, ERROR: " <<
      "Teko is required for computing custom dot product or secant"<<std::endl);
#endif

Teuchos::RCP<ROL::ThyraVector<double>> rol_p_primal = Teuchos::rcp(new ROL::ThyraVector<double>(p));
if(useCustomDotProduct) {
  if(lumpHessianMatrix)
    rol_p_primal = Teuchos::rcp(new ROL::PrimalScaledThyraVector<double>(p, scaling_vector_p));
  else
    rol_p_primal = Teuchos::rcp(new ROL::PrimalLinearOpScaledThyraVector<double>(p, H_dotP, invH_dotP));
}

  //! check correctness of Derivatives prvided by Model Evaluator
  if(rolParams.get<bool>("Check Derivatives", false) && useCustomDotProduct) {
    Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec1 = p->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec2 = p->clone_v();

    ::Thyra::seed_randomize<double>( seed );

    int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
    double norm_p = rol_p_primal->norm();

    for(int i=0; i< num_checks; i++) {

      *out << "\nPiro::PerformROLAnalysis: Performing gradient check with user defined dot-product" << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

      // compute direction 1
      ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec1.ptr());
      

      auto rol_p_direction1 = rol_p_primal->clone();
      rol_p_direction1->set(ROL::ThyraVector<double>(p_rand_vec1));

      double norm_d = rol_p_direction1->norm();
      if(norm_d*norm_p > 0.0)
        rol_p_direction1->scale(norm_p/norm_d);

      // compute direction 2
      ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec2.ptr());


      auto rol_p_direction2 = rol_p_primal->clone();
      rol_p_direction2->set(ROL::ThyraVector<double>(p_rand_vec2));

      norm_d = rol_p_direction2->norm();
      if(norm_d*norm_p > 0.0)
        rol_p_direction2->scale(norm_p/norm_d);

      int num_steps = 10;
      int order = 2;

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Reduced Derivative Checks", false)) {
        *out << "Piro::PerformROLAnalysis: Checking Reduced Gradient Accuracy" << std::endl;
        reduced_obj.checkGradient(*rol_p_primal, *rol_p_direction1, true, *out);

        *out << "Piro::PerformROLAnalysis: Checking Symmetry of reduced objective Hessian" << std::endl;
        reduced_obj.update(*rol_p_primal,ROL::UpdateType::Temp);
        auto hsymCheck = reduced_obj.checkHessSym(*rol_p_primal, *rol_p_direction1, *rol_p_direction2, false,*out);
        *out << "Piro::PerformROLAnalysis: Checking Symmetry of reduced objective Hessian - output:" << std::endl;
        *out << std::right
                << std::setw(20) << "<w, H(x)v>"
                << std::setw(20) << "<v, H(x)w>"
                << std::setw(20) << "abs error"
                << "\n";
        *out << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << hsymCheck[0]
                << std::setw(20) << hsymCheck[1]
                << std::setw(20) << hsymCheck[2]
                << "\n";
        *out << "Piro::PerformROLAnalysis: Checking Accuracy of reduced objective Hessian" << std::endl;
        reduced_obj.checkHessVec(*rol_p_primal, *rol_p_direction1,true,*out,num_steps,order);
      }
    }
  }

  int return_status = 0;
  if(rolParams.get<bool>("Perform Optimization", true)) {
    // Run Algorithm
    Teuchos::RCP<ROL::BoundConstraint<double> > boundConstraint;
    bool boundConstrained = rolParams.get<bool>("Bound Constrained", false);

    if(boundConstrained) {
      Teuchos::RCP<Thyra::VectorBase<double>> p_lo = model->getLowerBounds().get_p(0)->clone_v();
      Teuchos::RCP<Thyra::VectorBase<double>> p_up = model->getUpperBounds().get_p(0)->clone_v();

      //ROL::Thyra_BoundConstraint<double> boundConstraint(p_lo->clone_v(), p_up->clone_v(), eps_bound);
      boundConstraint = rcp( new ROL::Bounds<double>(ROL::makePtr<ROL::ThyraVector<double> >(p_lo), ROL::makePtr<ROL::ThyraVector<double> >(p_up)));
    }

      RolOutputBuffer<char> rolOutputBuffer;
      auto rolOutputStream = Teuchos::rcp(new std::ostream (&rolOutputBuffer));
      Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(rolOutputStream);
      rolOutput->setOutputToRootOnly(0);

      

      if ( useFullSpace ) {
        //using default dot product for x
        auto sopt_vec = ROL::makePtr<ROL::Vector_SimOpt<double>>(rol_x_ptr,rol_p_primal);
        auto r_ptr = rol_x_ptr->clone();
        double tol = 1e-5;
        constr_ptr->solve(*r_ptr,*rol_x_ptr,*rol_p_ptr,tol);
        if(boundConstrained) {
          *out << "Piro::PerformROLSteadyAnalysis: Solving Full Space Bound Constrained Optimization Problem" << std::endl;
          auto u_bnd = ROL::makePtr<ROL::BoundConstraint<double>>(*rol_x_ptr);
          ROL::Ptr<ROL::BoundConstraint<double> > bnd = ROL::makePtr<ROL::BoundConstraint_SimOpt<double> >(u_bnd,boundConstraint);
          auto prob = ROL::makePtr<ROL::Problem<double>>(obj_ptr, sopt_vec);
          prob->addBoundConstraint(bnd);
          prob->addConstraint("Constraint", constr_ptr,r_ptr);
          bool lumpConstraints(false), printToStream(true);
          prob->finalize(lumpConstraints, printToStream, *rolOutput);
          ROL::Solver<double> optSolver(prob, rolParams.sublist("ROL Options"));
          optSolver.solve(*out);
          return_status = optSolver.getAlgorithmState()->statusFlag;
        } else {
          *out << "Piro::PerformROLSteadyAnalysis: Solving Full Space Unconstrained Optimization Problem" << std::endl;
          auto prob = ROL::makePtr<ROL::Problem<double>>(obj_ptr, sopt_vec);
          prob->addConstraint("Constraint", constr_ptr,r_ptr);
          bool lumpConstraints(false), printToStream(true);
          prob->finalize(lumpConstraints, printToStream, *rolOutput);
          ROL::Solver<double> optSolver(prob, rolParams.sublist("ROL Options"));
          optSolver.solve(*out);
          return_status = optSolver.getAlgorithmState()->statusFlag;
        }
      } else {
        Teuchos::RCP<CustomLBFGSSecant<double>> customSecant = useCustomSecant ? Teuchos::rcp(new CustomLBFGSSecant<double> (H_sec, invH_sec, secantMaxStorage, secantScaling)) : Teuchos::null;
        if(boundConstrained) {
          *out << "Piro::PerformROLSteadyAnalysis: Solving Reduced Space Bound Constrained Optimization Problem" << std::endl;
          auto algo = ROL::TypeB::AlgorithmFactory<double>(rolParams.sublist("ROL Options"),customSecant);
          algo->run(*rol_p_primal, reduced_obj, *boundConstraint, *rolOutput); 
          return_status = algo->getState()->statusFlag;
        }  else {
          *out << "Piro::PerformROLSteadyAnalysis: Solving Reduced Space Unconstrained Optimization Problem" << std::endl;
          auto algo = ROL::TypeU::AlgorithmFactory<double>(rolParams.sublist("ROL Options"),customSecant);
          algo->run(*rol_p_primal, reduced_obj, *rolOutput);
          return_status = algo->getState()->statusFlag;
        }
      }
      if(analysisVerbosity > 1)  //write recap of optimization convergence
        *out << rolOutputBuffer.getStringStream().str();
  }
  return return_status;
#else
  (void)piroModel;
  (void)p;
  out = Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "ERROR: Trilinos/Piro was not configured to include ROL analysis."
       << "\nYou must enable ROL." << endl;
  return 0;  // should not fail tests
#endif
}

int
Piro::PerformROLTransientAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& piroParams,
    RCP< Thyra::VectorBase<double> >& p,
    RCP< Piro::ROL_ObserverBase<double> > observer)
{
  auto analysisParams = piroParams.sublist("Analysis");
  int analysisVerbosity = analysisParams.get<int>("Output Level",2);

  RCP<std::ostream> out;
  if(analysisVerbosity > 0)
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  else // no output
    out = Teuchos::rcp(new Teuchos::oblackholestream());


#if defined(HAVE_PIRO_ROL) && defined(HAVE_PIRO_TEMPUS)

  using std::string;
  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double>> model, adjointModel;
  Teuchos::RCP<Piro::TransientSolver<double>> piroTSolver;
  auto rolParams = analysisParams.sublist("ROL");  
  int num_parameters = rolParams.get<int>("Number Of Parameters", 1);

  auto piroTempusSolver = Teuchos::rcp_dynamic_cast<Piro::TempusSolver<double>>(Teuchos::rcpFromRef(piroModel));
  if(Teuchos::nonnull(piroTempusSolver)) {
    piroTSolver = Teuchos::rcp_dynamic_cast<Piro::TransientSolver<double>>(piroTempusSolver);

    std::vector<int> p_indices(num_parameters);

    for(int i=0; i<num_parameters; ++i) {
      std::ostringstream ss; ss << "Parameter Vector Index " << i;
      p_indices[i] = rolParams.get<int>(ss.str(), i);
    }


    Teuchos::RCP<const Thyra::ProductVectorBase<double> > prodvec_p 
      = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double>>(piroTempusSolver->getSubModel()->getNominalValues().get_p(0));

    if ( prodvec_p.is_null()) {
      model = Teuchos::rcp(new Piro::ProductModelEvaluator<double>(
        Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroTempusSolver->getSubModel()),
        p_indices));

      if (!piroTempusSolver->getAdjointSubModel().is_null()) {
        adjointModel = Teuchos::rcp(new Piro::ProductModelEvaluator<double>(
          Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroTempusSolver->getAdjointSubModel()),
          p_indices));
      }
    }
    else {
      model = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroTempusSolver->getSubModel());
      adjointModel = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroTempusSolver->getAdjointSubModel());
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
        std::endl << "Piro::PerformROLTransientAnalysis, ERROR: " <<
        "only Piro::TempusSolver is currently supported for piroModel"<<std::endl);
  }

  rolParams.validateParameters(*Piro::getValidPiroAnalysisROLParameters(num_parameters),0);

  int g_index = rolParams.get<int>("Response Vector Index", 0);  
  auto p_names = Teuchos::rcp(new std::vector<std::string>());

  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    const auto names_array = *piroTSolver->getModel().get_p_names(0);
    for (int k=0; k<names_array.size(); k++) {
      p_names->push_back(names_array[k]);
    }
  }

  //set names of parameters in the "Optimization Status" sublist
  piroParams.sublist("Optimization Status").set("Parameter Names", p_names);

  if(rolParams.isParameter("Objective Recovery Value"))
    piroParams.sublist("Optimization Status").set("Objective Recovery Value", rolParams.get<double>("Objective Recovery Value"));

  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space = model->get_p_space(0);
  p = model->getNominalValues().get_p(0)->clone_v();

  ROL::Ptr<ROL::Vector<double> > rol_p_ptr = ROL::makePtr<ROL::ThyraVector<double>>(p);
  //Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space;
  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> x_space = model->get_x_space();

  Teuchos::RCP<Thyra::VectorBase<double>> x = Thyra::createMember(x_space);
  Thyra::copy(*model->getNominalValues().get_x(), x.ptr());

  ROL::Ptr<ROL::Vector<double> > rol_x_ptr = ROL::makePtr<ROL::ThyraVector<double>>(x);
  Teuchos::RCP<Thyra::VectorBase<double>> lambda_vec = Thyra::createMember(x_space);
  ROL::Ptr<ROL::Vector<double> > rol_lambda_ptr = ROL::makePtr<ROL::ThyraVector<double>>(lambda_vec);

  Teuchos::EVerbosityLevel analysisVerbosityLevel;
  switch(analysisVerbosity) {
    case 1: analysisVerbosityLevel= Teuchos::VERB_LOW; break;
    case 2: analysisVerbosityLevel= Teuchos::VERB_MEDIUM; break;
    case 3: analysisVerbosityLevel= Teuchos::VERB_HIGH; break;
    case 4: analysisVerbosityLevel= Teuchos::VERB_EXTREME; break;
    default: analysisVerbosityLevel= Teuchos::VERB_NONE;
  }

  auto tempus_params = Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList(piroParams.sublist("Tempus")));

  std::string integratorName = tempus_params->get<std::string>("Integrator Name");
  double t_0 = tempus_params->sublist(integratorName).sublist("Time Step Control").get<double>("Initial Time");
  double t_f = tempus_params->sublist(integratorName).sublist("Time Step Control").get<double>("Final Time");
  double dt = tempus_params->sublist(integratorName).sublist("Time Step Control").get<double>("Initial Time Step");
  int nt = (t_f-t_0)/dt;
  auto timeStamps = ROL::TimeStamp<double>::make_uniform(t_0,t_f,{0.0,1.0},nt);


  // Create FORWARD Tempus Integrator from ModelEvaluator.
  Teuchos::RCP<Tempus::Integrator<double>> forward_integrator =
    Tempus::createIntegratorBasic<double>(tempus_params, model);

  // Create ADJOINT Tempus Integrator from ModelEvaluator.
  Teuchos::RCP<Tempus::Integrator<double>> adjoint_integrator =
    Tempus::createIntegratorBasic<double>(tempus_params, adjointModel);

  int seed = rolParams.get<int>("Seed For Thyra Randomize", 42);

  //! set initial guess (or use the one provided by the Model Evaluator)
  std::string init_guess_type = rolParams.get<string>("Parameter Initial Guess Type", "From Model Evaluator");
  if(init_guess_type == "Uniform Vector")
    rol_p_ptr->setScalar(rolParams.get<double>("Uniform Parameter Guess", 1.0));
  else if(init_guess_type == "Random Vector") {
    Teuchos::Array<double> minmax(2); minmax[0] = -1; minmax[1] = 1;
    minmax = rolParams.get<Teuchos::Array<double> >("Min And Max Of Random Parameter Guess", minmax);
    ::Thyra::randomize<double>( minmax[0], minmax[1], ROL::dynamicPtrCast<ROL::ThyraVector<double>>(rol_p_ptr)->getVector().ptr());
  }
  else if(init_guess_type != "From Model Evaluator") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
              std::endl << "Piro::PerformROLSteadyAnalysis, ERROR: " <<
              "Parameter Initial Guess Type \"" << init_guess_type << "\" is not Known.\nValid options are: \"Parameter Scalar Guess\", \"Uniform Vector\" and \"Random Vector\""<<std::endl);
  }

  // Three options are curently used to determine how Piro
  // should solve a transient optimization problem and not all of
  // the possible combinations are implemented at this time.
  //
  // 1. The first option is whether a full space or reduced space approach
  // should be used. Currently, Piro only supports reduced space 
  // approaches for transient problem. The option is kept just to be consistent
  // with steady state optimization problems.
  //
  // 2. The second option is whether we use Tempus to compute the response and its
  // total derivatives with respect to the parameters or if we use ROL to compute the
  // response and the derivative.
  // In the first approach ROL is not aware that we solve a transient problem.
  //
  // 3. The third option is whether the response depends only on the final time step or
  // if it is integrated over time.
  // Currently, both cases of the option 2 support a response that depends only on 
  // the final time step. However, only the "Tempus Driver" set to false support the 
  // time integrated response for now.
  //
  bool useFullSpace = rolParams.get("Full Space",false);
  bool useTempusDriver = rolParams.get("Tempus Driver",false);
  bool useFinalTimeStepResponse = rolParams.get<bool>("Response Depends Only On Final Time");

  if(analysisVerbosity >= 3) {
    *out << "\nPiro PerformAnalysis: ROL options:" << std::endl;
    rolParams.sublist("ROL Options").print(*out);
    *out << std::endl;
  }

  Teuchos::RCP<ROL::BoundConstraint<double> > boundConstraint;
  bool boundConstrained = rolParams.get<bool>("Bound Constrained", false);

  if(boundConstrained) {
    Teuchos::RCP<Thyra::VectorBase<double>> p_lo = model->getLowerBounds().get_p(0)->clone_v();
    Teuchos::RCP<Thyra::VectorBase<double>> p_up = model->getUpperBounds().get_p(0)->clone_v();

    //ROL::Thyra_BoundConstraint<double> boundConstraint(p_lo->clone_v(), p_up->clone_v(), eps_bound);
    boundConstraint = rcp( new ROL::Bounds<double>(ROL::makePtr<ROL::ThyraVector<double> >(p_lo), ROL::makePtr<ROL::ThyraVector<double> >(p_up)));
  }

  int return_status = 0;

  RolOutputBuffer<char> rolOutputBuffer;
  auto rolOutputStream = Teuchos::rcp(new std::ostream(&rolOutputBuffer));
  Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(rolOutputStream);
  rolOutput->setOutputToRootOnly(0);

  Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_p = p->clone_v();
  ::Thyra::put_scalar<double>( 1.0, scaling_vector_p.ptr());
  ROL::PrimalScaledThyraVector<double> rol_p_primal(p, scaling_vector_p);

  if(useTempusDriver) {

    ROL::Ptr<ROL::Objective<double> > obj_ptr = ROL::makePtr<Piro::ThyraProductME_TempusFinalObjective<double>>(model, piroTempusSolver, forward_integrator, adjoint_integrator, adjointModel, g_index, piroParams, nt, analysisVerbosityLevel, observer);

    if ( useFullSpace ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLTransientAnalysis, ERROR: " <<
          "full space approach is currently not supported."<<std::endl);
    }

    if ( !useFinalTimeStepResponse ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLTransientAnalysis, ERROR: " <<
          "integrated response is currently not supported with the Tempus driver."<<std::endl);
    }
    
    if(rolParams.get<bool>("Perform Optimization", true)) {
      if(boundConstrained) {
        *out << "Piro::PerformROLTransientAnalysis: Solving Reduced Space Bound Constrained Optimization Problem" << std::endl;
        auto algo = ROL::TypeB::AlgorithmFactory<double>(rolParams.sublist("ROL Options"));

        std::streambuf *coutbuf = NULL;
        std::ofstream out_file;
        if(rolParams.get<bool>("Redirect Tempus Output", true)) {
          out_file.open(rolParams.get<string>("Tempus Output Filename", "log_tempus.txt"));
          coutbuf = std::cout.rdbuf();
          std::cout.rdbuf(out_file.rdbuf());
        }
        algo->run(rol_p_primal, *obj_ptr, *boundConstraint, *rolOutput);
        if(rolParams.get<bool>("Redirect Tempus Output", true)) {
          std::cout.rdbuf(coutbuf);
        }

        return_status = algo->getState()->statusFlag;
      }  else {
        *out << "Piro::PerformROLTransientAnalysis: Solving Reduced Space Unconstrained Optimization Problem" << std::endl;
        auto algo = ROL::TypeU::AlgorithmFactory<double>(rolParams.sublist("ROL Options"));
        
        std::streambuf *coutbuf = NULL;
        std::ofstream out_file;
        if(rolParams.get<bool>("", true)) {
          out_file.open(rolParams.get<string>("Tempus Output Filename", "log_tempus.txt"));
          coutbuf = std::cout.rdbuf();
          std::cout.rdbuf(out_file.rdbuf());
        }
        algo->run(rol_p_primal, *obj_ptr, *rolOutput);
        if(rolParams.get<bool>("Redirect Tempus Output", true)) {
          std::cout.rdbuf(coutbuf);
        }

        return_status = algo->getState()->statusFlag;
      }
      if (return_status == ROL::EExitStatus::EXITSTATUS_STEPTOL) return_status = 0;
      if(analysisVerbosity > 1)  //write recap of optimization convergence
        *out << rolOutputBuffer.getStringStream().str();
    }

    //! check correctness of Gradient prvided by Model Evaluator
    if(rolParams.get<bool>("Check Derivatives", false)) {
      Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec1 = p->clone_v();
      Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec1 = x->clone_v();
      Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec2 = p->clone_v();
      Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec2 = x->clone_v();

      ::Thyra::seed_randomize<double>( seed );

      auto rol_x_zero = rol_x_ptr->clone(); rol_x_zero->zero();
      auto rol_p_zero = rol_p_ptr->clone(); rol_p_zero->zero();

      int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
      double norm_p = rol_p_ptr->norm();

      ROL::Vector_SimOpt<double> sopt_vec(rol_x_ptr,rol_p_ptr);

      for(int i=0; i< num_checks; i++) {

        *out << "\nPiro::PerformROLTransientAnalysis: Performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

        // compute direction 1
        ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec1.ptr());

        ROL::Ptr<ROL::Vector<double> > rol_p_direction1 = ROL::makePtr<ROL::ThyraVector<double>>(p_rand_vec1);

        double norm_d = rol_p_direction1->norm();
        if(norm_d*norm_p > 0.0)
          rol_p_direction1->scale(norm_p/norm_d);

        *out << "Piro::PerformROLTransientAnalysis: Checking Reduced Gradient Accuracy" << std::endl;

        const double ten(10);
        std::vector<double> steps(ROL_NUM_CHECKDERIV_STEPS);
        for(int li=0;li<ROL_NUM_CHECKDERIV_STEPS;++li) {
          steps[li] = pow(ten,static_cast<double>(-li-1));
        }

        std::streambuf *coutbuf = NULL;
        std::ofstream out_file;
        if(rolParams.get<bool>("Redirect Tempus Output", true)) {
          out_file.open(rolParams.get<string>("Tempus Output Filename", "log_tempus.txt"));
          coutbuf = std::cout.rdbuf();
          std::cout.rdbuf(out_file.rdbuf());
        }
        obj_ptr->checkGradient(rol_p_primal, rol_p_primal.dual(), *rol_p_direction1, steps, true, *rolOutput, 1);
        if(rolParams.get<bool>("Redirect Tempus Output", true)) {
          std::cout.rdbuf(coutbuf);
        }

        if(analysisVerbosity > 1)  //write recap of optimization convergence
          *out << rolOutputBuffer.getStringStream().str();
      }
    }

    return return_status;
  }
  else {
    ROL::Ptr<ROL::DynamicObjective<double> > obj_ptr = ROL::makePtr<Piro::ThyraProductME_ROL_DynamicObjective<double>>(model, forward_integrator, adjoint_integrator, adjointModel, g_index, piroParams, nt, useFinalTimeStepResponse, true, analysisVerbosityLevel, observer);
    ROL::Ptr<ROL::DynamicConstraint<double> > constr_ptr = ROL::makePtr<Piro::ThyraProductME_ROL_DynamicConstraint<double>>(forward_integrator, adjoint_integrator, adjointModel, piroParams, analysisVerbosityLevel, observer);

    constr_ptr->setSolveParameters(rolParams.sublist("ROL Options"));
    ROL::dynamicPtrCast<Piro::ThyraProductME_ROL_DynamicConstraint<double>>(constr_ptr)->setNumResponses(piroTSolver->num_g());


    if ( useFullSpace ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLTransientAnalysis, ERROR: " <<
          "full space approach is currently not supported."<<std::endl);
    }
    else {
      auto reduced_obj_ptr = ROL::makePtr<ROL::ReducedDynamicObjective<double>>(obj_ptr,constr_ptr,rol_x_ptr,rol_p_ptr,rol_lambda_ptr, *timeStamps, piroParams);

      ROL::ReducedDynamicStationaryControlsObjective<double> reduced_stationarycontrols_obj(reduced_obj_ptr, rol_p_ptr, nt);

      if(rolParams.get<bool>("Perform Optimization", true)) {
        if(boundConstrained) {
          *out << "Piro::PerformROLTransientAnalysis: Solving Reduced Space Bound Constrained Optimization Problem" << std::endl;
          auto algo = ROL::TypeB::AlgorithmFactory<double>(rolParams.sublist("ROL Options"));
          algo->run(rol_p_primal, reduced_stationarycontrols_obj, *boundConstraint, *rolOutput); 
          return_status = algo->getState()->statusFlag;
        }  else {
          *out << "Piro::PerformROLTransientAnalysis: Solving Reduced Space Unconstrained Optimization Problem" << std::endl;
          auto algo = ROL::TypeU::AlgorithmFactory<double>(rolParams.sublist("ROL Options"));
          algo->run(rol_p_primal, reduced_stationarycontrols_obj, *rolOutput);
          return_status = algo->getState()->statusFlag;
        }
        if (return_status == ROL::EExitStatus::EXITSTATUS_STEPTOL) return_status = 0;
        if(analysisVerbosity > 1)  //write recap of optimization convergence
          *out << rolOutputBuffer.getStringStream().str();
      }

      //! check correctness of Gradient prvided by Model Evaluator
      if(rolParams.get<bool>("Check Derivatives", false)) {
        Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec1 = p->clone_v();
        Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec1 = x->clone_v();
        Teuchos::RCP<Thyra::VectorBase<double> > p_rand_vec2 = p->clone_v();
        Teuchos::RCP<Thyra::VectorBase<double> > x_rand_vec2 = x->clone_v();

        ::Thyra::seed_randomize<double>( seed );

        auto rol_x_zero = rol_x_ptr->clone(); rol_x_zero->zero();
        auto rol_p_zero = rol_p_ptr->clone(); rol_p_zero->zero();

        int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
        double norm_p = rol_p_ptr->norm();

        ROL::Vector_SimOpt<double> sopt_vec(rol_x_ptr,rol_p_ptr);

        for(int i=0; i< num_checks; i++) {

          *out << "\nPiro::PerformROLTransientAnalysis: Performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

          // compute direction 1
          ::Thyra::randomize<double>( -1.0, 1.0, p_rand_vec1.ptr());

          ROL::Ptr<ROL::Vector<double> > rol_p_direction1 = ROL::makePtr<ROL::ThyraVector<double>>(p_rand_vec1);

          double norm_d = rol_p_direction1->norm();
          if(norm_d*norm_p > 0.0)
            rol_p_direction1->scale(norm_p/norm_d);

          *out << "Piro::PerformROLTransientAnalysis: Checking Reduced Gradient Accuracy" << std::endl;

          const double ten(10);
          std::vector<double> steps(ROL_NUM_CHECKDERIV_STEPS);
          for(int li=0;li<ROL_NUM_CHECKDERIV_STEPS;++li) {
            steps[li] = pow(ten,static_cast<double>(-li-1));
          }
          reduced_stationarycontrols_obj.checkGradient(rol_p_primal, rol_p_primal.dual(), *rol_p_direction1, steps, true, *rolOutput, 1);

          if(analysisVerbosity > 1)  //write recap of optimization convergence
            *out << rolOutputBuffer.getStringStream().str();
        }
      }
    }

    return return_status;
  }
#else
  (void)piroModel;
  (void)p;
  out = Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "ERROR: Trilinos/Piro was not configured to include ROL analysis."
       << "\nYou must enable ROL and Tempus." << endl;
  return 0;  // should not fail tests
#endif
}

int
Piro::PerformROLAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& piroParams,
    RCP< Thyra::VectorBase<double> >& p,
    RCP< Piro::ROL_ObserverBase<double> > observer)
{
  auto analysisParams = piroParams.sublist("Analysis");
  bool transient = analysisParams.get<bool>("Transient", false);
  if ( transient )
    return PerformROLTransientAnalysis(piroModel, piroParams, p, observer);
  return PerformROLSteadyAnalysis(piroModel, piroParams, p, observer);
}


#if defined(HAVE_PIRO_ROL) && defined(HAVE_PIRO_HDSALIB)

HDSA::Ptr<HDSA::Dense_Matrix<double>> Vector_To_Dense(const HDSA::Vector<double>& x) {
  HDSA::Ptr<HDSA::Dense_Matrix<double>> y = HDSA::makePtr<HDSA::Dense_Matrix<double>>(x.Dimension(), 1);

  y->Zeros();

  for (int i = 0; i < x.Dimension(); ++i) {
    y->Set_Entry(i, 0, x.Get_Entry(i));
  }

  return y;
}


HDSA::Ptr<HDSA::Dense_Matrix<double>> Append_Betas(const HDSA::Dense_Matrix<double>& old_betas,
                                                  const HDSA::Dense_Matrix<double>& new_betas) {

  const int old_len = old_betas.Number_of_Rows() * old_betas.Number_of_Columns();
  const int new_len = new_betas.Number_of_Rows() * new_betas.Number_of_Columns();

  HDSA::Ptr<HDSA::Dense_Matrix<double>> all_betas = HDSA::makePtr<HDSA::Dense_Matrix<double>>(old_len + new_len, 1);

  all_betas->Zeros();
  {
    int nrows = old_betas.Number_of_Rows();
    for (int i = 0; i < old_len; ++i) 
      all_betas->Set_Entry(i, 0, old_betas(i % nrows, i / nrows));
  }

  {
    int nrows = new_betas.Number_of_Rows();
    for (int i = 0; i < new_len; ++i)
      all_betas->Set_Entry(old_len + i, 0, new_betas(i % nrows, i / nrows));
  }
  
  return all_betas;
}

HDSA::Ptr<HDSA::Dense_Matrix<double>>
Compute_Beta_From_z(
    const HDSA::Vector<double>& z,
    const HDSA::Vector<double>& z_opt,
    const HDSA::MD_z_Prior_Interface<double>& z_prior_interface,
    const HDSA::MD_OED<double>::Offline_Data& offline)
{
  const int r = offline.r;

  // dz = z - z_opt
  HDSA::Ptr<HDSA::Vector<double>> dz = z_opt.Clone();
  dz->Set(z);
  dz->Scaled_Plus(-1.0, z_opt);

  // M_z dz
  HDSA::Ptr<HDSA::Vector<double>> Mz_dz = z_opt.Clone();
  z_prior_interface.Apply_M_z(*Mz_dz, *dz);

  // rhs = V^T M_z dz
  HDSA::Dense_Matrix<double> rhs(r, 1);
  rhs.Zeros();

  for (int i = 0; i < r; ++i) {
    rhs.Set_Entry(i, 0, (*offline.V)[i]->Dot(*Mz_dz));
  }

  // Solve
  //
  //   (V^T M_z V) beta = V^T M_z (z - z_opt)
  //
  HDSA::Ptr<HDSA::Dense_Matrix<double>> beta =
      HDSA::makePtr<HDSA::Dense_Matrix<double>>(r, 1);
  beta->Zeros();

  HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<double>(
      *offline.Vt_Mz_V, *beta, rhs);

  return beta;
}

double M_z_Norm_Difference(const HDSA::Vector<double>& a, const HDSA::Vector<double>& b,
                          const HDSA::MD_z_Prior_Interface<double>& z_prior_interface) {
  HDSA::Ptr<HDSA::Vector<double>> diff = b.Clone();
  diff->Set(a);
  diff->Scaled_Plus(static_cast<double>(-1), b);

  HDSA::Ptr<HDSA::Vector<double>> Mz_diff = b.Clone();
  z_prior_interface.Apply_M_z(*Mz_diff, *diff);

  const double norm_sq = diff->Dot(*Mz_diff);

  return std::sqrt(std::max(static_cast<double>(0), norm_sq));
}

double Estimate_Trace_Wz_Inverse_Mz(const HDSA::Vector<double>& prototype,
    const HDSA::MD_z_Prior_Interface<double>& z_prior_interface, const int num_samples)
{
  auto xi = prototype.Clone();
  auto Mz_xi = prototype.Clone();
  auto Wz_inv_Mz_xi = prototype.Clone();

  double trace_val = 0.0;

  for (int j = 0; j < num_samples; ++j) {

    // Gaussian N(0,1)
    xi->Randomize_Standard_Normal();

    // Convert Gaussian entries to Rademacher {-1, +1}.
    auto xi_rol = Teuchos::rcp_dynamic_cast<HDSA::ROL_Vector<double>>(xi);
    xi_rol->rol_vec->applyUnary(ROL::Elementwise::Sign<double>());

    z_prior_interface.Apply_M_z(*Mz_xi, *xi);

    z_prior_interface.Apply_W_z_Inverse(*Wz_inv_Mz_xi,*Mz_xi);

    trace_val += xi->Dot(*Wz_inv_Mz_xi);
  }
  return trace_val / num_samples;
}

#endif

struct OperatorInfo {
  bool use_identity = true;
  int response_index = -1;
  int param_index = -1;
};

OperatorInfo
getOperatorInfo(
    const Teuchos::ParameterList& hdsaParams,
    const std::string& operatorName,
    const bool param_index_required = false)
{
  OperatorInfo info;

  if (!hdsaParams.isSublist(operatorName))
    return info;

  const Teuchos::ParameterList& operatorList = hdsaParams.sublist(operatorName);

  const auto operatorType = operatorList.get<std::string>("Operator Type");

  if (operatorType == "Hessian Of Response") {
    info.use_identity = false;
    info.response_index = operatorList.get<int>("Response Index");
    if(param_index_required)
      info.param_index =  operatorList.get<int>("Parameter Index");
  }
  else if (operatorType == "Identity") {
    info.use_identity = true;
  }  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,   Teuchos::Exceptions::InvalidParameter, 
        std::endl  << "Piro::PerformHDSASteadyAnalysis, ERROR: " << 
        "Operator Type for \"" << operatorName << "\" not recognized. Available options are:\n" <<
        "\"Identity\" and \"Hessian Of Response\""  << std::endl);
  }

  return info;
}


int
Piro::PerformHDSAAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& piroParams,
    RCP< Thyra::VectorBase<double> >& p,
    RCP< Piro::ROL_ObserverBase<double> > observer,
    const std::vector<Teuchos::RCP< Thyra::VectorBase<double> > >& x_diff_at_samples,
    const std::vector<Teuchos::RCP< Thyra::VectorBase<double> > >& p_samples)
{

  auto analysisParams = piroParams.sublist("Analysis");
  int analysisVerbosity = analysisParams.get<int>("Output Level",2);

  RCP<std::ostream> out;
  if(analysisVerbosity > 0)
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  else // no output
    out = Teuchos::rcp(new Teuchos::oblackholestream());

#if defined(HAVE_PIRO_ROL) && defined(HAVE_PIRO_HDSALIB)

  Teuchos::EVerbosityLevel analysisVerbosityLevel;
  switch(analysisVerbosity) {
    case 1: analysisVerbosityLevel= Teuchos::VERB_LOW; break;
    case 2: analysisVerbosityLevel= Teuchos::VERB_MEDIUM; break;
    case 3: analysisVerbosityLevel= Teuchos::VERB_HIGH; break;
    case 4: analysisVerbosityLevel= Teuchos::VERB_EXTREME; break;
    default: analysisVerbosityLevel= Teuchos::VERB_NONE;
  }

  using std::string;
  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double>> model, adjointModel;
  Teuchos::RCP<Piro::SteadyStateSolver<double>> piroSSSolver;

  auto hdsaParams = analysisParams.sublist("HDSA");  
  if(analysisVerbosity >= 3) {
    *out << "\nPiro PerformAnalysis: HDSA options:" << std::endl;
    hdsaParams.sublist("HDSA Options").print(*out);
    *out << std::endl;
  }

  int num_parameters = hdsaParams.get<int>("Number Of Parameters", 1);
  std::vector<int> p_indices(num_parameters);
  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    p_indices[i] = hdsaParams.get<int>(ss.str(), i);
  }

#ifdef HAVE_PIRO_NOX
  auto piroNOXSolver = Teuchos::rcp_dynamic_cast<Piro::NOXSolver<double>>(Teuchos::rcpFromRef(piroModel));
  if(Teuchos::nonnull(piroNOXSolver)) {
    piroSSSolver = Teuchos::rcp_dynamic_cast<Piro::SteadyStateSolver<double>>(piroNOXSolver);

    Teuchos::RCP<const Thyra::ProductVectorBase<double> > prodvec_p 
      = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double>>(piroNOXSolver->getSubModel()->getNominalValues().get_p(0));

    if ( prodvec_p.is_null()) {
      model = Teuchos::rcp(new Piro::ProductModelEvaluator<double>(
        Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getSubModel()),
        p_indices));

      if (!piroNOXSolver->getAdjointSubModel().is_null()) {
        adjointModel = Teuchos::rcp(new Piro::ProductModelEvaluator<double>(
          Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getAdjointSubModel()),
          p_indices));
      }
    }
    else {
      model = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getSubModel());
      adjointModel = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(piroNOXSolver->getAdjointSubModel());
    }
  } else
#endif
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
        std::endl << "Piro::PerformROLSteadyAnalysis, ERROR: " <<
        "only Piro::NOXSolver is currently supported for piroModel"<<std::endl);
  }

  //hdsaParams.validateParameters(*Piro::getValidPiroAnalysisHDSAParameters(num_parameters),0);
  int g_index = hdsaParams.get<int>("Response Vector Index", 0);  
  auto p_names = Teuchos::rcp(new std::vector<std::string>());

  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    const auto names_array = *piroSSSolver->getModel().get_p_names(0);
    for (int k=0; k<names_array.size(); k++) {
      p_names->push_back(names_array[k]);
    }
  }

  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space = model->get_p_space(0);
  p = model->getNominalValues().get_p(0)->clone_v();
  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> x_space = model->get_x_space();
  Teuchos::RCP<Thyra::VectorBase<double>> x = model->getNominalValues().get_x()->clone_v();

  if(analysisVerbosity > 2)
    *out << "\nPiro::PerformHDSAAnalysis: Solve Constraint" << std::endl;

  ROL::Ptr<ROL::Constraint_SimOpt<double> > constr_ptr = ROL::makePtr<Piro::ThyraProductME_Constraint_SimOpt<double>>(model, adjointModel, piroParams, analysisVerbosityLevel, observer);
  constr_ptr->setSolveParameters(hdsaParams.sublist("HDSA Options"));
  {
    auto thyra_constr_ptr = ROL::dynamicPtrCast<Piro::ThyraProductME_Constraint_SimOpt<double>>(constr_ptr);
    if(hdsaParams.isParameter("Use NOX Solver") && hdsaParams.get<bool>("Use NOX Solver"))
      thyra_constr_ptr->setExternalSolver(Teuchos::rcpFromRef(piroModel));
    thyra_constr_ptr->setNumResponses(piroSSSolver->num_g());
  }
  ROL::Ptr<ROL::Vector<double> > rol_p_ptr = ROL::makePtr<ROL::ThyraVector<double>>(p);
  ROL::Ptr<ROL::Vector<double> > rol_x_ptr = ROL::makePtr<ROL::ThyraVector<double>>(x);
  auto c_ptr = rol_x_ptr->clone();
  double tol = 1e-5;
  constr_ptr->solve(*c_ptr,*rol_x_ptr,*rol_p_ptr,tol); 

  Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_p = Teuchos::null;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > priorEllOp(Teuchos::null), invPriorEllOp(Teuchos::null), invEllOp(Teuchos::null), massOp(Teuchos::null), priorMassOp(Teuchos::null), invPriorMassOp(Teuchos::null), priorMassCholOp(Teuchos::null);

  #ifdef HAVE_PIRO_TEKO
  {
    Teuchos::RCP<const Piro::ProductModelEvaluator<double>> model_PME = getProductModelEvaluator(model);
    {
      if(analysisVerbosity > 2)
        *out << "\nPiro::PerformHDSAAnalysis: Start the computation of the Prior Elliptic Operator" << std::endl;

      const auto priorEllOpInfo = getOperatorInfo(hdsaParams, "Prior Elliptic Operator");
      Teko::BlockedLinearOp bH;

      if (!priorEllOpInfo.use_identity && !model_PME.is_null()) {
        bH = Teko::createBlockedOp();
        model_PME->block_diagonal_hessian_22(bH, *rol_x_ptr, *rol_p_ptr, priorEllOpInfo.response_index);
      }
     
      if (!priorEllOpInfo.use_identity) {
        int numBlocks = bH->productRange()->numBlocks();
        std::vector<Teko::LinearOp> diag(numBlocks);
        for (int i=0; i<numBlocks; ++i) {
          auto linOp = Teuchos::rcp_dynamic_cast<Thyra::LinearOpWithSolveBase<double>>(
          Teuchos::rcp_const_cast<Thyra::LinearOpBase<double>>(Teko::getBlock(i, i, bH)));
          diag[i] = Thyra::nonconstInverse(linOp);
        }
        priorEllOp = bH;
        invPriorEllOp = Teko::createBlockUpperTriInverseOp(bH, diag, "invPriorEllOp");
      }

      if(analysisVerbosity > 2)
        *out << "Piro::PerformHDSAAnalysis: End of the computation of the Prior Elliptic Operator" << std::endl;
    }

    {
      if(analysisVerbosity > 2)
        *out << "\nPiro::PerformHDSAAnalysis: Start the computation of the Prior Mass Operator" << std::endl;
      
      const auto priorMassOpInfo = getOperatorInfo(hdsaParams, "Prior Mass Operator");
      Teko::BlockedLinearOp bH;

      if (!priorMassOpInfo.use_identity && !model_PME.is_null()) {
        bH = Teko::createBlockedOp();
        model_PME->block_diagonal_hessian_22(bH, *rol_x_ptr, *rol_p_ptr, priorMassOpInfo.response_index);
      }
     
      if (!priorMassOpInfo.use_identity) {
        int numBlocks = bH->productRange()->numBlocks();  
        std::vector<Teko::LinearOp> diag(numBlocks);
        for (int i=0; i<numBlocks; ++i) {
          auto linOp = Teuchos::rcp_dynamic_cast<Thyra::LinearOpWithSolveBase<double>>(
          Teuchos::rcp_const_cast<Thyra::LinearOpBase<double>>(Teko::getBlock(i, i, bH)));
          diag[i] = Thyra::nonconstInverse(linOp);
        }
        priorMassOp = bH;
        invPriorMassOp = Teko::createBlockUpperTriInverseOp(bH, diag, "invPriorMassOp");
      }

      if(analysisVerbosity > 2)
        *out << "Piro::PerformHDSAAnalysis: End of the computation of the Prior Mass Operator" << std::endl;
    }

    {
      if(analysisVerbosity > 2)
        *out << "\nPiro::PerformHDSAAnalysis: Start the computation of the Prior Mass Cholesky Operator" << std::endl;

      const auto priorMassCholOpInfo = getOperatorInfo(hdsaParams, "Prior Mass Cholesky Operator");
      Teko::BlockedLinearOp bH;

      if (!priorMassCholOpInfo.use_identity && !model_PME.is_null()) {
        bH = Teko::createBlockedOp();
        model_PME->block_diagonal_hessian_22(bH, *rol_x_ptr, *rol_p_ptr, priorMassCholOpInfo.response_index);
        priorMassCholOp = bH;
      }

      if(analysisVerbosity > 2)
        *out << "Piro::PerformHDSAAnalysis: End of the computation of the Prior Mass Cholesky Operator" << std::endl;
    }

    {
      if(analysisVerbosity > 2)
        *out << "\nPiro::PerformHDSAAnalysis: Start the computation of the State Elliptic Operator" << std::endl;

      const auto stateEllOpInfo = getOperatorInfo(hdsaParams, "State Elliptic Operator", true);
      Teuchos::RCP<Thyra::LinearOpBase<double>> H;
      if (!stateEllOpInfo.use_identity && !model_PME.is_null()) {        
        model_PME->block_diagonal_hessian_22(H, *rol_x_ptr, *rol_x_ptr, stateEllOpInfo.response_index, stateEllOpInfo.param_index);
      }
     
      if (!stateEllOpInfo.use_identity) {
        auto linOp = Teuchos::rcp_dynamic_cast<Thyra::LinearOpWithSolveBase<double>>(
        Teuchos::rcp_const_cast<Thyra::LinearOpBase<double>>(H));
        invEllOp = Thyra::nonconstInverse(linOp);
      }

      if(analysisVerbosity > 2)
        *out << "Piro::PerformHDSAAnalysis: End of the computation of the State Elliptic Operator" << std::endl;
    }

    {
      if(analysisVerbosity > 2)
        *out << "\nPiro::PerformHDSAAnalysis: Start the computation of the State Mass Operator" << std::endl;

      const auto stateMassOpInfo = getOperatorInfo(hdsaParams, "State Mass Operator", true);
      Teuchos::RCP<Thyra::LinearOpBase<double>> H;
      if (!stateMassOpInfo.use_identity && !model_PME.is_null()) {
        model_PME->block_diagonal_hessian_22(H, *rol_x_ptr, *rol_x_ptr, stateMassOpInfo.response_index, stateMassOpInfo.param_index);
        massOp = H;
      }

      if(analysisVerbosity > 2)
        *out << "Piro::PerformHDSAAnalysis: End of the computation of the State Mass Operator" << std::endl;
    }
  }

#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
      std::endl << "Piro::PerformHDSAAnalysis, ERROR: Teko is required for HDSA analysis"<<std::endl);
#endif

  // Run Algorithm
  int return_status = 0;

  RolOutputBuffer<char> rolOutputBuffer;
  auto rolOutputStream = Teuchos::rcp(new std::ostream(&rolOutputBuffer));
  Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(rolOutputStream);
  rolOutput->setOutputToRootOnly(0);

  ROL::Ptr<ROL::Objective_SimOpt<double> > obj_ptr = ROL::makePtr<Piro::ThyraProductME_Objective_SimOpt<double>>(model, g_index, piroParams, analysisVerbosityLevel, observer);
  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<double> > opt_prob_interface = HDSA::makePtr<HDSA::MD_ROL_Opt_Prob_Interface<double> >(obj_ptr, constr_ptr, rol_x_ptr, rol_p_ptr);
  ROL::Ptr<ROL::Vector<double> > rol_opt_p_ptr = rol_p_ptr->clone(); rol_opt_p_ptr->set(*rol_p_ptr);
  ROL::Ptr<ROL::Vector<double> > rol_opt_x_ptr = rol_x_ptr->clone(); rol_opt_x_ptr->set(*rol_x_ptr);
  HDSA::Ptr<Piro::HDSA_MD_ROL_Data_Interface<double> > data_interface = HDSA::makePtr<Piro::HDSA_MD_ROL_Data_Interface<double> >(rol_opt_x_ptr, rol_opt_p_ptr);

  double alpha_u = hdsaParams.sublist("MD Prior").get<double>("alpha_u", 0.25);
  double alpha_z = hdsaParams.sublist("MD Prior").get<double>("alpha_z", 0.0001);
  
  bool OED = hdsaParams.get("Perform HDSA Analysis With OED", false);
  auto seed_oed_init = static_cast<unsigned int>(hdsaParams.sublist("MD OED").get("Seed For Random Initialization", 0));
  ::Thyra::seed_randomize<double>( seed_oed_init );
  HDSA::Ptr<HDSA::Vector<double>> dz0;
  if (OED) {
    dz0 = data_interface->Get_z_opt()->Clone();
    dz0->Randomize_Standard_Normal();
  }

  auto hdsa_seed = static_cast<unsigned int>(hdsaParams.get("Normal Random Generator Seed", 42));
  HDSA::Ptr<HDSA::Random_Number_Generator<double> > random_number_generator = HDSA::makePtr<HDSA::Random_Number_Generator<double> >(hdsa_seed);
  ::Thyra::seed_randomize<double>(hdsa_seed+1);

  HDSA::Ptr<HDSA::MD_u_Prior_Interface<double> > u_prior_interface = HDSA::makePtr<Piro::HDSA_MD_ROL_Elliptic_u_Prior_Interface<double> >(alpha_u,random_number_generator,invEllOp,massOp);
  HDSA::Ptr<HDSA::MD_z_Prior_Interface<double> > z_prior_interface = HDSA::makePtr<Piro::HDSA_MD_ROL_Elliptic_z_Prior_Interface<double> >(alpha_z,priorEllOp,invPriorEllOp,priorMassOp,invPriorMassOp, priorMassCholOp);

  HDSA::Ptr<HDSA::Vector<double> > u_vec = data_interface->Get_u_opt()->Clone();
  HDSA::Ptr<HDSA::MD_Elliptic_u_Prior_Interface<double> > elliptic_u_prior_interface = HDSA::dynamicPtrCast<HDSA::MD_Elliptic_u_Prior_Interface<double> >(u_prior_interface);
  int num_sing_vals = hdsaParams.sublist("MD Prior").get<int>("Number Of Singular Values", 50);
  int oversampling = hdsaParams.sublist("MD Prior").get<int>("Oversampling Factor", 1);
  int num_subspace_iters = hdsaParams.sublist("MD Prior").get<int>("Number Of Subspace Iterations", 2);
  int num_prior_samples = hdsaParams.sublist("MD Prior").get("Number Of Prior Samples",100);

  elliptic_u_prior_interface->Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, *u_vec);

  HDSA::Ptr<HDSA::MD_Hessian_Analysis<double> > hessian_analysis = HDSA::makePtr<HDSA::MD_Hessian_Analysis<double> >(opt_prob_interface,z_prior_interface);
  int num_evals = hdsaParams.sublist("MD Hessian Analysis").get<int>("Rank", 20);
  oversampling = hdsaParams.sublist("MD Hessian Analysis").get<int>("Oversampling Factor", 10);
  hessian_analysis->Compute_Hessian_GEVP(data_interface->Get_z_opt(),num_evals,oversampling);
  
  int num_continuation_steps = hdsaParams.sublist("MD Continuation Update").get("Number Of Continuation Steps", 0); 
  double grad_tol = hdsaParams.sublist("MD Continuation Update").get("Gradient Tolerance", 1e-4); 
  if(!OED) {
    *out << "Performing HDSA Analysis. Sample size: " << p_samples.size() << " " << x_diff_at_samples.size() << std::endl;
    *out << "||p_opt||: " << data_interface->Get_z_opt()->Norm() << ", ||sol||: " << data_interface->Get_u_opt()->Norm() <<std::endl;
    for (int k=0; k<p_samples.size(); k++) {
      auto p_sample = createProductVector(p_samples[k]);
      ROL::Ptr<ROL::Vector<double> > rol_p_samples_ptr = ROL::makePtr<ROL::ThyraVector<double>>(p_sample);
      ROL::Ptr<ROL::Vector<double> > rol_x_diffs_ptr = ROL::makePtr<ROL::ThyraVector<double>>(x_diff_at_samples[k]);
      data_interface->Z_Data_push_back(rol_p_samples_ptr);
      data_interface->Y_Data_push_back(rol_x_diffs_ptr);
      *out << "||param_" << k << "||: " << rol_p_samples_ptr->norm() << ", ||sol_diff_: " << k << "||: " << rol_x_diffs_ptr->norm() <<std::endl;
    }

    if(num_prior_samples > 0) {
      HDSA::Ptr<HDSA::MD_Prior_Sampling<double> > prior_sampling = HDSA::makePtr<HDSA::MD_Prior_Sampling<double> >(data_interface,u_prior_interface,z_prior_interface);
      HDSA::Ptr<HDSA::MultiVector<double> > prior_samples_at_z_opt = prior_sampling->Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
    }
    
    HDSA::Ptr<HDSA::MD_Posterior_Sampling<double> > post_sampling = HDSA::makePtr<HDSA::MD_Posterior_Sampling<double> >(data_interface,u_prior_interface,z_prior_interface);
    int num_post_samples = hdsaParams.sublist("MD Posterior").get("Number Of Posterior Samples", num_prior_samples);
    double alpha_d = hdsaParams.sublist("MD Posterior").get<double>("alpha_d", 1.0e-5);
    post_sampling->Compute_Posterior_Data(alpha_d,num_post_samples);

    HDSA::Ptr<HDSA::MD_Update<double> > update = HDSA::makePtr<HDSA::MD_Update<double> >(data_interface,u_prior_interface,z_prior_interface,opt_prob_interface,post_sampling,hessian_analysis,random_number_generator,num_continuation_steps,grad_tol);

    HDSA::Ptr<HDSA::Vector<double> > z_update;
    if(num_post_samples > 0) {
      HDSA::Ptr<HDSA::MD_Posterior_Vectors<double> > posterior_update_samples = update->Posterior_Update_Samples();
      std::string name = "posterior_samples";
      posterior_update_samples->samples->Write_to_File(name);
      z_update = posterior_update_samples->mean;
    } else {
      //Alternatively just compute the mean, no sampling needed
      z_update = update->Posterior_Update_Mean(); 
    }

    HDSA::ROL_Vector<double>& z_update_rol = dynamic_cast<HDSA::ROL_Vector<double>&>(*z_update);
    rol_p_ptr->set(*z_update_rol.rol_vec);

    constr_ptr->solve(*c_ptr,*rol_x_ptr,*rol_p_ptr,tol);
  } else {

    /*
      OED setup and offline reduced matrix construction.
    */
    HDSA::Ptr<HDSA::MD_OED<double>> md_oed = HDSA::makePtr<HDSA::MD_OED<double>>(data_interface, u_prior_interface, z_prior_interface, hessian_analysis);
    md_oed->Offline_Computation();
    const int r = md_oed->Get_Reduced_Dimension();
    const typename HDSA::MD_OED<double>::Offline_Data& offline = md_oed->Get_Offline_Data();

    TEUCHOS_TEST_FOR_EXCEPTION(!offline.Is_Initialized(), std::logic_error, "Piro::PerformHDSAAnalysis, ERROR: OED offline data were not initialized." << std::endl);


    const int num_data_samples = static_cast<int>(p_samples.size());

    TEUCHOS_TEST_FOR_EXCEPTION(num_data_samples == 0, std::logic_error,
        "Piro::PerformHDSAAnalysis, ERROR: OED requires at least the initial data sample." << std::endl);

    int num_trace_samples = hdsaParams.sublist("MD OED").get<int>("Number Of Randomized Trace Samples", 50);
    const double alpha_k_denom = Estimate_Trace_Wz_Inverse_Mz(*data_interface->Get_z_opt(), *z_prior_interface, num_trace_samples);

    typename HDSA::MD_OED<double>::SPG_Options spg_options;
    spg_options.max_iter = hdsaParams.sublist("MD OED").get<int>("Max Number Of OED Iterations", 300);
    spg_options.pg_tol = hdsaParams.sublist("MD OED").get("PG Tolerance", 1e-8);
    spg_options.armijo_c = hdsaParams.sublist("MD OED").get("Armijo Coefficient", 1e-4);
    spg_options.backtrack_factor = hdsaParams.sublist("MD OED").get("Backtrack Factor", 0.5);
    spg_options.max_backtracks = hdsaParams.sublist("MD OED").get("Max Number Of Backtrack Steps", 30);
    spg_options.nonmonotone_window = hdsaParams.sublist("MD OED").get<int>("Nonmonotone Window", 5);
    spg_options.verbosity = hdsaParams.sublist("MD OED").get("Verbosity", false);

    const int num_post_samples = hdsaParams.sublist("MD Posterior").get("Number Of Posterior Samples", num_prior_samples);
    const double alpha_d = hdsaParams.sublist("MD Posterior").get<double>("alpha_d", 1.0e-5);

    /*
    * Initial guess for the OED optimization.
    */
    auto Mz_dz0 = dz0->Clone();
    z_prior_interface->Apply_M_z(*Mz_dz0, *dz0);

    const double dz0_Mz_norm = std::sqrt(dz0->Dot(*Mz_dz0));

    TEUCHOS_TEST_FOR_EXCEPTION(!(dz0_Mz_norm > 0.0), std::logic_error, "Piro::PerformHDSAAnalysis, ERROR: zero M_z norm for the OED initial direction." << std::endl);

    const double scale = 1e-2 / dz0_Mz_norm;

    HDSA::Dense_Matrix<double> beta_0(r, 1);
    beta_0.Zeros();
    for (int i = 0; i < r; ++i) {
      const double beta_i =  scale * (*hessian_analysis->Get_Evecs())[i]->Dot(*Mz_dz0);
      beta_0.Set_Entry(i, 0, beta_i);
    }

    /*
    * betas contains the reduced coordinates of the previously proposed OED parameters: 
    *   p_samples[1], p_samples[2], ...
    *   p_samples[0] is z_opt and corresponds to beta = 0, which is included implicitly by the OED formulation.
    */
    HDSA::Ptr<HDSA::Dense_Matrix<double>> betas = HDSA::makePtr<HDSA::Dense_Matrix<double>>(0, 1);

    std::vector<HDSA::Ptr<HDSA::Vector<double>>> z_bars;

    HDSA::Ptr<HDSA::Vector<double>> z_lofi = data_interface->Get_z_opt()->Clone();
    z_lofi->Set(*data_interface->Get_z_opt());

    HDSA::Ptr<HDSA::Dense_Matrix<double>> current_beta_bar =  HDSA::nullPtr;

    *out << "\n=====================================================" << std::endl;
    *out << "Beginning sequential OED workflow" << std::endl;
    *out << "Reduced dimension r = " << r << std::endl;
    *out << "Number of available data samples = " << num_data_samples << std::endl;
    *out << "====================================================="  << std::endl;

    HDSA::Ptr<HDSA::Vector<double>> hdsa_rol_p_ptr = HDSA::makePtr<HDSA::ROL_Vector<double>>(rol_p_ptr);

    /*
    * Replay only the posterior updates using the available data.
    */
    for (int step = 0; step < num_data_samples; ++step) {

      auto p_sample = createProductVector(p_samples[step]);

      ROL::Ptr<ROL::Vector<double>> rol_p_samples_ptr = ROL::makePtr<ROL::ThyraVector<double>>(p_sample);

      ROL::Ptr<ROL::Vector<double>> rol_x_diffs_ptr = ROL::makePtr<ROL::ThyraVector<double>>(x_diff_at_samples[step]);

      data_interface->Z_Data_push_back(rol_p_samples_ptr);
      data_interface->Y_Data_push_back(rol_x_diffs_ptr);

      /*
      * Recover beta for a previously proposed OED parameter.
      *
      * sample 0 is z_opt, so beta_0_sample = 0.
      */
      if (step > 0) {

        HDSA::Ptr<HDSA::Vector<double>> z_sample =  HDSA::makePtr<HDSA::ROL_Vector<double>>(rol_p_samples_ptr);

        HDSA::Ptr<HDSA::Dense_Matrix<double>> beta_sample = Compute_Beta_From_z(*z_sample, *data_interface->Get_z_opt(), *z_prior_interface, offline);

        betas = Append_Betas(*betas, *beta_sample);

        /*
        * Optional diagnostic: check how well the current reduced
        * space represents the actual previously proposed parameter.
        */
        {
          HDSA::Ptr<HDSA::Vector<double>> z_reconstructed = data_interface->Get_z_opt()->Clone();

          z_reconstructed->Set(*data_interface->Get_z_opt());

          for (int i = 0; i < r; ++i) {
            z_reconstructed->Scaled_Plus((*beta_sample)(i, 0), *(*offline.V)[i]);
          }

          const double projection_error =  M_z_Norm_Difference(*z_sample, *z_reconstructed, *z_prior_interface);

          *out << "Reduced-space projection error for sample "  << step << " = "  << projection_error << std::endl;
        }
      }

      /*
      * Recompute the discrepancy posterior using all samples available through this step.
      */
      HDSA::Ptr<HDSA::MD_Posterior_Sampling<double>> post_sampling =
          HDSA::makePtr<HDSA::MD_Posterior_Sampling<double>>(data_interface, u_prior_interface, z_prior_interface);

      post_sampling->Compute_Posterior_Data(alpha_d, num_post_samples);

      HDSA::Ptr<HDSA::Vector<double>> u_k = data_interface->Get_u_opt()->Clone();
      HDSA::Ptr<HDSA::Vector<double>> z_k = data_interface->Get_z_opt()->Clone();
      HDSA::Ptr<HDSA::Vector<double>> beta_k = HDSA::makePtr<HDSA::Std_Vector<double>>(r);

      HDSA::Ptr<HDSA::MD_Continuation_Update<double>> cont_update = HDSA::makePtr<HDSA::MD_Continuation_Update<double>>(
        data_interface, z_prior_interface, opt_prob_interface, post_sampling, hessian_analysis, random_number_generator, num_continuation_steps, grad_tol);

      *out << "\nPosterior update using "  << step + 1  << " data sample(s)"  << std::endl;

      cont_update->Posterior_Update_Mean(*u_k, *z_k, *beta_k);
      current_beta_bar = Vector_To_Dense(*beta_k);
      z_bars.push_back(z_k);
    }

    /*
    * Determine the next OED radius and alpha_k.
    */
    HDSA::Ptr<HDSA::Vector<double>> radius_reference;

    if (num_data_samples == 1) {
      radius_reference = z_lofi;
    } else {
      radius_reference = z_bars[num_data_samples - 2];
    }

    const double prev_z_distance = M_z_Norm_Difference(*z_bars.back(), *radius_reference, *z_prior_interface);
    const double alpha_k = prev_z_distance * prev_z_distance / alpha_k_denom;
    const double constr_radius = prev_z_distance;

    *out << "\nNext OED proposal" << std::endl;
    *out << "-----------------------------------------------------" << std::endl;
    *out << "Previous posterior movement radius = " << constr_radius << std::endl;
    *out << "OED covariance coefficient alpha_k = " << alpha_k << std::endl;

    md_oed->Set_Covariance_Coefficient(alpha_k);
    
    typename HDSA::MD_OED<double>::Seq_Design_Result seq_result =
        md_oed->Generate_Seq_Optimal_Design(beta_0, alpha_d, *betas, *current_beta_bar, constr_radius, spg_options);

    *out << "Sequential OED final objective = " << seq_result.optimizer_info.final_objective  << std::endl;
    *out << "Sequential OED projected-gradient norm = " << seq_result.optimizer_info.projected_gradient_norm  << std::endl;

    /*
    * Return the newly proposed parameter.
    */
    hdsa_rol_p_ptr->Set(*(*seq_result.Z_new)[0]);
  } 

  if(Teuchos::nonnull(observer)) {
    const ROL::ThyraVector<double>  & thyra_x = dynamic_cast<const ROL::ThyraVector<double>&>(*rol_x_ptr);
    observer->parametersChanged();
    observer->observeSolution(1, *(thyra_x.getVector()), Teuchos::null, Teuchos::null, Teuchos::null);
  }

  return 0;

#else
  (void)piroModel;
  (void)p;
  out = Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "ERROR: Trilinos/Piro was not configured to include HDSA analysis."
       << "\nYou must enable HDSA and ROL." << endl;
  return 0;  // should not fail tests
#endif
}


RCP<const Teuchos::ParameterList>
Piro::getValidPiroAnalysisParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("Valid Piro Analysis Params"));;

  validPL->set<std::string>("Analysis Package", "","Must be: Solve or ROL.");
  validPL->set<bool>("Output Final Parameters", false, "");
  validPL->set<bool>("Transient", false, "");
  validPL->sublist("Solve",     false, "");
  validPL->sublist("ROL",       false, "");
  validPL->sublist("HDSA",       false, "");
  validPL->set<int>("Output Level", 2, "Verbosity level, ranges from 0 (no output) to 4 (extreme output)");
  validPL->set<int>("Write Interval", 1, "Iterval between writes to mesh");

  return validPL;
}


RCP<const Teuchos::ParameterList>
Piro::getValidPiroAnalysisROLParameters(int num_parameters)
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("Valid Piro Analysis ROL Params"));

  validPL->set<int>("Response Vector Index", 0, "Index of the response to be used as objective for ROL optimization");
  validPL->set<int>("Number Of Parameters", 1, "Number of the parameters to use as control for ROL optimization");
  
  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    validPL->set<int>(ss.str(), 0, "Index Of the Parameter Vector to be used as control for ROL optimization");
  }

  validPL->set<std::string>("Parameter Initial Guess Type", "From Model Evaluator", "How to initialize parameters, options: \"Uniform Vector\", \"Random Vector\", \"From Model Evaluator\" (use the value stored in model evaluator)");
  validPL->set<double>("Uniform Parameter Guess", 2.0, "Value to use to uniformly intialize the parameter");
  Teuchos::Array<double> range = Teuchos::tuple(1.0,3.0);
  validPL->set<Teuchos::Array<double>>("Min And Max Of Random Parameter Guess", range, "Array providing the range of values of values to randomply intialize the parameter");  
  validPL->set<int>("Seed For Thyra Randomize", 42, "Seed of Thyra random generator");

  validPL->set<bool>("Check Derivatives", false, "Whether to perform derivatives check");
  validPL->set<bool>("Perform Optimization", false, "Whether to perform the optimization (this allows to check the derivative without solving the optimization problem)");
  validPL->set<bool>("Test Vector", false, "Whether to check the implmentation of ROL Thyra Vector");
  validPL->set<int>("Number Of Vector Tests", 1, "Number of vectors to use when testing the implmentation of ROL Thyra Vector");

  validPL->set<bool>("Bound Constrained", true, "Whether to enforce bounds to the parameters during the optimization");
  validPL->set<bool>("Full Space", true, "Whether to use a full-space or a reduced-space optimization approach");
  validPL->set<bool>("Tempus Driver", false, "Whether to use Tempus to compute the derivative");
  validPL->set<bool>("Redirect Tempus Output", true, "Whether to redirect Tempus output to file");
  validPL->set<string>("Tempus Output Filename", "log_tempus.txt", "Filename for the Tempus output");
  validPL->set<bool>("Response Depends Only On Final Time", true, "Whether the response depends only on the solution and parameters at the final time");
  validPL->set<bool>("Use NOX Solver", true, "Whether to use NOX for solving the state equation or the native ROL solver");

  validPL->set<double>("Objective Recovery Value", 1.0e10, "Objective value used when the state solver does not converge. If not defined, the objective will be computed using the unconverged state");

  validPL->sublist("Derivative Checks",  false, "Options for derivative checks");
  validPL->sublist("ROL Options",  false, "Options to pass to ROL");
  validPL->sublist("Matrix Based Dot Product",  false, "Sublist to define a Matrix based dot product (instead of the l2 one) to define gradient in ROL");
  validPL->sublist("Custom Secant", false, "Sublist to define a custom secant");

  return validPL;
}
