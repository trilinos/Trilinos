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
#include "Piro_CustomLBFGSSecant.hpp"
#include "ROL_LinearOpScaledThyraVector.hpp"
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

int
Piro::PerformAnalysis(
    Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
    Teuchos::ParameterList& piroParams,
    RCP< Thyra::VectorBase<double> >& result,
    RCP< Piro::ROL_ObserverBase<double> > observer)
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

  }
#endif
  else {
    if (analysis == "ROL")
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
  std::vector<std::string> p_names;

  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    const auto names_array = *piroSSSolver->getModel().get_p_names(0);
    for (int k=0; k<names_array.size(); k++) {
      p_names.push_back(names_array[k]);
    }
  }

  //set names of parameters in the "Optimization Status" sublist
  piroParams.sublist("Optimization Status").set("Parameter Names", Teuchos::rcpFromRef(p_names));

  if(rolParams.isParameter("Objective Recovery Value"))
    piroParams.sublist("Optimization Status").set("Objective Recovery Value", rolParams.get<double>("Objective Recovery Value"));

  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space = model->get_p_space(0);
  p = model->getNominalValues().get_p(0)->clone_v();

  ROL::ThyraVector<double> rol_p(p);
  //Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space;
  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> x_space = model->get_x_space();

  Teuchos::RCP<Thyra::VectorBase<double>> x = Thyra::createMember(x_space);
  Thyra::copy(*model->getNominalValues().get_x(), x.ptr());

  ROL::ThyraVector<double> rol_x(x);
  Teuchos::RCP<Thyra::VectorBase<double>> lambda_vec = Thyra::createMember(x_space);
  ROL::ThyraVector<double> rol_lambda(lambda_vec);

  Teuchos::EVerbosityLevel analysisVerbosityLevel;
  switch(analysisVerbosity) {
    case 1: analysisVerbosityLevel= Teuchos::VERB_LOW; break;
    case 2: analysisVerbosityLevel= Teuchos::VERB_MEDIUM; break;
    case 3: analysisVerbosityLevel= Teuchos::VERB_HIGH; break;
    case 4: analysisVerbosityLevel= Teuchos::VERB_EXTREME; break;
    default: analysisVerbosityLevel= Teuchos::VERB_NONE;
  }  
  Piro::ThyraProductME_Objective_SimOpt<double> obj(model, g_index, piroParams, analysisVerbosityLevel, observer);
  Piro::ThyraProductME_Constraint_SimOpt<double> constr(model, adjointModel, piroParams, analysisVerbosityLevel, observer);

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

    auto rol_x_zero = rol_x.clone(); rol_x_zero->zero();
    auto rol_p_zero = rol_p.clone(); rol_p_zero->zero();

    int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
    double norm_p = rol_p.norm();
    double norm_x = rol_x.norm();

    ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x),ROL::makePtrFromRef(rol_p));

    for(int i=0; i< num_checks; i++) {

      *out << "\nPiro::PerformROLSteadyAnalysis: Performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

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

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Reduced Derivative Checks", false)) {
        *out << "Piro::PerformROLSteadyAnalysis: Checking Reduced Gradient Accuracy" << std::endl;
        reduced_obj.checkGradient(rol_p, rol_p_direction1, true, *out);
      }
      // Check derivatives.

      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Objective Gradient " << std::endl;
      obj.checkGradient(sopt_vec,sopt_vec_direction1,true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Objective Gradient in x direction" << std::endl;
      obj.checkGradient(sopt_vec,sopt_vec_direction1_x,true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Objective Gradient in p direction" << std::endl;
      obj.checkGradient(sopt_vec,sopt_vec_direction1_p,true,*out,num_steps,order);


      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Constraint Gradient " << std::endl;
      constr.checkApplyJacobian(sopt_vec,sopt_vec_direction1,rol_x_direction1, true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Constraint Gradient in x direction (Jacobian) " << std::endl;
      constr.checkApplyJacobian(sopt_vec,sopt_vec_direction1_x,rol_x_direction1,true,*out,num_steps,order);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of Constraint Gradient in p direction" << std::endl;
      constr.checkApplyJacobian(sopt_vec,sopt_vec_direction1_p,rol_x_direction1,true,*out,num_steps,order);

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Expensive Derivative Checks", false))
        constr.checkApplyAdjointJacobian(sopt_vec,rol_x_direction1,rol_x_direction1,sopt_vec,true,*out,num_steps);

      *out << "Piro::PerformROLSteadyAnalysis: Checking Consistency of Constraint Gradient and its adjoint" << std::endl;
      constr.checkAdjointConsistencyJacobian(rol_x_direction1, sopt_vec_direction2, sopt_vec,true,*out);

      obj.update(rol_x,rol_p,ROL::UpdateType::Temp);
      constr.update(rol_x,rol_p,ROL::UpdateType::Temp);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian" << std::endl;
      obj.checkHessSym(sopt_vec,sopt_vec_direction1, sopt_vec_direction2, true,*out);

      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian (H_xx = H_xx^T)" << std::endl;
      obj.checkHessSym(sopt_vec,sopt_vec_direction1_x, sopt_vec_direction2_x, true,*out);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian (H_xp = H_px^T)" << std::endl;
      obj.checkHessSym(sopt_vec,sopt_vec_direction1_x, sopt_vec_direction2_p, true,*out);
      *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of objective Hessian (H_pp = H_pp^T)" << std::endl;
      obj.checkHessSym(sopt_vec,sopt_vec_direction1_p, sopt_vec_direction2_p, true,*out);

      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of objective Hessian" << std::endl;
      obj.checkHessVec(sopt_vec,sopt_vec_direction1,true,*out,num_steps,order);

      if(rolParams.sublist("Derivative Checks").get<bool>("Perform Reduced Derivative Checks", false)) {
        *out << "Piro::PerformROLSteadyAnalysis: Checking Symmetry of reduced objective Hessian" << std::endl;
        reduced_obj.update(rol_p,ROL::UpdateType::Temp);
        auto hsymCheck = reduced_obj.checkHessSym(rol_p, rol_p_direction1, rol_p_direction2, false,*out);
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
        reduced_obj.checkHessVec(rol_p, rol_p_direction1,true,*out,num_steps,order);
      }

      *out << "Piro::PerformROLSteadyAnalysis: Checking Accuracy of constraint Hessian" << std::endl;
      constr.checkApplyAdjointHessian(sopt_vec, rol_x_direction1, sopt_vec_direction2, sopt_vec_direction2, true,*out,num_steps,order);

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
  int reponse_index_dotProd = -1;
  if(rolParams.isSublist("Matrix Based Dot Product")) {
    const Teuchos::ParameterList& matrixDotProductList = rolParams.sublist("Matrix Based Dot Product");
    auto matrixType = matrixDotProductList.get<std::string>("Matrix Type");
    if(matrixType == "Hessian Of Response") {
      useCustomDotProduct = true;
      reponse_index_dotProd = matrixDotProductList.sublist("Matrix Types").sublist("Hessian Of Response").get<int>("Response Index");
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
  int reponse_index_secant = -1;
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
      reponse_index_secant = customSecantList.sublist("Initialization Types").sublist("Hessian Of Response").get<int>("Response Index");
    }
    else if(initializationType == "Identity") {
      reponse_index_secant = -1;
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
      model_PME->block_diagonal_hessian_22(bH_dotP, rol_x, rol_p, reponse_index_dotProd);
    }
    if(useCustomSecant && (reponse_index_secant != -1 ) && !model_PME.is_null()) {

      if (reponse_index_dotProd == reponse_index_secant)
        bH_sec = bH_dotP;
      else {
        bH_sec = Teko::createBlockedOp();
        model_PME->block_diagonal_hessian_22(bH_sec, rol_x, rol_p, reponse_index_secant);
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
      if(reponse_index_secant == -1 ) {// identity initialization 
        invH_sec = H_sec = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(p_space));
      } else if ((reponse_index_dotProd == reponse_index_secant) && Teuchos::nonnull(H_dotP) && Teuchos::nonnull(invH_dotP)) {
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
  (void)reponse_index_dotProd;
  (void)reponse_index_secant;
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
      std::ostream rolOutputStream(&rolOutputBuffer);
      Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(Teuchos::rcpFromRef(rolOutputStream));
      rolOutput->setOutputToRootOnly(0);

      

      if ( useFullSpace ) {
        //using default dot product for x
        ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x),rol_p_primal);
        auto r_ptr = rol_x.clone();
        double tol = 1e-5;
        constr.solve(*r_ptr,rol_x,rol_p,tol);
        if(boundConstrained) {
          *out << "Piro::PerformROLSteadyAnalysis: Solving Full Space Bound Constrained Optimization Problem" << std::endl;
          ROL::BoundConstraint<double> u_bnd(rol_x);
          ROL::Ptr<ROL::BoundConstraint<double> > bnd = ROL::makePtr<ROL::BoundConstraint_SimOpt<double> >(ROL::makePtrFromRef(u_bnd),boundConstraint);
          ROL::Problem<double> prob(ROL::makePtrFromRef(obj), ROL::makePtrFromRef(sopt_vec));
          prob.addBoundConstraint(bnd);
          prob.addConstraint("Constraint", ROL::makePtrFromRef(constr),r_ptr);
          bool lumpConstraints(false), printToStream(true);
          prob.finalize(lumpConstraints, printToStream, *rolOutput);
          ROL::Solver<double> optSolver(ROL::makePtrFromRef(prob), rolParams.sublist("ROL Options"));
          optSolver.solve(*out);
          return_status = optSolver.getAlgorithmState()->statusFlag;
        } else {
          *out << "Piro::PerformROLSteadyAnalysis: Solving Full Space Unconstrained Optimization Problem" << std::endl;
          ROL::Problem<double> prob(ROL::makePtrFromRef(obj), ROL::makePtrFromRef(sopt_vec));//, ROL::makePtrFromRef(constr), r_ptr);
          prob.addConstraint("Constraint", ROL::makePtrFromRef(constr),r_ptr);
          bool lumpConstraints(false), printToStream(true);
          prob.finalize(lumpConstraints, printToStream, *rolOutput);
          ROL::Solver<double> optSolver(ROL::makePtrFromRef(prob), rolParams.sublist("ROL Options"));
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
  std::vector<std::string> p_names;

  for(int i=0; i<num_parameters; ++i) {
    std::ostringstream ss; ss << "Parameter Vector Index " << i;
    const auto names_array = *piroTSolver->getModel().get_p_names(0);
    for (int k=0; k<names_array.size(); k++) {
      p_names.push_back(names_array[k]);
    }
  }

  //set names of parameters in the "Optimization Status" sublist
  piroParams.sublist("Optimization Status").set("Parameter Names", Teuchos::rcpFromRef(p_names));

  if(rolParams.isParameter("Objective Recovery Value"))
    piroParams.sublist("Optimization Status").set("Objective Recovery Value", rolParams.get<double>("Objective Recovery Value"));

  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space = model->get_p_space(0);
  p = model->getNominalValues().get_p(0)->clone_v();

  ROL::ThyraVector<double> rol_p(p);
  //Teuchos::RCP<Thyra::VectorSpaceBase<double> const> p_space;
  Teuchos::RCP<Thyra::VectorSpaceBase<double> const> x_space = model->get_x_space();

  Teuchos::RCP<Thyra::VectorBase<double>> x = Thyra::createMember(x_space);
  Thyra::copy(*model->getNominalValues().get_x(), x.ptr());

  ROL::ThyraVector<double> rol_x(x);
  Teuchos::RCP<Thyra::VectorBase<double>> lambda_vec = Thyra::createMember(x_space);
  ROL::ThyraVector<double> rol_lambda(lambda_vec);

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

  ROL::Ptr<ROL::Vector<double> > rol_p_ptr = ROL::makePtrFromRef(rol_p);
  ROL::Ptr<ROL::Vector<double> > rol_x_ptr = ROL::makePtrFromRef(rol_x);
  ROL::Ptr<ROL::Vector<double> > rol_lambda_ptr = ROL::makePtrFromRef(rol_lambda);

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
  std::ostream rolOutputStream(&rolOutputBuffer);
  Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(Teuchos::rcpFromRef(rolOutputStream));
  rolOutput->setOutputToRootOnly(0);

  Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_p = p->clone_v();
  ::Thyra::put_scalar<double>( 1.0, scaling_vector_p.ptr());
  ROL::PrimalScaledThyraVector<double> rol_p_primal(p, scaling_vector_p);

  if(useTempusDriver) {

    Piro::ThyraProductME_TempusFinalObjective<double> tempus_obj(model, piroTempusSolver, forward_integrator, adjoint_integrator, adjointModel, g_index, piroParams, nt, analysisVerbosityLevel, observer);

    ROL::Ptr<ROL::Objective<double> > obj_ptr = ROL::makePtrFromRef(tempus_obj);

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

        std::streambuf *coutbuf;
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
        
        std::streambuf *coutbuf;
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

      auto rol_x_zero = rol_x.clone(); rol_x_zero->zero();
      auto rol_p_zero = rol_p.clone(); rol_p_zero->zero();

      int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
      double norm_p = rol_p.norm();
      double norm_x = rol_x.norm();

      ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x),ROL::makePtrFromRef(rol_p));

      for(int i=0; i< num_checks; i++) {

        *out << "\nPiro::PerformROLTransientAnalysis: Performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

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

        *out << "Piro::PerformROLTransientAnalysis: Checking Reduced Gradient Accuracy" << std::endl;

        RolOutputBuffer<char> rolOutputBuffer;
        std::ostream rolOutputStream(&rolOutputBuffer);
        Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(Teuchos::rcpFromRef(rolOutputStream));
        rolOutput->setOutputToRootOnly(0);

        const double ten(10);
        std::vector<double> steps(ROL_NUM_CHECKDERIV_STEPS);
        for(int i=0;i<ROL_NUM_CHECKDERIV_STEPS;++i) {
          steps[i] = pow(ten,static_cast<double>(-i-1));
        }

        std::streambuf *coutbuf;
        std::ofstream out_file;
        if(rolParams.get<bool>("Redirect Tempus Output", true)) {
          out_file.open(rolParams.get<string>("Tempus Output Filename", "log_tempus.txt"));
          coutbuf = std::cout.rdbuf();
          std::cout.rdbuf(out_file.rdbuf());
        }
        obj_ptr->checkGradient(rol_p_primal, rol_p_primal.dual(), rol_p_direction1, steps, true, *rolOutput, 1);
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
    Piro::ThyraProductME_ROL_DynamicObjective<double> obj(model, forward_integrator, adjoint_integrator, adjointModel, g_index, piroParams, nt, useFinalTimeStepResponse, true, analysisVerbosityLevel, observer);
    Piro::ThyraProductME_ROL_DynamicConstraint<double> constr(forward_integrator, adjoint_integrator, adjointModel, piroParams, analysisVerbosityLevel, observer);

    constr.setSolveParameters(rolParams.sublist("ROL Options"));
    constr.setNumResponses(piroTSolver->num_g());

    ROL::Ptr<ROL::DynamicObjective<double> > obj_ptr = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::DynamicConstraint<double> > constr_ptr = ROL::makePtrFromRef(constr);

    if ( useFullSpace ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Piro::PerformROLTransientAnalysis, ERROR: " <<
          "full space approach is currently not supported."<<std::endl);
    }
    else {
      ROL::ReducedDynamicObjective<double> reduced_obj(obj_ptr,constr_ptr,rol_x_ptr,rol_p_ptr,rol_lambda_ptr, *timeStamps, piroParams);

      ROL::Ptr<ROL::ReducedDynamicObjective<double> > reduced_obj_ptr = ROL::makePtrFromRef(reduced_obj);
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

        auto rol_x_zero = rol_x.clone(); rol_x_zero->zero();
        auto rol_p_zero = rol_p.clone(); rol_p_zero->zero();

        int num_checks = rolParams.sublist("Derivative Checks").get<int>("Number Of Derivative Checks", 1);
        double norm_p = rol_p.norm();
        double norm_x = rol_x.norm();

        ROL::Vector_SimOpt<double> sopt_vec(ROL::makePtrFromRef(rol_x),ROL::makePtrFromRef(rol_p));

        for(int i=0; i< num_checks; i++) {

          *out << "\nPiro::PerformROLTransientAnalysis: Performing gradient check " << i+1 << " of " << num_checks << ", at parameter initial guess" << std::endl;

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

          *out << "Piro::PerformROLTransientAnalysis: Checking Reduced Gradient Accuracy" << std::endl;

          RolOutputBuffer<char> rolOutputBuffer;
          std::ostream rolOutputStream(&rolOutputBuffer);
          Teuchos::RCP<Teuchos::FancyOStream> rolOutput = Teuchos::getFancyOStream(Teuchos::rcpFromRef(rolOutputStream));
          rolOutput->setOutputToRootOnly(0);

          const double ten(10);
          std::vector<double> steps(ROL_NUM_CHECKDERIV_STEPS);
          for(int i=0;i<ROL_NUM_CHECKDERIV_STEPS;++i) {
            steps[i] = pow(ten,static_cast<double>(-i-1));
          }
          reduced_stationarycontrols_obj.checkGradient(rol_p_primal, rol_p_primal.dual(), rol_p_direction1, steps, true, *rolOutput, 1);

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
