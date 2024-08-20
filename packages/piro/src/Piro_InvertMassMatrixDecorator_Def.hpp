// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>


#include "Piro_InvertMassMatrixDecorator.hpp"
#include "Thyra_VectorBase.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
//#include "Thyra_EpetraThyraWrappers.hpp"
//#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"

#ifdef HAVE_PIRO_MUELU
#include "Teuchos_AbstractFactoryStd.hpp"
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#endif


template <typename Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::InvertMassMatrixDecorator(
                          Teuchos::RCP<Teuchos::ParameterList> stratParams,
                          Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model_,
                          bool massMatrixIsConstant_, bool lumpMassMatrix_,
                          bool massMatrixIsCoeffOfSecondDeriv_) :
  model(model_),
  massMatrixIsConstant(massMatrixIsConstant_),
  lumpMassMatrix(lumpMassMatrix_),
  massMatrixIsCoeffOfSecondDeriv(massMatrixIsCoeffOfSecondDeriv_),
  calcMassMatrix(true)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Create x_dot vector, fill with 0.0 so implicit fill gives
  // correct results for explicit fill
  x_dot = Thyra::createMember<Scalar>(model->get_x_space());
  Thyra::put_scalar<Scalar>(0.0, x_dot.ptr());

  // get allocated space for Mass Matrix
  massMatrix = model->create_W_op();
  if (lumpMassMatrix) invDiag = Thyra::createMember<Scalar>(model->get_x_space());

  Teuchos::RCP<Teuchos::FancyOStream> out
     = Teuchos::VerboseObjectBase::getDefaultOStream();

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef HAVE_PIRO_MUELU
  Stratimikos::enableMueLu(linearSolverBuilder);
#endif

   if (stratParams->isParameter("Linear Solver Type")) {
     const std::string lin_solver_type =  stratParams->get<std::string>("Linear Solver Type");
     if (lin_solver_type == "Amesos") {
       *out << "WARNING: Amesos solver does not work with Piro::InvertMassMatrix; switching to Amesos2 to avoid cast error.\n";
       stratParams->set<std::string>("Linear Solver Type", "Amesos2");  
     }
     else if (lin_solver_type == "AztecOO") {
       *out << "WARNING: AztecOO solver does not work with Piro::InvertMassMatrix; switching to Belos to avoid cast error.\n";
       stratParams->set<std::string>("Linear Solver Type", "Belos"); 
       if (stratParams->get<std::string>("Preconditioner Type") == "Ifpack") { 
         stratParams->set<std::string>("Preconditioner Type", "Ifpack2");  
         *out << "WARNING: Ifpack preconditioner does not work with Piro::InvertMassMatrix; switching to Ifpack2 to avoid cast error.\n";
       }
       else if (stratParams->get<std::string>("Preconditioner Type") ==  "ML") {
         *out << "WARNING: ML preconditioner does not work with Piro::InvertMassMatrix; switching to MueLu to avoid cast error.\n";
         stratParams->set<std::string>("Preconditioner Type", "MueLu");  
       }
     }
   }
   linearSolverBuilder.setParameterList(stratParams);

  // Create a linear solver factory given information read from the
  // parameter list.
  //lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");
  lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

  // Setup output stream and the verbosity level
  lowsFactory->setOStream(out);
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);
}

template<typename Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::~InvertMassMatrixDecorator()
{
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_x_space() const
{
  return model->get_x_space();
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_f_space() const
{
  return model->get_f_space();
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_f_multiplier_space() const
{
  return model->get_f_multiplier_space();
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_p_space(int l) const
{
  return model->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
Piro::InvertMassMatrixDecorator<Scalar>::get_p_names(int l) const
{
  return model->get_p_names(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_g_space(int j) const
{
  return model->get_g_space(j);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_g_multiplier_space(int j) const
{
  return model->get_g_multiplier_space(j);
}

template<typename Scalar>
Teuchos::ArrayView<const std::string>
Piro::InvertMassMatrixDecorator<Scalar>::get_g_names(int j) const
{
  return model->get_g_names(j);
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues = this->createInArgsImpl();
  nominalValues.setArgs(
      model->getNominalValues(),
      /* ignoreUnsupported = */ true,
      /* cloneObjects = */ false);
  return nominalValues;
}

template<typename Scalar>
Teuchos::RCP< Thyra::LinearOpBase< Scalar > >
Piro::InvertMassMatrixDecorator<Scalar>::create_W_op () const
{
  return model->create_W_op();
}

template<typename Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > 
Piro::InvertMassMatrixDecorator<Scalar>::create_W () const
{
  Teuchos::RCP<const Thyra::LinearOpBase<double> > A = massMatrix; 
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > W = Thyra::linearOpWithSolve(*lowsFactory, A);
  return W; 
}

template<typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_W_factory() const
{
  return model->get_W_factory();
}

template<typename Scalar>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::create_W_prec() const
{
  return model->create_W_prec();
}


template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::getLowerBounds() const
{
  return Thyra::ModelEvaluatorBase::InArgs<Scalar>(); // Default value
}


template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::getUpperBounds() const
{
  return Thyra::ModelEvaluatorBase::InArgs<Scalar>(); // Default value
}

template<typename Scalar>
void
Piro::InvertMassMatrixDecorator<Scalar>::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint,
    const bool wasSolved)
{
  model->reportFinalPoint(finalPoint,wasSolved);
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::createInArgs() const
{
  return this->createInArgsImpl();
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> result = model->createInArgs();
  result.setModelEvalDescription(this->description());
  return result; 
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> result = model->createOutArgs();
  result.setModelEvalDescription(this->description());
  return result; 
}

template<typename Scalar>
void
Piro::InvertMassMatrixDecorator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //IKT 9/22/2021: to circumvent the following, need to implement M*DfDp in the code.
  //Talk to Eric Phipps if you have questions about this.
  if (outArgs.Np()>0) {
    if (outArgs.get_DfDp(0).getMultiVector() != Teuchos::null)
      TEUCHOS_TEST_FOR_EXCEPTION(true,
          Teuchos::Exceptions::InvalidParameter,
          "\n Error! Piro::InvertMassMatrixDecorator, needed for explicit time-stepping, does not have DfDp implemented!\n" <<
          "Forward sensitivities will not work; please re-run without them or run with an implicit time-stepper.\n";)
  }

  if (outArgs.get_f() == Teuchos::null) {
    // Probably just getting g -- pass through
    model->evalModel(inArgs, outArgs);
  }
  else {
    // Add massMatrix to outargs, with appropriate alpha and beta
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs(outArgs);
    Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs(inArgs);

    if (!massMatrixIsCoeffOfSecondDeriv) {
      modelInArgs.set_x_dot(x_dot);
      modelInArgs.set_alpha(-1.0);
      modelInArgs.set_beta(0.0);
    }
    else {  // Mass Matrix is coeff of Second deriv
      modelInArgs.set_x_dot_dot(x_dot);
      modelInArgs.set_alpha(0.0);
      modelInArgs.set_beta(0.0);
      modelInArgs.set_W_x_dot_dot_coeff(-1.0);
    }

    if (calcMassMatrix) {
      modelOutArgs.set_W_op(massMatrix);
      //The following 2 lines were added to prevent Jacobian and residual 
      //from being set at the same time...  see below.
      modelOutArgs.set_f(Teuchos::null); 
      model->evalModel(modelInArgs, modelOutArgs); 
    }

    //Create mass matrix 
    if (calcMassMatrix) {
      if (!lumpMassMatrix) {
        // Create a linear solver based on the forward operator massMatrix
	lows = this->create_W(); 
      }
      else { // Lump matrix into inverse of diagonal
        Thyra::put_scalar<Scalar>(1.0, invDiag.ptr());
        Thyra::apply<Scalar>(*massMatrix, Thyra::NOTRANS, *invDiag, invDiag.ptr(), 1.0, 0.0);
        Thyra::reciprocal<Scalar>(*invDiag, invDiag.ptr());
        //IKT, 5/31/17: adding the following logic which checks invDiag vector for nans
        //and throws an exception if nans are found.
        typedef Teuchos::ScalarTraits<typename Thyra::ModelEvaluator<Scalar>::ScalarMag> SMT;
        const Scalar sumInvDiag = sum(*invDiag);
        bool isNanInvDiag = SMT::isnaninf(sumInvDiag);
        if (isNanInvDiag) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,
              Teuchos::Exceptions::InvalidParameter,
              "\n Error! Piro::InvertMassMatrixDecorator: lumped mass has 0 values, \n" <<
              "so its inverse is not defined! \n";)
        } 
      }
    }

    //set f and W_op in modelOutArgs
    modelOutArgs.set_f(outArgs.get_f()); 
    modelOutArgs.set_W_op(outArgs.get_W());
    
    //Evaluate the underlying model
    model->evalModel(modelInArgs, modelOutArgs);
    
    // Invert the mass matrix:   f_exp = M^{-1} f_imp
    if (!lumpMassMatrix) { //standard solve
      // Solve the linear system for x, given b
      ::Thyra::solve<double>(*lows, ::Thyra::NOTRANS, *modelOutArgs.get_f(), outArgs.get_f().ptr());
    }
    else { //diagonal lumped solve
      Teuchos::RCP<Thyra::VectorBase<Scalar> > f = outArgs.get_f();
      Thyra::ele_wise_prod_update<Scalar>(1.0, *invDiag, f.ptr());
    }

    // Do not recompute mass matrix in future if it is a constant
    if (massMatrixIsConstant) calcMassMatrix = false;
  }
}
