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

#include <cmath>


#include "Piro_InvertMassMatrixDecorator.hpp"
#include "Thyra_VectorBase.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
//#include "Thyra_EpetraThyraWrappers.hpp"
//#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"

#ifdef Piro_ENABLE_Ifpack2
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#endif

#ifdef Piro_ENABLE_MueLu
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueluTpetraHelpers.hpp"
#endif


template<typename Scalar>
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

#ifdef Piro_ENABLE_Ifpack2
  typedef Thyra::PreconditionerFactoryBase<double> Base;
  typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double> > Impl;
  linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef Piro_ENABLE_MueLu
  Stratimikos::enableMueLuTpetra(linearSolverBuilder);
#endif
  
   linearSolverBuilder.setParameterList(stratParams);

  // Create a linear solver factory given information read from the
  // parameter list.
  lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

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
  // TODO
  TEUCHOS_TEST_FOR_EXCEPTION(true,
         Teuchos::Exceptions::InvalidParameter,
         "Calling reportFinalPoint in Piro_InvertMassMatrixDecorator_Def.hpp line 215" << std::endl);
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
  return model->createInArgs();
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> 
Piro::InvertMassMatrixDecorator<Scalar>::createOutArgsImpl() const 
{
  return model->createOutArgs();  
}

template<typename Scalar>
void Piro::InvertMassMatrixDecorator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  if (outArgs.Np()>0) {
   if (outArgs.get_DfDp(0).getMultiVector() != Teuchos::null) 
     std::cout << "InvertMassMatrixDecorator:: NOT IMPLEMENTED FOR dfdp!! " << std::endl;
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
    //FIXME! this would not work in Thyra::ModelEvaluator!
    /*else {  // Mass Matrix is coeff of Second deriv
      modelInArgs.set_x_dotdot(x_dot);
      modelInArgs.set_alpha(0.0); 
      modelInArgs.set_beta(0.0);
      modelInArgs.set_omega(-1.0);
    }*/
    
    if (calcMassMatrix) {
      modelOutArgs.set_W_op(massMatrix);
    }

    //Evaluate the underlying model
    model->evalModel(modelInArgs, modelOutArgs);

    // Invert the mass matrix:   f_exp = M^{-1} f_imp

    if (!lumpMassMatrix) {
      // Create a linear solver based on the forward operator massMatrix
      if (calcMassMatrix) {
        A = massMatrix; //Teuchos::rcp(new Thyra::LinearOpWithSolveBase<double>(*massMatrix ));
        lows = Thyra::linearOpWithSolve(*lowsFactory, A);
      }

      // Solve the linear system for x, given b 
      ::Thyra::solve<double>(*lows, ::Thyra::NOTRANS, *modelOutArgs.get_f(), outArgs.get_f().ptr());
    }
    else { // Lump matrix into inverse of diagonal
    
      if (calcMassMatrix) {
        Thyra::put_scalar<Scalar>(1.0, invDiag.ptr());
        Thyra::apply<Scalar>(*massMatrix, Thyra::NOTRANS, *invDiag, invDiag.ptr(), 1.0, 0.0);
        Thyra::reciprocal<Scalar>(*invDiag, invDiag.ptr());
      }
      Teuchos::RCP<Thyra::VectorBase<Scalar> > f = outArgs.get_f(); 
      Thyra::Vp_StVtV<Scalar>(f.ptr(), 1.0, *invDiag, *modelOutArgs.get_f()); 
    }

    // Do not recompute mass matrix in future if it is a constant
    if (massMatrixIsConstant) calcMassMatrix = false;
  }
}
