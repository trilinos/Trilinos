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

#ifdef HAVE_PIRO_IFPACK2
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

#ifdef HAVE_PIRO_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#endif


#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::InvertMassMatrixDecorator(
#else
template <typename Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::InvertMassMatrixDecorator(
#endif
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

#ifdef ALBANY_BUILD
#ifdef HAVE_PIRO_IFPACK2
  typedef Thyra::PreconditionerFactoryBase<double> Base;
  typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double, LocalOrdinal, GlobalOrdinal, Node> > Impl;
  linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef HAVE_PIRO_MUELU
  Stratimikos::enableMueLu<LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);
#endif
#else
#ifdef HAVE_PIRO_IFPACK2
  typedef Thyra::PreconditionerFactoryBase<double> Base;
  typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double> > Impl;
  linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef HAVE_PIRO_MUELU
  Stratimikos::enableMueLu(linearSolverBuilder);
#endif
#endif

   linearSolverBuilder.setParameterList(stratParams);

  // Create a linear solver factory given information read from the
  // parameter list.
  lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

  // Setup output stream and the verbosity level
  lowsFactory->setOStream(out);
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~InvertMassMatrixDecorator()
#else
template<typename Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::~InvertMassMatrixDecorator()
#endif
{
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_x_space() const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_x_space() const
#endif
{
  return model->get_x_space();
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_f_space() const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_f_space() const
#endif
{
  return model->get_f_space();
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_p_space(int l) const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_p_space(int l) const
#endif
{
  return model->get_p_space(l);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Teuchos::Array<std::string> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_p_names(int l) const
#else
template<typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
Piro::InvertMassMatrixDecorator<Scalar>::get_p_names(int l) const
#endif
{
  return model->get_p_names(l);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_g_space(int j) const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_g_space(int j) const
#endif
{
  return model->get_g_space(j);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::ArrayView<const std::string>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_g_names(int j) const
#else
template<typename Scalar>
Teuchos::ArrayView<const std::string>
Piro::InvertMassMatrixDecorator<Scalar>::get_g_names(int j) const
#endif
{
  return model->get_g_names(j);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNominalValues() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::getNominalValues() const
#endif
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues = this->createInArgsImpl();
  nominalValues.setArgs(
      model->getNominalValues(),
      /* ignoreUnsupported = */ true,
      /* cloneObjects = */ false);
  return nominalValues;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP< Thyra::LinearOpBase< Scalar > >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::create_W_op () const
#else
template<typename Scalar>
Teuchos::RCP< Thyra::LinearOpBase< Scalar > >
Piro::InvertMassMatrixDecorator<Scalar>::create_W_op () const
#endif
{
  return model->create_W_op();
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_W_factory() const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::get_W_factory() const
#endif
{
  return model->get_W_factory();
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::create_W_prec() const
#else
template<typename Scalar>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar> >
Piro::InvertMassMatrixDecorator<Scalar>::create_W_prec() const
#endif
{
  return model->create_W_prec();
}


#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLowerBounds() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::getLowerBounds() const
#endif
{
  return Thyra::ModelEvaluatorBase::InArgs<Scalar>(); // Default value
}


#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getUpperBounds() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::getUpperBounds() const
#endif
{
  return Thyra::ModelEvaluatorBase::InArgs<Scalar>(); // Default value
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::reportFinalPoint(
#else
template<typename Scalar>
void
Piro::InvertMassMatrixDecorator<Scalar>::reportFinalPoint(
#endif
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint,
    const bool wasSolved)
{
  model->reportFinalPoint(finalPoint,wasSolved);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createInArgs() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::createInArgs() const
#endif
{
  return this->createInArgsImpl();
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createInArgsImpl() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::createInArgsImpl() const
#endif
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> result = model->createInArgs();
  result.setModelEvalDescription(this->description());
  return result; 
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createOutArgsImpl() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::InvertMassMatrixDecorator<Scalar>::createOutArgsImpl() const
#endif
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> result = model->createOutArgs();
  result.setModelEvalDescription(this->description());
  return result; 
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::evalModelImpl(
#else
template<typename Scalar>
void
Piro::InvertMassMatrixDecorator<Scalar>::evalModelImpl(
#endif
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
        A = massMatrix; //Teuchos::rcp(new Thyra::LinearOpWithSolveBase<double>(*massMatrix ));
        lows = Thyra::linearOpWithSolve(*lowsFactory, A);
      }
      else { // Lump matrix into inverse of diagonal
        Thyra::put_scalar<Scalar>(1.0, invDiag.ptr());
        Thyra::apply<Scalar>(*massMatrix, Thyra::NOTRANS, *invDiag, invDiag.ptr(), 1.0, 0.0);
        Thyra::reciprocal<Scalar>(*invDiag, invDiag.ptr());
        //IKT, 5/31/17: adding the following logic which checks invDiag vector for nans
        //and throws an exception if nans are found.
        typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType ScalarMag;
        typedef Teuchos::ScalarTraits<ScalarMag> SMT;
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

    //set f and unset W_op in modelOutArgs
    //This is to avoid calling getting Jacobian and residual at the same time, 
    //which we need to do for Aeras problems calling this function from Albany. 
    modelOutArgs.set_f(outArgs.get_f()); 
    modelOutArgs.set_W_op(Teuchos::null); 
    
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
