// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MODEL_EVALUATOR_HELPERS_HPP
#define THYRA_MODEL_EVALUATOR_HELPERS_HPP


#include "Thyra_ModelEvaluator.hpp"


namespace Thyra {


/** \brief Create a clone of an InArgs object.
 *
 * Warning!  This function only creates a shallow copy of the underlying input
 * objects.  Therefore, be careful if you try to modify any of these.
 *
 * \relates ModelEvaluatorDefaultBase
 */
template<class Scalar>
RCP<ModelEvaluatorBase::InArgs<Scalar> >
clone( const ModelEvaluatorBase::InArgs<Scalar> &inArgs );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
derivativeGradient(
  const RCP<MultiVectorBase<Scalar> > &grad
  );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DfDp_mv(
  const ModelEvaluator<Scalar>& model,
  int l,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DgDx_dot_mv(
  const ModelEvaluator<Scalar>& model,
  int j,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DgDx_mv(
  const ModelEvaluator<Scalar>& model,
  int j,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DgDp_mv(
  const ModelEvaluator<Scalar>& model,
  int j,
  int l,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
get_dmv(
  const ModelEvaluatorBase::Derivative<Scalar> &deriv
  ,const std::string &derivName
  );


/** \relates ModelEvaluatorDefaultBase */
template<class Scalar>
RCP<MultiVectorBase<Scalar> >
get_mv(
  const ModelEvaluatorBase::Derivative<Scalar> &deriv
  ,const std::string &derivName
  ,ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


/** \brief Assert that that Thyra objects imbedded in a Derivative object
 * matches its function and variable spaces.
 *
 * \relates ModelEvaluatorDefaultBase
 */
template<class Scalar>
void assertDerivSpaces(
  const std::string &modelEvalDescription,
  const ModelEvaluatorBase::Derivative<Scalar> &deriv,
  const std::string &deriv_name,
  const VectorSpaceBase<Scalar> &fnc_space,
  const std::string &fnc_space_name,
  const VectorSpaceBase<Scalar> &var_space,
  const std::string &var_space_name
  );


/** \brief Assert that an InArgs and OutArgs object are setup consistently.
 *
 * \relates ModelEvaluatorDefaultBase
 */
template<class Scalar>
void assertInArgsOutArgsSetup(
  const std::string &modelEvalDescription,
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  );


/** \brief Assert that the objects in an InArgs object match a given model.
 *
 * \relates ModelEvaluatorDefaultBase
 */
template<class Scalar>
void assertInArgsEvalObjects(
  const ModelEvaluator<Scalar> &model,
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs
  );


/** \brief Assert that the objects in an OutArgs object match a given model.
 *
 * \relates ModelEvaluatorDefaultBase
 */
template<class Scalar>
void assertOutArgsEvalObjects(
  const ModelEvaluator<Scalar> &model,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs,
  const ModelEvaluatorBase::InArgs<Scalar> *inArgs = 0
  );


/** \brief Evaluate <tt>f(x)</tt>. */
template<class Scalar>
void eval_f(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x
  ,VectorBase<Scalar> *f
  );


/** \brief Evaluate <tt>f(x)</tt> and <tt>W(x) = DfDx(x)</tt>. */
template<class Scalar>
void eval_f_W(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x
  ,VectorBase<Scalar> *f
  ,LinearOpWithSolveBase<Scalar> *W
  );


/** \brief Evaluate <tt>f(x,t)</tt>. */
template<class Scalar>
void eval_f(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x
  ,const Scalar &t
  ,VectorBase<Scalar> *f
  );


/** \brief Evaluate <tt>g(j)(p))</tt>. */
template<class Scalar>
void eval_g(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const int j,
  const Ptr<VectorBase<Scalar> > &g_j
  );

#ifndef THYRA_HIDE_DEPRECATED_CODE
/** \brief Deprecated . */
template<class Scalar>
THYRA_DEPRECATED
void eval_g(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const int j,
  VectorBase<Scalar> *g_j
  )
{
  eval_g(model, l, p_l, j, Teuchos::ptr(g_j));
}
#endif // THYRA_HIDE_DEPRECATED_CODE

/** \brief Evaluate <tt>g(j)(p,t))</tt>. */
template<class Scalar>
void eval_g(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const Scalar &t,
  const int j,
  VectorBase<Scalar> *g_j
  );


/** \brief Evaluate <tt>g(j)(p))</tt> and/or D(g)/D(p). */
template<class Scalar>
void eval_g_DgDp(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const int j,
  const Ptr<VectorBase<Scalar> > &g_j,
  const ModelEvaluatorBase::Derivative<Scalar> &DgDp_j_l
  );


/** \brief Evaluate <tt>f(x_dot,x,t)</tt>. */
template<class Scalar>
void eval_f(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x_dot
  ,const VectorBase<Scalar> &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,VectorBase<Scalar> *f
  );


/** \brief Evaluate <tt>f(x_dot,x,t)</tt> and <tt>W(x_dot,x,t,alpha,beta) =
 * alpha*DfDx_dot(x_dot,x,t) + beta*DfDx(x_dot,x,t)</tt>. */
template<class Scalar>
void eval_f_W(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x_dot
  ,const VectorBase<Scalar> &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,const Scalar &alpha
  ,const Scalar &beta
  ,VectorBase<Scalar> *f
  ,LinearOpWithSolveBase<Scalar> *W
  );


#ifdef HAVE_THYRA_ME_POLYNOMIAL


/** \brief . */
template<class Scalar>
void eval_f_poly(
  const ModelEvaluator<Scalar> &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> > &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,Teuchos::Polynomial< VectorBase<Scalar> > *f_poly
  );


/** \brief . */
template<class Scalar>
void eval_f_poly(
  const ModelEvaluator<Scalar> &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> > &x_dot_poly
  ,const VectorBase<Scalar> &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,Teuchos::Polynomial< VectorBase<Scalar> > *f_poly
  );


#endif // HAVE_THYRA_ME_POLYNOMIAL


} // namespace Thyra


//
// Implementations
//


#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"


template<class Scalar>
Teuchos::RCP<Thyra::ModelEvaluatorBase::InArgs<Scalar> >
Thyra::clone( const ModelEvaluatorBase::InArgs<Scalar> &inArgs )
{
  RCP<ModelEvaluatorBase::InArgs<Scalar> >
    newInArgs = Teuchos::rcp(new ModelEvaluatorBase::InArgs<Scalar>);
  *newInArgs = inArgs;
  return newInArgs;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::Derivative<Scalar>
Thyra::derivativeGradient(
  const RCP<MultiVectorBase<Scalar> > &grad
  )
{
  return ModelEvaluatorBase::Derivative<Scalar>(
    grad,
    ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM
    );
}


template<class Scalar>
Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
Thyra::create_DfDp_mv(
  const ModelEvaluator<Scalar>& model,
  int l,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(!(orientation==ModelEvaluatorBase::DERIV_MV_BY_COL));
  return createMembers( model.get_f_space(), model.get_p_space(l)->dim() );
}


template<class Scalar>
Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
Thyra::create_DgDx_dot_mv(
  const ModelEvaluator<Scalar>& model,
  int j,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  typedef ModelEvaluatorBase MEB;
  switch(orientation) {
    case MEB::DERIV_MV_BY_COL:
      return
        MEB::DerivativeMultiVector<Scalar>(
          createMembers( model.get_g_space(j), model.get_x_space()->dim() )
          ,MEB::DERIV_MV_BY_COL
          );
    case MEB::DERIV_TRANS_MV_BY_ROW:
      return
        MEB::DerivativeMultiVector<Scalar>(
          createMembers( model.get_x_space(), model.get_g_space(j)->dim() )
          ,MEB::DERIV_TRANS_MV_BY_ROW
          );
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return MEB::DerivativeMultiVector<Scalar>(); // Never executed!
}


template<class Scalar>
Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
Thyra::create_DgDx_mv(
  const ModelEvaluator<Scalar>& model,
  int j,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  return create_DgDx_dot_mv(model,j,orientation);
}


template<class Scalar>
Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
Thyra::create_DgDp_mv(
  const ModelEvaluator<Scalar>& model,
  int j,
  int l,
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  typedef ModelEvaluatorBase MEB;
  switch(orientation) {
    case MEB::DERIV_MV_BY_COL:
      return
        MEB::DerivativeMultiVector<Scalar>(
          createMembers( model.get_g_space(j), model.get_p_space(l)->dim() )
          ,MEB::DERIV_MV_BY_COL
          );
    case MEB::DERIV_TRANS_MV_BY_ROW:
      return
        MEB::DerivativeMultiVector<Scalar>(
          createMembers( model.get_p_space(l), model.get_g_space(j)->dim() )
          ,MEB::DERIV_TRANS_MV_BY_ROW
          );
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return MEB::DerivativeMultiVector<Scalar>(); // Never executed!
}


template<class Scalar>
Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
Thyra::get_dmv(
  const ModelEvaluatorBase::Derivative<Scalar> &deriv
  ,const std::string &derivName
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get()!=NULL, std::logic_error
    ,"Error, LinearOpBase type not expected for " << derivName <<"!"
    );
  return deriv.getDerivativeMultiVector();
}


template<class Scalar>
Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
Thyra::get_mv(
  const ModelEvaluatorBase::Derivative<Scalar> &deriv
  ,const std::string &derivName
  ,ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  typedef ModelEvaluatorBase MEB;
  TEUCHOS_TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get()!=NULL, std::logic_error
    ,"Error, LinearOpBase type not expected for " << derivName <<"!"
    );
  MEB::DerivativeMultiVector<Scalar>
    dmv = deriv.getDerivativeMultiVector();
  RCP<MultiVectorBase<Scalar> >
    mv = dmv.getMultiVector();
  if( mv.get() ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      dmv.getOrientation() != orientation, std::logic_error
      ,"Error, the orientation " << toString(dmv.getOrientation()) << " is not the"
      " expected orientation of " << toString(orientation)
      << " for " << derivName << "!"
      );
  }
  return mv;
}


template<class Scalar>
void Thyra::assertDerivSpaces(
  const std::string &modelEvalDescription,
  const ModelEvaluatorBase::Derivative<Scalar> &deriv,
  const std::string &deriv_name,
  const VectorSpaceBase<Scalar> &fnc_space,
  const std::string &fnc_space_name,
  const VectorSpaceBase<Scalar> &var_space,
  const std::string &var_space_name
  )
{
  typedef ModelEvaluatorBase MEB;
  if (!is_null(deriv.getLinearOp())) {
    const RCP<const LinearOpBase<Scalar> > lo = deriv.getLinearOp();
    if (!is_null(lo->range())) {
      THYRA_ASSERT_VEC_SPACES_NAMES(
        modelEvalDescription,
        *lo->range(), deriv_name + ".range()",
        fnc_space, fnc_space_name
        );
      THYRA_ASSERT_VEC_SPACES_NAMES(
        modelEvalDescription,
        *lo->domain(), deriv_name + ".domain()",
        var_space, var_space_name
        );
    }
  }
  else if(!is_null(deriv.getMultiVector())) {
    const RCP<const LinearOpBase<Scalar> > mv = deriv.getMultiVector();
    switch(deriv.getMultiVectorOrientation()) {
      case MEB::DERIV_MV_BY_COL: {
        THYRA_ASSERT_VEC_SPACES_NAMES(
          modelEvalDescription,
          *mv->range(), deriv_name + ".range()",
          fnc_space, fnc_space_name
          );
        THYRA_ASSERT_VEC_SPACES_NAMES(
          modelEvalDescription,
          *mv->domain(), deriv_name + ".domain()",
          var_space, var_space_name
          );
        break;
      }
      case MEB::DERIV_TRANS_MV_BY_ROW: {
        THYRA_ASSERT_VEC_SPACES_NAMES(
          modelEvalDescription,
          *mv->range(), deriv_name + "^T.range()",
          var_space, var_space_name
          );
        THYRA_ASSERT_VEC_SPACES_NAMES(
          modelEvalDescription,
          *mv->domain(), deriv_name + "^T.domain()",
          fnc_space, fnc_space_name
          );
        break;
      }
#ifdef TEUCHOS_DEBUG
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
    }
  }
}


template<class Scalar>
void Thyra::assertInArgsOutArgsSetup(
  const std::string &modelEvalDescription,
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  )
{

  typedef ModelEvaluatorBase MEB;

  const int Ng = outArgs.Ng();
  const int Np = outArgs.Np();

  // Description
  TEUCHOS_ASSERT_EQUALITY(inArgs.modelEvalDescription(), modelEvalDescription);
  TEUCHOS_ASSERT_EQUALITY(outArgs.modelEvalDescription(), modelEvalDescription);

  // Np
  TEUCHOS_TEST_FOR_EXCEPTION(
    inArgs.Np() != outArgs.Np(), std::logic_error,
    "Error: The underlying model " << modelEvalDescription << " incorrectly\n"
    "set inArgs.Np() = "<<inArgs.Np()<<" != outArgs.Np() = "
    <<outArgs.Np()<<"!"
    );

  // x_dot
  TEUCHOS_TEST_FOR_EXCEPTION(
    inArgs.supports(MEB::IN_ARG_x_dot) && !inArgs.supports(MEB::IN_ARG_x),
    std::logic_error,
    "Error: The underlying model " << modelEvalDescription << " supports\n"
    "x_dot but does not support x!"
    );

  // t
  TEUCHOS_TEST_FOR_EXCEPTION(
    inArgs.supports(MEB::IN_ARG_x_dot) && !inArgs.supports(MEB::IN_ARG_t),
    std::logic_error,
    "Error: The underlying model " << modelEvalDescription << " supports\n"
    "x_dot but does not support t!"
    );

  // W and W_op
  TEUCHOS_TEST_FOR_EXCEPTION(
    (
      ( outArgs.supports(MEB::OUT_ARG_W) || outArgs.supports(MEB::OUT_ARG_W_op) )
      &&
      !inArgs.supports(MEB::IN_ARG_x)
      ),
    std::logic_error,
    "Error: The underlying model " << modelEvalDescription << " says that\n"
    "it supports W and/or W_op but it does not support x!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    (
      ( outArgs.supports(MEB::OUT_ARG_W) || outArgs.supports(MEB::OUT_ARG_W_op) )
      &&
      inArgs.supports(MEB::IN_ARG_x_dot)
      &&
      !( inArgs.supports(MEB::IN_ARG_alpha) && inArgs.supports(MEB::IN_ARG_beta) )
      ),
    std::logic_error,
    "Error: The underlying model " << modelEvalDescription << " supports W and/or W_op\n"
    "and x_dot but it does not support alpha and beta as InArgs!"
    );

  for ( int l = 0; l < Np; ++l ) {

    // DfDp(l): OutArgs checks this automatically!

    for ( int j = 0; j < Ng; ++j ) {

      // DgDx_dot(j)
      TEUCHOS_TEST_FOR_EXCEPTION(
        ( !outArgs.supports(MEB::OUT_ARG_DgDx_dot,j).none()
          && !inArgs.supports(MEB::IN_ARG_x_dot) ),
        std::logic_error,
        "Error: The underlying model " << modelEvalDescription << " says that\n"
        "it supports DgDx_dot("<<j<<") but it does not support x_dot!"
        );

      // DgDx(j)
      TEUCHOS_TEST_FOR_EXCEPTION(
        ( !outArgs.supports(MEB::OUT_ARG_DgDx,j).none()
          && !inArgs.supports(MEB::IN_ARG_x) ),
        std::logic_error,
        "Error: The underlying model " << modelEvalDescription << " says that\n"
        "it supports DgDx("<<j<<") but it does not support x!"
        );

      // DgDp(j,l): OutArgs checks this automatically!

    }

  }

}


template<class Scalar>
void Thyra::assertInArgsEvalObjects(
  const ModelEvaluator<Scalar> &model,
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs
  )
{
  
  typedef ModelEvaluatorBase MEB;

  const std::string description = model.description();
  const int Np = inArgs.Np();

  model.createInArgs().assertSameSupport(inArgs);

  // x_dot
  if ( inArgs.supports(MEB::IN_ARG_x_dot) && !is_null(inArgs.get_x_dot()) ) {
    THYRA_ASSERT_VEC_SPACES(
      description, *inArgs.get_x_dot()->space(), *model.get_x_space() );
  }

  // x
  if ( inArgs.supports(MEB::IN_ARG_x) && !is_null(inArgs.get_x()) ) {
    THYRA_ASSERT_VEC_SPACES(
      description, *inArgs.get_x()->space(), *model.get_x_space() );
  }
    
  // p(l)
  for ( int l = 0; l < Np; ++l ) {
    if (!is_null(inArgs.get_p(l))) {
      THYRA_ASSERT_VEC_SPACES(
        description, *inArgs.get_p(l)->space(), *model.get_p_space(l) );
    }
  }

}


template<class Scalar>
void Thyra::assertOutArgsEvalObjects(
  const ModelEvaluator<Scalar> &model,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs,
  const ModelEvaluatorBase::InArgs<Scalar> *inArgs
  )
{

  typedef ScalarTraits<Scalar> ST;
  typedef Teuchos::Utils TU;
  typedef ModelEvaluatorBase MEB;

  const std::string description = model.description();
  const int Ng = outArgs.Ng();
  const int Np = outArgs.Np();

  if (inArgs) {
    TEUCHOS_ASSERT_EQUALITY(outArgs.Np(), inArgs->Np());
  }

  model.createOutArgs().assertSameSupport(outArgs);

  // f
  if ( outArgs.supports(MEB::OUT_ARG_f) && !is_null(outArgs.get_f()) ) {
    THYRA_ASSERT_VEC_SPACES(
      description, *outArgs.get_f()->space(), *model.get_f_space() );
  }

  // W
  if ( outArgs.supports(MEB::OUT_ARG_W) && !is_null(outArgs.get_W()) ) {
    if (!is_null(outArgs.get_W()->range())) {
      THYRA_ASSERT_VEC_SPACES(
        description, *outArgs.get_W()->range(), *model.get_f_space() );
      THYRA_ASSERT_VEC_SPACES(
        description, *outArgs.get_W()->domain(), *model.get_x_space() );
    }
  }
    
  // W_op
  if ( outArgs.supports(MEB::OUT_ARG_W_op) && !is_null(outArgs.get_W_op()) ) {
    if (!is_null(outArgs.get_W_op()->range())) {
      THYRA_ASSERT_VEC_SPACES(
        description, *outArgs.get_W_op()->range(), *model.get_f_space() );
      THYRA_ASSERT_VEC_SPACES(
        description, *outArgs.get_W_op()->domain(), *model.get_x_space() );
    }
  }

  // alpha and beta (not really in outArgs but can only be validated if W or
  // W_op is set)
  if (
    inArgs
    &&
    (
      ( outArgs.supports(MEB::OUT_ARG_W) && !is_null(outArgs.get_W()) )
      ||
      ( outArgs.supports(MEB::OUT_ARG_W_op) && !is_null(outArgs.get_W_op()) )
      )
    )
  {
    if ( inArgs->supports(MEB::IN_ARG_alpha) && inArgs->supports(MEB::IN_ARG_beta) ) {
      // 08/25/08 tscoffe:  In the block-composed linear operator case for
      // Rythmos::ImplicitRKModelEvaluator, I need to specify that a given
      // block is all zeros and I'm depending on the underlying model to
      // intelligently fill the block with zeros if both alpha and beta are
      // zero.  
      //TEUCHOS_TEST_FOR_EXCEPT( inArgs->get_alpha() == ST::zero() && inArgs->get_beta() == ST::zero() );
    }
    else if ( inArgs->supports(MEB::IN_ARG_beta) ) {
      TEUCHOS_TEST_FOR_EXCEPT( inArgs->get_beta() == ST::zero() );
    }
  }

  // DfDp(l)
  if (outArgs.supports(MEB::OUT_ARG_f)) {
    for ( int l = 0; l < Np; ++l ) {
      if (!outArgs.supports(MEB::OUT_ARG_DfDp,l).none()) {
        assertDerivSpaces(
          description,
          outArgs.get_DfDp(l), "DfDp("+TU::toString(l)+")",
          *model.get_f_space(), "f_space",
          *model.get_p_space(l), "p_space("+TU::toString(l)+")"
          );
      }
    }
  }
    
  // g(l)
  for ( int j = 0; j < Ng; ++j ) {
    if (!is_null(outArgs.get_g(j))) {
      THYRA_ASSERT_VEC_SPACES(
        description, *outArgs.get_g(j)->space(), *model.get_g_space(j) );
    }
  }

  // DgDx_dot(j)
  for ( int j = 0; j < Ng; ++j ) {
    if (!outArgs.supports(MEB::OUT_ARG_DgDx_dot,j).none()) {
      assertDerivSpaces(
        description,
        outArgs.get_DgDx_dot(j), "DgDx_dot("+TU::toString(j)+")",
        *model.get_g_space(j), "g_space("+TU::toString(j)+")",
        *model.get_x_space(), "x_space"
        );
    }
  }

  // DgDx(j)
  for ( int j = 0; j < Ng; ++j ) {
    if (!outArgs.supports(MEB::OUT_ARG_DgDx,j).none()) {
      assertDerivSpaces(
        description,
        outArgs.get_DgDx(j), "DgDx("+TU::toString(j)+")",
        *model.get_g_space(j), "g_space("+TU::toString(j)+")",
        *model.get_x_space(), "x_space"
        );
    }
  }

  // Assert DgDp(j,l)
  for ( int j = 0; j < Ng; ++j ) {
    for ( int l = 0; l < Np; ++l ) {
      if (!outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none()) {
        const std::string j_str = TU::toString(j);
        const std::string l_str = TU::toString(l);
        assertDerivSpaces(
          description,
          outArgs.get_DgDp(j,l), "DgDp("+j_str+","+l_str+")",
          *model.get_g_space(j), "g_space("+j_str+")",
          *model.get_p_space(l), "p_space("+l_str+")"
          );
      }
    }
  }

}


template<class Scalar>
void Thyra::eval_f(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x
  ,VectorBase<Scalar> *f
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();
  inArgs.set_x(Teuchos::rcp(&x,false));
  outArgs.set_f(Teuchos::rcp(f,false));
  model.evalModel(inArgs,outArgs);
}


template<class Scalar>
void Thyra::eval_f_W(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x
  ,VectorBase<Scalar> *f
  ,LinearOpWithSolveBase<Scalar> *W
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();

  inArgs.set_x(Teuchos::rcp(&x,false));

  if (f) outArgs.set_f(Teuchos::rcp(f,false));
  if (W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}


template<class Scalar>
void Thyra::eval_f(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x
  ,const Scalar &t
  ,VectorBase<Scalar> *f
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(t);
  outArgs.set_f(Teuchos::rcp(f,false));
  model.evalModel(inArgs,outArgs);
}


template<class Scalar>
void Thyra::eval_g(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const int j,
  const Ptr<VectorBase<Scalar> > &g_j
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs= model.createOutArgs();
  inArgs.set_p(l, Teuchos::rcpFromRef(p_l));
  outArgs.set_g(j, Teuchos::rcpFromRef(*g_j));
  model.evalModel(inArgs,outArgs);
}


template<class Scalar>
void Thyra::eval_g(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const Scalar &t,
  const int j,
  VectorBase<Scalar> *g_j
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs= model.createOutArgs();
  inArgs.set_p(l,Teuchos::rcp(&p_l,false));
  inArgs.set_t(t);
  outArgs.set_g(j,Teuchos::rcp(g_j,false));
  model.evalModel(inArgs,outArgs);
}


template<class Scalar>
void Thyra::eval_g_DgDp(
  const ModelEvaluator<Scalar> &model,
  const int l,
  const VectorBase<Scalar> &p_l,
  const int j,
  const Ptr<VectorBase<Scalar> > &g_j,
  const ModelEvaluatorBase::Derivative<Scalar> &DgDp_j_l
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs= model.createOutArgs();
  inArgs.set_p(l, Teuchos::rcpFromRef(p_l));
  if (!is_null(g_j)) {
    outArgs.set_g(j, Teuchos::rcpFromPtr(g_j));
  }
  if (!DgDp_j_l.isEmpty()) {
    outArgs.set_DgDp(j, l, DgDp_j_l);
  }
  model.evalModel(inArgs,outArgs);
}


template<class Scalar>
void Thyra::eval_f(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x_dot
  ,const VectorBase<Scalar> &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,VectorBase<Scalar> *f
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();

  inArgs.set_x_dot(Teuchos::rcp(&x_dot,false));
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f(Teuchos::rcp(f,false));

  model.evalModel(inArgs,outArgs);

}


template<class Scalar>
void Thyra::eval_f_W(
  const ModelEvaluator<Scalar> &model
  ,const VectorBase<Scalar> &x_dot
  ,const VectorBase<Scalar> &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,const Scalar &alpha
  ,const Scalar &beta
  ,VectorBase<Scalar> *f
  ,LinearOpWithSolveBase<Scalar> *W
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();

  inArgs.set_x_dot(Teuchos::rcp(&x_dot,false));
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);

  if(f) outArgs.set_f(Teuchos::rcp(f,false));
  if(W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}


#ifdef HAVE_THYRA_ME_POLYNOMIAL


template<class Scalar>
void Thyra::eval_f_poly(
  const ModelEvaluator<Scalar> &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> > &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,Teuchos::Polynomial< VectorBase<Scalar> > *f_poly
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();

  inArgs.set_x_poly(Teuchos::rcp(&x_poly,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f_poly(Teuchos::rcp(f_poly,false));

  model.evalModel(inArgs,outArgs);

}


template<class Scalar>
void Thyra::eval_f_poly(
  const ModelEvaluator<Scalar> &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> > &x_dot_poly
  ,const VectorBase<Scalar> &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t
  ,Teuchos::Polynomial< VectorBase<Scalar> > *f_poly
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inArgs = model.createInArgs();
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();

  inArgs.set_x_dot_poly(Teuchos::rcp(&x_dot_poly,false));
  inArgs.set_x_poly(Teuchos::rcp(&x_poly,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f_poly(Teuchos::rcp(f_poly,false));

  model.evalModel(inArgs,outArgs);

}


#endif // HAVE_THYRA_ME_POLYNOMIAL


#endif // THYRA_MODEL_EVALUATOR_HELPERS_HPP
