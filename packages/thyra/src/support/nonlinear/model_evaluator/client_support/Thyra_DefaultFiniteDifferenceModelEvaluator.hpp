// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Array.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object and computes certain
 * derivatives using finite differences.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class DefaultFiniteDifferenceModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, fdStepLen )

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, fdStepLenIsRelative )

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultFiniteDifferenceModelEvaluator(
    const ScalarMag       fdStepLen            = -1.0
    ,const bool           fdStepLenIsRelative  = false
    );

  /** \brief . */
  DefaultFiniteDifferenceModelEvaluator(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel = NULL
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::DefaultFiniteDifferenceModelEvaluator(
  const ScalarMag       fdStepLen
  ,const bool           fdStepLenIsRelative
  )
  :fdStepLen_(fdStepLen)
  ,fdStepLenIsRelative_(fdStepLenIsRelative)
{}

template<class Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::DefaultFiniteDifferenceModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
  )
{
  initialize(thyraModel);
}

template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
}

template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::uninitialize(
  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
}

// Overridden from ModelEvaulator.

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np,Ng);
  outArgs.setSupports(wrappedOutArgs);
  // Just support derivatives of DgDp for now!
  for( int j = 0; j < Ng; ++j ) {
    for( int l = 0; l < Np; ++l ) {
      outArgs.setSupports(MEB::OUT_ARG_DgDp,j,l,MEB::DERIV_TRANS_MV_BY_ROW);
    }
  }
  // ToDo: Add support for more derivatives as needed!
  return outArgs;
}

template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  typedef typename ST::magnitudeType ScalarMag;

  typedef RefCountPtr<VectorBase<Scalar> >         V_ptr;
  typedef RefCountPtr<const VectorBase<Scalar> >   CV_ptr;
  typedef RefCountPtr<MultiVectorBase<Scalar> >    MV_ptr;

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar>::evalModel(...) ...\n";

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\ninArgs =\n" << Teuchos::describe(inArgs,verbLevel)
      << "\noutArgs on input =\n" << Teuchos::describe(outArgs,Teuchos::VERB_LOW);

  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  typedef Teuchos::VerboseObjectTempState<ModelEvaluatorBase> VOTSME;
  VOTSME thyraModel_outputTempState(thyraModel,out,verbLevel);

  //
  // Just do the g_0(p_0) case for now!
  //

  const RefCountPtr<const VectorSpaceBase<Scalar> >
    p_space = thyraModel->get_p_space(0),
    g_space = thyraModel->get_g_space(0);

  //
  // Compute the base point
  //

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the base point ...\n";

  MEB::InArgs<Scalar>  wrappedInArgs = inArgs;
  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  V_ptr g_base_ptr = ( outArgs.get_g(0).get() ? outArgs.get_g(0) : createMember(g_space) );
  wrappedOutArgs.set_g(0,g_base_ptr);
  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
  ConstDetachedVectorView<Scalar> g(*wrappedOutArgs.get_g(0));
  
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\ng="<<g[0]<<"\n";

  bool failed = false;
  if( ST::isnaninf(g[0]) )
    failed = true;
  
  //
  // Compute the variations
  //

  MV_ptr DgDp_trans = outArgs.get_DgDp(0,0).getDerivativeMultiVector().getMultiVector();
  if(DgDp_trans.get()) {
    if(!failed) {
      const V_ptr Dg0Dp_trans_ptr = DgDp_trans->col(0);
      DetachedVectorView<Scalar> Dg0Dp_trans(*Dg0Dp_trans_ptr);
      const ScalarMag eps = ( fdStepLen() > 0.0 ? fdStepLen() : ST::squareroot(ST::eps()) );
      V_ptr g_perturbed_ptr = createMember(g_space);
      wrappedOutArgs.set_g(0,g_perturbed_ptr);
      V_ptr p_perturbed_ptr = createMember(p_space);
      V_V(&*p_perturbed_ptr,*inArgs.get_p(0));
      wrappedInArgs.set_p(0,p_perturbed_ptr);
      ConstDetachedVectorView<Scalar> p(*inArgs.get_p(0));
      const int np = p_space->dim();
      for( int i = 0; i < np && !failed; ++i ) {
        const Scalar
          p_i = p[i],
          p_perturbed_i = p_i + eps;
        if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
          *out << "\nComputing perturbation i="<<i<<", p_i="<<p_i<<", eps="<<eps<<", p_perturbed_i=p_i+eta="<<p_perturbed_i<<" ...\n";
        set_ele(i,p_perturbed_i,&*p_perturbed_ptr);
        thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
        set_ele(i,p_i,&*p_perturbed_ptr);
        ConstDetachedVectorView<Scalar> g_perturbed(*wrappedOutArgs.get_g(0));
        const Scalar dgdp_i = (g_perturbed[0] - g[0])/eps;
        if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
          *out << "\ng_perturbed_i="<<g_perturbed[0]<<", dgdp_i=(g_perturbed_i-g)/eps="<<dgdp_i<<"\n";
        Dg0Dp_trans[i] = dgdp_i;
        if( ST::isnaninf(dgdp_i) )
          failed = true;
      }
    }
  }
  
/*

  //
  // Get the output objects
  //

  const int Ng = outArgs.Ng(), Np = outArgs.Np();
  typedef Teuchos::Array<bool>    bool_t;
  typedef Teuchos::Array<V_ptr>   V_ptr_array;
  typedef Teuchos::Array<MV_ptr>  MV_ptr_array;
  bool_t  perturb_p(Np,false);
  bool_t  perturb_g(Ng,false);
  V_ptr_array     g(Ng);
  MV_ptr_array_t  DgDp(Ng*Np);
  for( int j = 0; j < Ng; ++j ) {
    for( int l = 0; l < Np; ++l ) {
      if( outArgs.supports(MEB::OUT_ARG_DgDp,j,l) && outArgs.get_DgDp(j,l).get() ) {
        perturb_p[l] = true;
        perturb_g[j] = true;
      }
    }
    const V_ptr g_j = outArgs.get_g(j);
    g[j] = (
      g_j.get()
      ? g_j
      : ( perturb_g[j]
          ? createMember(thyraModel->get_g_space(j))
          : Teuchos::null
        )
      );
  }

  //
  // Compute the output functions at the base point 
  //
  
  MEB::InArgs<Scalar>  wrappedInArgs = inArgs;
  MEB::OutArgs<Scalar> wrappedOutArgs = outArgs;
  for( int j = 0; j < Ng; ++j ) {
    V_ptr = g_j;
    if((g_j=g[j].get())
       wrappedOutArgs.set_g(j,g_j));
  }
  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
  
  //
  // Do the perturbations of p!
  //
  
  V_ptr_array g_perturbed(Ng);
  for( int j = 0; j < Ng; ++j ) {
    if(perturb_g[j]) {
      g_perturbed[j] = createMember(thyraModel->get_g_space(j));
    }
  }
  for( int l = 0; l < Np; ++l ) {
    if(!perturb_p[l]) continue;
    RefCountPtr<const VectorBase<Scalar> >
      p_l = inArgs.get_p(l);
    RefCountPtr<const VectorBase<Scalar> >
      p_perturbed_l = createMember(p_l->space());
    assign(&*p_perturbed_l,*p_l);
    wrappedInArgs.set_p(l,p_perturbed_l);
    const int np_l = p_l->space()->dim();
    for( int l_i = 0; l_i < np_l; ++l_i ) {
      const Scalar p_l_i = get_ele(l_i,*p_l);
      const Scalar eps = fdStepLen(); // ToDo: Also consider relative step also
      const Scalar p_perturbed_l_i = p_l_i + eta;
      set_ele(l_i,p_perturbed_l_i,&*p_perturbed_l);
      for( int j = 0; j < Ng; ++j ) {
        if( outArgs.supports(MEB::OUT_ARG_DgDp,j,l) &&outArgs.get_DgDp(j,l).get() ) {
          wrappedOutArgs.set_g(j,g_perturbed[j]);
        }
      }
      thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
      for( int j = 0; j < Ng; ++j ) {
        wrappedOutArgs.set_g(j,Teuchos::null);
      }
      set_ele(l_i,p_l_i,&*p_perturbed_l);
      for( int j = 0; j < Ng; ++j ) {
        MV_ptr_t DgDp_j_l;
        if( outArgs.supports(MEB::OUT_ARG_DgDp,j,l) && (DgDp_j_j=outArgs.get_DgDp(j,l)).get() ) {
          V_VmV(&*DgDp_j_l->col(l_i),*g_perturbed[j],*g[j]);
        }
      }
    }
    wrappedInArgs.set_p(l,p_l);
  }
*/

  if(failed) {
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *out
        << "\nEvaluation failed, returning NaNs ...\n";
    outArgs.setFailed();
  }

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\noutArgs on output =\n" << Teuchos::describe(outArgs,verbLevel);
  
  totalTimer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal evaluation time = "<<totalTimer.totalElapsedTime()<<" sec\n"
      << "\nLeaving Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar>::evalModel(...) ...\n";
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultFiniteDifferenceModelEvaluator<Scalar>::description() const
{
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultFiniteDifferenceModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << ",fdStepLen="<<fdStepLen();
  oss << "}";
  return oss.str();
}

} // namespace Thyra

#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_HPP
