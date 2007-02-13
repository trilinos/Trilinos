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

#ifndef THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP
#define THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_ParameterDrivenMultiVectorInput.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object and adds a simple, but
 * fairly general, inverse-type of response function.
 *
 * The following response function is added to the end of the supported
 * response functions:

 \verbatim

  g_(getUnderlyingModel()->Ng())(x,p,...)
    = observationMultiplier*(0.5/no)*sum((os(ko)*(o(ko)(x,p)-ot(lo)))^2,lo=0...no-1)
    + parameterMultiplier*(0.5/np)*sum((ps(kp)*(p(kp)-pt(kp)))^2,i=kp...np-1)

 \endverbatim 

 * where <tt>ot</tt> is the target vector for some observation (see below),
 * <tt>o</tt> is a scaling vector (where some components may be zero),
 * <tt>observationMultiplier</tt> is a scalar for the observation matching
 * term, <tt>p</tt> is one of the parameter subvectors supported by the
 * underlying model, <tt>pt</tt> is a nomial parameter vector for which
 * violations are penalized, <tt>ps</tt> is the scaling vector for the
 * parameter change violations, and <tt>parameterMultiplier</tt> is a scalar for
 * the parameter change violations.
 *
 * The observation function <tt>o(x,p)</tt> can be the state vector itself
 * <tt>o(x,p) = x</tt> for <tt>obs_idx < 0</tt>, or can be any of the built-in
 * response functions <tt>o(x,p) = g(obs_idx)(x,p)</tt> where <tt>0 <= obs_idx
 * < getUnderlyingModel()->Ng()</tt>.
 *
 * The first set of terms multiplying <tt>observationMultiplier</tt> is called
 * the <em>observation matching term</tt>.  The second set of terms
 * multiplying <tt>parameterMultiplier</tt> is called the <em>parameter
 * regularization term</tt>.
 *
 * Note that <tt>this->Ng() == getUnderlyingModel()->Ng() + 1</tt>.
 *
 * Let's consider the derivatives of this inverse function, which will just
 * refer to here as the scalar function <tt>g(x,p)</tt> (leaving out scaling
 * vectors <tt>os</tt> and <tt>ps</tt> for now) which takes the form:

 \verbatim

  g(x,p)
    = observationMultiplier*(0.5/no)*sum(diff_o(ko)(x,p))^2,i=ko...no-1)
      + parameterMultiplier*(0.5/np)*sum(diff_p(kp)(p))^2,kp=0...np-1)

 \endverbatim 

 * where <tt>diff_o(ko)(x,p) = (o(ko)(x,p) - ot(ko))</tt> and
 * <tt>diff_p(kp)(p) = (p(kp) - pt(kp))</tt>.
 *
 * The derivatives of <tt>g(x,p)</tt> with respect to <tt>x(kx)</tt> and
 * <tt>p(kp)</tt> are the scalar values:

 \verbatim

   DgDx(ix) = observationMultiplier*(1/no)*sum( (diff_o(ko)*D(o(ko))/D(x(ix))), ko=0...no-1 )

            = (observationMultiplier*(1/no)) * diff_o^T * D(o)/D(x(ix))

   DgDp(ip) = observationMultiplier*(1/no)*sum( (diff_o(ko)*D(o(ko))/D(p(ip))), ko=0...no-1 )
              + parameterMultiplier*(1.0/np) * diff_p(ip)

            = (observationMultiplier*(1/no)) * diff_o^T * D(o)/D(p(ip))
              + (parameterMultiplier*(1.0/np)) * diff_p(ip)

 \endverbatim 

 * The full derivatives are then:

 \verbatim

   DgDx = (observationMultiplier*(1/no)) * diff_o^T * DoDx

   DgDp = (observationMultiplier*(1/no)) * diff_o^T * DoDp
          + (parameterMultiplier*(1.0/np)) * diff_p^T

 \endverbatim 

 * In the gradient form:

 \verbatim

   DgDx^T = (observationMultiplier*(1/no)) * (DoDx^T) * diff_o

   DgDp^T = (observationMultiplier*(1/no)) * (DoDp^T) * diff_o
            + (parameterMultiplier*(1.0/np)) * diff_p

 \endverbatim

 * Note that when <tt>obs_idx < -1</tt> than <tt>o(x,p) = x</tt> and therefore
 * <tt>DoDx = I</tt> and <tt>DoDp = 0</tt>. and there derivatives simplify to:

 \verbatim

   DgDx^T = (observationMultiplier*(1/no)) * diff_o

   DgDp^T = (parameterMultiplier*(1.0/np)) * diff_p

 \endverbatim

 * When <tt>obs_idx >= 0</tt>, this implementation currently requires that
 * <tt>(DoDx^T)</tt> and <tt>(DoDp^T)</tt> be computed and returned by the
 * underlying model as multi-vector objects.  In the future, we really only
 * need the action of <tt>DoDx^T</tt> and <tt>DoDp^T</tt> onto vectors as
 * shown above.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class DefaultInverseModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
  , virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief Observation target vector <tt>ot</tt>. */
  STANDARD_COMPOSITION_MEMBERS( VectorBase<Scalar>, observationTarget );

  /** \brief Parameter base vector <tt>pt</tt>. */
  STANDARD_COMPOSITION_MEMBERS( VectorBase<Scalar>, parameterBase );

  /** \brief MultiVectorFileIOBase object used to read the observation target
   * vector <tt>ot</tt> as directed by the parameter list. */
  STANDARD_NONCONST_COMPOSITION_MEMBERS( MultiVectorFileIOBase<Scalar>, observationTargetIO );

  /** \brief MultiVectorFileIOBase object used to read the parameter base
   * vector <tt>pt</tt> as directed by the parameter list. */
  STANDARD_NONCONST_COMPOSITION_MEMBERS( MultiVectorFileIOBase<Scalar>, parameterBaseIO );

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultInverseModelEvaluator();

  /** \brief . */
  DefaultInverseModelEvaluator(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >       &thyraModel
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >       &thyraModel
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> >             *thyraModel
    );

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief .
   *
   * Note that <tt>observationTargetIO()</tt> and <tt>parameterBaseIO()</tt>
   * must be set before calling this function in order to use the parameter
   * sublist to read in the vectors <tt>observationTarget()</tt> and
   * <tt>parameterBase()</tt>.
   */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief .
   *
   * Note that <tt>observationTargetIO()</tt> and <tt>parameterBaseIO()</tt>
   * must be set before calling this function in order to have the sublists
   * added that will allow the vectors <tt>observationTarget()</tt> and
   * <tt>parameterBase()</tt> to be read in latter when the parameter list is
   * set..
   */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
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

private:

  // ////////////////////////////////
  // Private data members

  mutable Teuchos::RefCountPtr<const Teuchos::ParameterList> validParamList_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>  paramList_;

  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > inv_g_space_;

  int            obs_idx_;
  int            p_idx_;
  int            Ng_;

  double         observationMultiplier_;
  double         parameterMultiplier_; 

  mutable ParameterDrivenMultiVectorInput<Scalar> observationTargetReader_;
  mutable ParameterDrivenMultiVectorInput<Scalar> parameterBaseReader_;

  static const std::string ObservationIndex_name_;
  static const int ObservationIndex_default_;

  static const std::string ParameterSubvectorIndex_name_;
  static const int ParameterSubvectorIndex_default_;

  static const std::string ObservationMultiplier_name_;
  static const double ObservationMultiplier_default_;

  static const std::string ObservationTargetVector_name_;

  static const std::string ParameterMultiplier_name_;
  static const double ParameterMultiplier_default_;

  static const std::string ParameterBaseVector_name_;

  // ////////////////////////////////
  // Private member functions

  void initializeDefaults();

};

// /////////////////////////////////
// Implementations

// Static data members

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ObservationIndex_name_
= "Observation Index";

template<class Scalar>
const int
DefaultInverseModelEvaluator<Scalar>::ObservationIndex_default_
= -1;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterSubvectorIndex_name_
= "Parameter Subvector Index";

template<class Scalar>
const int
DefaultInverseModelEvaluator<Scalar>::ParameterSubvectorIndex_default_
= 0;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ObservationMultiplier_name_
= "Observation Multiplier";

template<class Scalar>
const double
DefaultInverseModelEvaluator<Scalar>::ObservationMultiplier_default_
= 1.0;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ObservationTargetVector_name_
= "Observation Target Vector";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterMultiplier_name_
= "Parameter Multiplier";

template<class Scalar>
const double
DefaultInverseModelEvaluator<Scalar>::ParameterMultiplier_default_
= 1e-6;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterBaseVector_name_
= "Parameter Base Vector";

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultInverseModelEvaluator<Scalar>::DefaultInverseModelEvaluator()
  :obs_idx_(-1),p_idx_(0)
{}

template<class Scalar>
DefaultInverseModelEvaluator<Scalar>::DefaultInverseModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &thyraModel
  )
  :obs_idx_(-1),p_idx_(0)
{
  initialize(thyraModel);
}

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &thyraModel
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  Ng_ = thyraModel->Ng()+1;
  inv_g_space_= thyraModel->get_x_space()->smallVecSpcFcty()->createVecSpc(1);
}

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::uninitialize(
  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >  *thyraModel
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
}

// Overridden from Teuchos::ParameterListAcceptor

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  using Teuchos::Array;
  using Teuchos::getParameterPtr;
  using Teuchos::rcp;
  using Teuchos::sublist;
  TEST_FOR_EXCEPT(0==paramList.get());
  paramList->validateParameters(*getValidParameters(),0); // Just validate my params
  paramList_ = paramList;
  obs_idx_ = paramList_->get(
    ObservationIndex_name_,ObservationIndex_default_);
  p_idx_ = paramList_->get(
    ParameterSubvectorIndex_name_,ParameterSubvectorIndex_default_);
  observationMultiplier_ = paramList_->get(
    ObservationMultiplier_name_,ObservationMultiplier_default_);
  if(get_observationTargetIO().get()) {
    observationTargetReader_.set_vecSpc(
      obs_idx_ < 0 ? this->get_x_space() : this->get_g_space(obs_idx_)
      );
    Teuchos::VerboseObjectTempState<ParameterDrivenMultiVectorInput<Scalar> >
      vots_observationTargetReader(
        rcp(&observationTargetReader_,false)
        ,this->getOStream(),this->getVerbLevel()
        );
    observationTargetReader_.setParameterList(
      sublist(paramList_,ObservationTargetVector_name_)
      );
    observationTargetReader_.readVector(
      "observation target vector",&observationTarget_);
  }
  parameterMultiplier_ = paramList_->get(
    ParameterMultiplier_name_,ParameterMultiplier_default_);
  if(get_parameterBaseIO().get()) {
    parameterBaseReader_.set_vecSpc(this->get_p_space(p_idx_));
    Teuchos::VerboseObjectTempState<ParameterDrivenMultiVectorInput<Scalar> >
      vots_parameterBaseReader(
        rcp(&parameterBaseReader_,false)
        ,this->getOStream(),this->getVerbLevel()
        );
    parameterBaseReader_.setParameterList(
      sublist(paramList_,ParameterBaseVector_name_)
      );
    parameterBaseReader_.readVector(
      "parameter base vector",&parameterBase_);
  }
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getParameterList()
{
  return paramList_;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getValidParameters() const
{
  if(validParamList_.get()==NULL) {
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set(
      ObservationIndex_name_,ObservationIndex_default_
      ,"The index of the observation function, obs_idx.\n"
      "If obs_idx < 0, then the observation will be the state vector x.\n"
      "If obs_idx >= 0, then the observation will be the response function g(obs_idx)."
      );
    pl->set(
      ParameterSubvectorIndex_name_,ParameterSubvectorIndex_default_
      ,"The index of the parameter subvector that will be used in the\n"
      "regularization term."
      );
    pl->set(
      ObservationMultiplier_name_,ObservationMultiplier_default_
      ,"observationMultiplier"
      );
    if(this->get_observationTargetIO().get())
      observationTargetReader_.set_fileIO(this->get_observationTargetIO());
    pl->sublist(ObservationTargetVector_name_).setParameters(
      *observationTargetReader_.getValidParameters()
      );
    pl->set(
      ParameterMultiplier_name_,ParameterMultiplier_default_
      ,"parameterMultiplier"
      );
    if(this->get_parameterBaseIO().get())
      parameterBaseReader_.set_fileIO(this->get_parameterBaseIO());
    pl->sublist(ParameterBaseVector_name_).setParameters(
      *parameterBaseReader_.getValidParameters()
      );
    validParamList_ = pl;
  }
  return validParamList_;
}

// Overridden from ModelEvaulator.

template<class Scalar>
int DefaultInverseModelEvaluator<Scalar>::Ng() const
{
  return Ng_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultInverseModelEvaluator<Scalar>::get_g_space(int j) const
{
  if(j==Ng_-1)
    return inv_g_space_;
  return this->getUnderlyingModel()->get_g_space(j);
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultInverseModelEvaluator<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np,Ng+1);
  outArgs.setSupports(wrappedOutArgs);
  outArgs.setSupports(MEB::OUT_ARG_DgDx,Ng,MEB::DERIV_TRANS_MV_BY_ROW);
  outArgs.setSupports(MEB::OUT_ARG_DgDp,Ng,p_idx_,MEB::DERIV_TRANS_MV_BY_ROW);
  return outArgs;
}

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::evalModel(
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

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultInverseModelEvaluator",inArgs,outArgs
    );

  //
  // See what needs to be computed
  //

  VectorBase<Scalar>
    *g_inv_out = outArgs.get_g(outArgs.Ng()-1).get();
  MultiVectorBase<Scalar>
    *DgDx_inv_trans_out = get_mv(
      outArgs.get_DgDx(outArgs.Ng()-1),"DgDx",MEB::DERIV_TRANS_MV_BY_ROW
      ).get();
  MultiVectorBase<Scalar>
    *DgDp_inv_trans_out = get_mv(
      outArgs.get_DgDp(outArgs.Ng()-1,p_idx_),"DgDp",MEB::DERIV_TRANS_MV_BY_ROW
      ).get();
  
  //
  // Compute all of the desired functions in the base model
  //

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the base point ...\n";

  MEB::InArgs<Scalar>  wrappedInArgs = inArgs;
  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  wrappedOutArgs.setArgs(outArgs,true);
  RefCountPtr<VectorBase<Scalar> >    wrapped_o;
  MEB::DerivativeMultiVector<Scalar>  wrapped_DoDx_trans, wrapped_DoDp_trans;
  if( obs_idx_ >= 0 && ( g_inv_out || DgDx_inv_trans_out || DgDp_inv_trans_out ) )
  {
    wrapped_o = createMember(thyraModel->get_g_space(obs_idx_));
    wrappedOutArgs.set_g(obs_idx_,wrapped_o);
    if( DgDx_inv_trans_out ) {
      wrapped_DoDx_trans = create_DgDx_mv(
        *thyraModel, obs_idx_, MEB::DERIV_TRANS_MV_BY_ROW
        );
      wrappedOutArgs.set_DgDx(obs_idx_,wrapped_DoDx_trans);
    }
    if( DgDp_inv_trans_out ) {
      wrapped_DoDp_trans = create_DgDp_mv(
        *thyraModel, obs_idx_, p_idx_, MEB::DERIV_TRANS_MV_BY_ROW
        );
      wrappedOutArgs.set_DgDp(obs_idx_,p_idx_,wrapped_DoDp_trans);
    }
  }
  
  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
  
  bool failed = wrappedOutArgs.isFailed();
  
  if(!failed) {
    //
    // Compute the inverse response function and its derivatives if asked to
    // do so.
    //
    const VectorBase<Scalar>
      &x = *inArgs.get_x(),
      &p = *inArgs.get_p(p_idx_);
    //
    //
    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
      o_space = ( obs_idx_ >= 0 ? this->get_g_space(obs_idx_) : this->get_x_space() ),
      p_space = this->get_p_space(p_idx_);
    const Index
      no = o_space->dim(),
      np = p_space->dim();
    Teuchos::RefCountPtr<VectorBase<Scalar> > diff_o, diff_p;
    if( g_inv_out || DgDx_inv_trans_out ) {
      const VectorBase<Scalar>
        &o = ( obs_idx_ < 0 ? x : *wrapped_o );
      diff_o = createMember(o_space);
      V_VmV( &*diff_o, o, *observationTarget_ );
    }
    if( g_inv_out || DgDp_inv_trans_out ) {
      diff_p = createMember(p_space);
      V_VmV( &*diff_p, p, *parameterBase_ );
    }
    if(g_inv_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function g(Np-1) ...\n";
      const Scalar
        observationTerm
        = ( observationMultiplier_ != ST::zero()
            ? observationMultiplier_*(0.5/no)*dot(*diff_o,*diff_o)
            : ST::zero()
          ),
        parameterTerm
        = ( parameterMultiplier_ != ST::zero()
            ? parameterMultiplier_*(0.5/np)*dot(*diff_p,*diff_p)
            : ST::zero()
          );
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out
          << "\nObservation matching term of g(Np-1):"
          << "\n  observationMultiplier*(0.5/no)*dot(diff_o,diff_o) = "
          << observationTerm
          << "\nParameter regularization term of g(Np-1):"
          << "\n  parameterMultiplier*(0.5/np)*dot(diff_p,diff_p) = "
          << parameterTerm
          << "\n";
      set_ele(0,observationTerm+parameterTerm,g_inv_out);
    }
    if(DgDx_inv_trans_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function derivative DgDx_trans ...\n";
      if( obs_idx_ < 0 ) {
        // DgDx^T = observationMultiplier*(1/no)*diff_x
        V_StV( &*DgDx_inv_trans_out->col(0), Scalar(observationMultiplier_*(1.0/no)), *diff_o );
      }
      else {
        // DgDx^T = (observationMultiplier*(1/no)) * (DoDx^T) * diff_o
        apply(
          *wrapped_DoDx_trans.getMultiVector(), NOTRANS
          ,*diff_o
          ,&*DgDx_inv_trans_out->col(0)
          ,Scalar(observationMultiplier_*(1.0/no))
          );
      }
    }
    if(DgDp_inv_trans_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function derivative DgDp_trans ...\n";
      if( obs_idx_ < 0 ) {
        // DgDp^T = parameterMultiplier*(1/np)*diff_p
        V_StV(
          &*DgDp_inv_trans_out->col(0)
          ,Scalar(parameterMultiplier_*(1.0/np))
          ,*diff_p
          );
      }
      else {
        // DgDp^T = (observationMultiplier*(1/no)) * (DoDp^T) * diff_o
        //          + (parameterMultiplier*(1.0/np)) * diff_p
        apply(
          *wrapped_DoDp_trans.getMultiVector(), NOTRANS
          ,*diff_o
          ,&*DgDp_inv_trans_out->col(0)
          ,Scalar(observationMultiplier_*(1.0/no))
          );
        if( parameterMultiplier_ != ST::zero() ) {
          Vp_StV(
            &*DgDp_inv_trans_out->col(0)
            ,Scalar(parameterMultiplier_*(1.0/np))
            ,*diff_p
            );
        }
      }
    }
  }
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultInverseModelEvaluator<Scalar>::description() const
{
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultInverseModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

// private

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::initializeDefaults()
{
  obs_idx_ = ObservationIndex_default_;
  p_idx_ = ParameterSubvectorIndex_default_;
  Ng_ = 0;
  observationMultiplier_ = ObservationMultiplier_default_;
  parameterMultiplier_ = ParameterMultiplier_default_; 
}

} // namespace Thyra

#endif // THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP
