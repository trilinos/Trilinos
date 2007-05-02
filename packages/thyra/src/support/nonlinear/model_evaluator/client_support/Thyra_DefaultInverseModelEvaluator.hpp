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
#include "Thyra_AssertOp.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object and adds a simple, but
 * fairly general, inverse-type of response function.
 *
 * The following response function is added to the end of the supported
 * response functions:

 \verbatim

  g_(getUnderlyingModel()->Ng())(x,p,...)
    = observationMultiplier * observationMatch(x,p)
    + parameterMultiplier * parameterRegularization(p)

 \endverbatim

 * where <tt>observationMatch(x,p)</tt> is some scalar-valued function that
 * gives the match of some state observation,
 * <tt>parameterRegularization(p)</tt> is some scaled valued function that
 * regularizes the parameters, and <tt>observationMultiplier</tt> and
 * <tt>parameterMultiplier</tt> are scalar constant multipliers for the state
 * observation and the parameter regularization respectively.
 *
 * The state observation matching function and the parameter regularization
 * function can be defined in one of two ways.
 *
 * If a symmetric positive definite linear operator <tt>Q_o</tt> is defined,
 * then the state observation matching function is given as:

 \verbatim

   observationMatch(x,p) = 0.5 * diff_o(x,p)^T * Q_o * diff_o(x,p)

 \endverbatim

 * and if <tt>Q_o</tt> is not defined, then the state observation matching
 * function is given as:

 \verbatim

   observationMatch(x,p) = (0.5/no) * diff_o(x,p)^T * diff_o(x,p)

 \endverbatim

 * where

 \verbagtim

   diff_o(x,p) = o(x,p) - ot

 \endverbatim

 * where <tt>ot</tt> is the target vector for some observation (see below) and
 * <tt>p</tt> is one of the parameter subvectors supported by the underlying
 * model.
 *
 * The observation function <tt>o(x,p)</tt> can be the state vector itself
 * <tt>o(x,p) = x</tt> for <tt>obs_idx < 0</tt>, or can be any of the built-in
 * response functions <tt>o(x,p) = g(obs_idx)(x,p)</tt> when <tt>0 <= obs_idx
 * < getUnderlyingModel()->Ng()</tt>.
 *
 * The parameter regularization function also has one of two definitions.
 *
 * If a symmetric positive definite linear operator <tt>Q_p</tt> is defined,
 * then the parameter regularization function is given as:

 \verbatim

   parameterRegularization(p) = 0.5 * diff_p(p)^T * Q_p * diff_p(p)

 \endverbatim

 * and if <tt>Q_p</tt> is not defined, then the parameter regularization
 * function is given as:

 \verbatim

   parameterRegularization(p) = (0.5/np) * diff_p(p)^T * diff_p(p)

 \endverbatim

 * where

 \verbagtim

   diff_p(p) = p - pt

 \endverbatim

 * where <tt>pt</tt> is a nomial parameter vector for which violations are
 * penalized against.
 *
 * Since this decorator class adds a response function, then <tt>this->Ng() ==
 * getUnderlyingModel()->Ng() + 1</tt>.
 *
 * Let's consider the derivatives of this inverse function.
 *
 * The first derivatives are given by:

 \verbatim

  d(g)/d(x) = observationMultiplier * d(observationMatch)/d(x)

  d(g)/d(p) = observationMultiplier * d(observationMatch)/d(p)
            + parameterMultiplier * d(parameterRegularization)/d(p)

 \endverbatim
 
 * where the derivatives of <tt>observationMatch(x,p)</tt> and
 * <tt>parameterRegularization(p)</tt> are given by:

 \verbatim


                              /  diff_o(x,p)^T * Q_o * d(o)/d(x) : Q_o defined
   d(observationMatch)/d(x) = |
                              \  (1/no) * diff_o(x,p)^T * d(o)/d(x) : Q_o not defined


                              /  diff_o(x,p)^T * Q_o * d(o)/d(p) : Q_o defined
   d(observationMatch)/d(p) = |
                              \  (1/no) * diff_o(x,p)^T * d(o)/d(p) : Q_o not defined


                                     /  diff_p(p)^T * Q_p : Q_p defined
   d(parameterRegularization)/d(p) = |
                                     \  (1/np) * diff_p(p)^T : Q_p not defined


 \endverbatim

 * Of course when <tt>obs_idx < -1</tt> where <tt>o(x,p) = x</tt> then
 * <tt>d(o)/d(x) = I</tt> and <tt>d(o)/d(p) = 0</tt> which also gives
 * <tt>d(observationMatch)/d(p) = 0</tt>.
 *
 * Also, we typically want these derivatives in gradient form which gives:


 \verbatim


  d(g)/d(x)^T = observationMultiplier * d(observationMatch)/d(x)^T


  d(g)/d(p)^T = observationMultiplier * d(observationMatch)/d(p)^T
              + parameterMultiplier * d(parameterRegularization)/d(p)^T


                                /  d(o)/d(x)^T * Q_o * diff_o(x,p) : Q_o defined
   d(observationMatch)/d(x)^T = |
                                \  (1/no) * d(o)/d(x)^T * diff_o(x,p) : Q_o not defined


                                /  d(o)/d(p)^T * Q_o * diff_o(x,p) : Q_o defined
   d(observationMatch)/d(p)^T = |
                                \  (1/no) * d(o)/d(p)^T * diff_o(x,p) : Q_o not defined


                                       /  Q_p * diff_p(p) : Q_p defined
   d(parameterRegularization)/d(p)^T = |
                                       \  (1/np) * diff_p(p) : Q_p not defined


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
  STANDARD_CONST_COMPOSITION_MEMBERS( VectorBase<Scalar>, observationTarget );

  /** \brief Parameter base vector <tt>pt</tt>. */
  STANDARD_CONST_COMPOSITION_MEMBERS( VectorBase<Scalar>, parameterBase );

  /** \brief Observation match weighting operator <tt>Q_o</tt>. */
  STANDARD_CONST_COMPOSITION_MEMBERS( LinearOpBase<Scalar>, observationMatchWeightingOp );

  /** \brief Parameter regulization weighting operator <tt>Q_p</tt>. */
  STANDARD_CONST_COMPOSITION_MEMBERS( LinearOpBase<Scalar>, parameterRegularizationWeightingOp );

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
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> > &thyraModel
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> > &thyraModel
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> > *thyraModel
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
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
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

  int obs_idx_;
  int p_idx_;
  int Ng_;

  double observationMultiplier_;
  double parameterMultiplier_; 

  bool observationTargetAsParameter_;

  mutable ParameterDrivenMultiVectorInput<Scalar> observationTargetReader_;
  mutable ParameterDrivenMultiVectorInput<Scalar> parameterBaseReader_;

  static const std::string ObservationIndex_name_;
  static const int ObservationIndex_default_;

  static const std::string ParameterSubvectorIndex_name_;
  static const int ParameterSubvectorIndex_default_;

  static const std::string ObservationMultiplier_name_;
  static const double ObservationMultiplier_default_;

  static const std::string ObservationTargetVector_name_;

  static const std::string ObservationTargetAsParameter_name_;
  static const bool ObservationTargetAsParameter_default_;

  static const std::string ParameterMultiplier_name_;
  static const double ParameterMultiplier_default_;

  static const std::string ParameterBaseVector_name_;

  // ////////////////////////////////
  // Private member functions

  void initializeDefaults();

};


/** \brief Non-member constructor.
 *
 * \relates DefaultInverseModelEvaluator
 */
template<class Scalar>
Teuchos::RefCountPtr<DefaultInverseModelEvaluator<Scalar> >
inverseModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> > &thyraModel
  )
{
  return Teuchos::rcp(
    new DefaultInverseModelEvaluator<Scalar>(thyraModel)
    );
}


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
DefaultInverseModelEvaluator<Scalar>::ObservationTargetAsParameter_name_
= "Observation Target as Parameter";

template<class Scalar>
const bool
DefaultInverseModelEvaluator<Scalar>::ObservationTargetAsParameter_default_
= false;

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
  :obs_idx_(-1),p_idx_(0), observationTargetAsParameter_(false)
{}

template<class Scalar>
DefaultInverseModelEvaluator<Scalar>::DefaultInverseModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &thyraModel
  )
  :obs_idx_(-1),p_idx_(0), observationTargetAsParameter_(false)
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
    Teuchos::RefCountPtr<VectorBase<Scalar> >
      observationTarget;
    observationTargetReader_.readVector(
      "observation target vector",&observationTarget);
    observationTarget_ = observationTarget;
  }
  observationTargetAsParameter_ = paramList_->get(
    ObservationTargetAsParameter_name_, ObservationTargetAsParameter_default_ );
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
    Teuchos::RefCountPtr<VectorBase<Scalar> >
      parameterBase;
    parameterBaseReader_.readVector(
      "parameter base vector",&parameterBase);
    parameterBase_ = parameterBase;
  }
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
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
    pl->set( ObservationIndex_name_,ObservationIndex_default_,
      "The index of the observation function, obs_idx.\n"
      "If obs_idx < 0, then the observation will be the state vector x.\n"
      "If obs_idx >= 0, then the observation will be the response function g(obs_idx)."
      );
    pl->set( ParameterSubvectorIndex_name_,ParameterSubvectorIndex_default_,
      "The index of the parameter subvector that will be used in the\n"
      "regularization term."
      );
    pl->set( ObservationMultiplier_name_,ObservationMultiplier_default_,
      "observationMultiplier"
      );
    if(this->get_observationTargetIO().get())
      observationTargetReader_.set_fileIO(this->get_observationTargetIO());
    pl->sublist(ObservationTargetVector_name_).setParameters(
      *observationTargetReader_.getValidParameters()
      );
    pl->set( ObservationTargetAsParameter_name_, ObservationTargetAsParameter_default_,
      "If true, then a parameter will be accepted for the state observation vector\n"
      "to allow it to be set by an external client through the InArgs object."
      );
    pl->set( ParameterMultiplier_name_,ParameterMultiplier_default_,
      "parameterMultiplier" );
    if(this->get_parameterBaseIO().get())
      parameterBaseReader_.set_fileIO(this->get_parameterBaseIO());
    pl->sublist(ParameterBaseVector_name_).setParameters(
      *parameterBaseReader_.getValidParameters()
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
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
ModelEvaluatorBase::InArgs<Scalar>
DefaultInverseModelEvaluator<Scalar>::createInArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::InArgs<Scalar> wrappedInArgs = thyraModel->createInArgs();
  const int wrapped_Np = wrappedInArgs.Np();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  const bool supports_x = wrappedInArgs.supports(MEB::IN_ARG_x);
  inArgs.setSupports(
    wrappedInArgs,
    wrapped_Np + ( supports_x && observationTargetAsParameter_ ? 1 : 0 )
    );
  return inArgs;
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
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
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

  MEB::InArgs<Scalar>  wrappedInArgs = thyraModel->createInArgs();
  wrappedInArgs.setArgs(inArgs,true);
  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  wrappedOutArgs.setArgs(outArgs,true);
  RefCountPtr<VectorBase<Scalar> > wrapped_o;
  MEB::DerivativeMultiVector<Scalar> wrapped_DoDx_trans, wrapped_DoDp_trans;
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

    RefCountPtr<const VectorBase<Scalar> >
      x_in = inArgs.get_x(),
      p_in = inArgs.get_p(p_idx_);

    const MEB::InArgs<Scalar> nominalValues = this->getNominalValues();
    RefCountPtr<const VectorBase<Scalar> >
      x = ( !is_null(x_in) ? x_in : nominalValues.get_x().assert_not_null() ),
      p = ( !is_null(p_in) ? p_in : nominalValues.get_p(p_idx_).assert_not_null() );

    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
      o_space = ( obs_idx_ >= 0 ? this->get_g_space(obs_idx_) : this->get_x_space() ),
      p_space = this->get_p_space(p_idx_);

    const Index
      no = o_space->dim(),
      np = p_space->dim();

    Teuchos::RefCountPtr<VectorBase<Scalar> > diff_o;
    if( g_inv_out || DgDx_inv_trans_out  ) {
      const VectorBase<Scalar>
        &o = ( obs_idx_ < 0 ? *x : *wrapped_o );
      diff_o = createMember(o_space);
      Teuchos::RefCountPtr<const VectorBase<Scalar> >
        observationTarget
        = ( observationTargetAsParameter_
          ? inArgs.get_p(inArgs.Np()-1)
          : Teuchos::null
          );
      if (is_null(observationTarget) )
        observationTarget = observationTarget_;
      if (!is_null(observationTarget)) {
        V_VmV( &*diff_o, o, *observationTarget );
      }
      else {
        assign( &*diff_o, o );
      }
    }
    
    Teuchos::RefCountPtr<VectorBase<Scalar> > diff_p;
    if( g_inv_out || DgDp_inv_trans_out ) {
      diff_p = createMember(p_space);
      if (!is_null(parameterBase_) ) {
        V_VmV( &*diff_p, *p, *parameterBase_ );
      }
      else {
        assign( &*diff_p, *p );
      }
    }
    
    Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
      Q_o = this->get_observationMatchWeightingOp(),
      Q_p = this->get_parameterRegularizationWeightingOp();

#ifdef TEUCHOS_DEBUG
    if (!is_null(Q_o)) {
      THYRA_ASSERT_VEC_SPACES(
        "Thyra::DefaultInverseModelEvaluator::evalModel(...)",
        *Q_o->range(), *o_space
        );
      THYRA_ASSERT_VEC_SPACES(
        "Thyra::DefaultInverseModelEvaluator::evalModel(...)",
        *Q_o->domain(), *o_space
        );
    }
    if (!is_null(Q_p)) {
      THYRA_ASSERT_VEC_SPACES(
        "Thyra::DefaultInverseModelEvaluator::evalModel(...)",
        *Q_p->range(), *p_space
        );
      THYRA_ASSERT_VEC_SPACES(
        "Thyra::DefaultInverseModelEvaluator::evalModel(...)",
        *Q_p->domain(), *p_space
        );
    }
    // Note, we have not proved that Q_o and Q_p are s.p.d. but at least we
    // have established that that have the right range and domain spaces!
#endif

    Teuchos::RefCountPtr<VectorBase<Scalar> > Q_o_diff_o;
    if ( !is_null(Q_o) && !is_null(diff_o) ) {
      Q_o_diff_o = createMember(Q_o->range()); // Should be same as domain!
      apply( *Q_o, NOTRANS, *diff_o, &*Q_o_diff_o );
    }

    Teuchos::RefCountPtr<VectorBase<Scalar> > Q_p_diff_p;
    if ( !is_null(Q_p)  && !is_null(diff_p)  ) {
      Q_p_diff_p = createMember(Q_p->range()); // Should be same as domain!
      apply( *Q_p, NOTRANS, *diff_p, &*Q_p_diff_p );
    }

    if(g_inv_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function g(Np-1) ...\n";
      const Scalar observationTerm
        = ( observationMultiplier_ != ST::zero()
          ? ( !is_null(Q_o)
            ?  observationMultiplier_*0.5*dot(*diff_o,*Q_o_diff_o)
            : observationMultiplier_*(0.5/no)*dot(*diff_o,*diff_o)
            )
          : ST::zero()
          );
      const Scalar parameterTerm
        = ( parameterMultiplier_ != ST::zero()
          ? ( !is_null(Q_p)
            ?  parameterMultiplier_*0.5*dot(*diff_p,*Q_p_diff_p)
            : parameterMultiplier_*(0.5/np)*dot(*diff_p,*diff_p)
            )
          : ST::zero()
          );
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out
          << "\nObservation matching term of g(Np-1):"
          << "\n  observationMultiplier*observationMatch(x,p) = "
          << observationTerm
          << "\nParameter regularization term of g(Np-1):"
          << "\n  parameterMultiplier*parameterRegularization(p) = "
          << parameterTerm
          << "\n";
      set_ele(0,observationTerm+parameterTerm,g_inv_out);
    }

    if(DgDx_inv_trans_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function derivative DgDx_trans ...\n";
      if( obs_idx_ < 0 ) {
        if (!is_null(Q_o)) {
          // DgDx^T = observationMultiplier * Q_o * diff_x
          V_StV(
            &*DgDx_inv_trans_out->col(0),
            observationMultiplier_,
            *Q_o_diff_o
            );
        }
        else {
          // DgDx^T = observationMultiplier * (1/no) * diff_x
          V_StV(
            &*DgDx_inv_trans_out->col(0),
            Scalar(observationMultiplier_*(1.0/no)),
            *diff_o
            );
        }
      }
      else {
        if (!is_null(Q_o)) {
          // DgDx^T = observationMultiplier * (DoDx^T) * Q_o * diff_o
          apply(
            *wrapped_DoDx_trans.getMultiVector(), NOTRANS,
            *Q_o_diff_o,
            &*DgDx_inv_trans_out->col(0),
            observationMultiplier_
            );
        }
        else {
          // DgDx^T = (observationMultiplier*(1/no)) * (DoDx^T) * diff_o
          apply(
            *wrapped_DoDx_trans.getMultiVector(), NOTRANS,
            *diff_o,
            &*DgDx_inv_trans_out->col(0),
            Scalar(observationMultiplier_*(1.0/no))
            );
        }
      }
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\n||DgDx_trans||inf = " << norm_inf(*DgDx_inv_trans_out->col(0)) << "\n";
    }

    if(DgDp_inv_trans_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function derivative DgDp_trans ...\n";
      // DgDp^T = 0
      assign( &*DgDp_inv_trans_out->col(0), ST::zero() );
      // DgDp^T += observationMultiplier * d(observationMatch)/d(p)^T
      if ( obs_idx_ >= 0 ) {
        if ( !is_null(Q_o) ) {
          // DgDp^T += observationMultiplier* * (DoDp^T) * Q_o * diff_o
          apply(
            *wrapped_DoDp_trans.getMultiVector(), NOTRANS,
            *Q_o_diff_o,
            &*DgDp_inv_trans_out->col(0),
            Scalar(observationMultiplier_*(1.0/no)),
            ST::one()
            );
        }
        else {
          // DgDp^T += (observationMultiplier*(1/no)) * (DoDp^T) * diff_o
          apply(
            *wrapped_DoDp_trans.getMultiVector(), NOTRANS,
            *diff_o,
            &*DgDp_inv_trans_out->col(0),
            Scalar(observationMultiplier_*(1.0/no)),
            ST::one()
            );
        }
      }
      else {
        // d(observationMatch)/d(p)^T = 0, nothing to do!
      }
      // DgDp^T += parameterMultiplier * d(parameterRegularization)/d(p)^T
      if( parameterMultiplier_ != ST::zero() ) {
        if ( !is_null(Q_p) ) {
          // DgDp^T += parameterMultiplier * Q_p * diff_p
          Vp_StV(
            &*DgDp_inv_trans_out->col(0),
            parameterMultiplier_,
            *Q_p_diff_p
            );
        }
        else {
          // DgDp^T += (parameterMultiplier*(1.0/np)) * diff_p
          Vp_StV(
            &*DgDp_inv_trans_out->col(0),
            Scalar(parameterMultiplier_*(1.0/np)),
            *diff_p
            );
        }
      }
      else {
        // This term is zero so there is nothing to do!
      }
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\n||DgDp_trans||inf = " << norm_inf(*DgDp_inv_trans_out->col(0)) << "\n";
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
