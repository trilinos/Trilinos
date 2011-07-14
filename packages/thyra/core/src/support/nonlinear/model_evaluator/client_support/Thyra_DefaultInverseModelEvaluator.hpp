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
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Time.hpp"


namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and adds a simple, but
 * fairly general, inverse response function.
 *
 * The following response function is added to the end of the supported
 * response functions:

 \verbatim

  g_inv(x,p)
    = g_(getUnderlyingModel()->Ng())(x,p,...)
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

  d(g_inv)/d(x) = observationMultiplier * d(observationMatch)/d(x)

  d(g_inv)/d(p) = observationMultiplier * d(observationMatch)/d(p)
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


  d(g_inv)/d(x)^T = observationMultiplier * d(observationMatch)/d(x)^T


  d(g_inv)/d(p)^T = observationMultiplier * d(observationMatch)/d(p)^T
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
 * Another feature supported by this class is the ability to tack on parameter
 * regularization to an existing response function.  This mode is enabled by
 * setting the parameter "Observation Pass Through" to <tt>true</tt>.  This
 * results in the observation matching term to be defined as:

 \verbatim

   observationMatch(x,p) = o(x,p)

 \endverbatim

 * and has the derivatives:

 \verbatim

   d(observationMatch)/d(x)^T = d(o)/d(x)^T


   d(observationMatch)/d(p)^T = d(o)/d(p)^ T

 \endverbatim

 * Everything else about the above discussion.
 * 
 * <b>Note:</b> In this case, of course, the observation response function
 * must have dimension 1.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
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
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &thyraModel
    );

  /** \brief . */
  void uninitialize(
    RCP<ModelEvaluator<Scalar> > *thyraModel
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
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief .
   *
   * Note that <tt>observationTargetIO()</tt> and <tt>parameterBaseIO()</tt>
   * must be set before calling this function in order to have the sublists
   * added that will allow the vectors <tt>observationTarget()</tt> and
   * <tt>parameterBase()</tt> to be read in latter when the parameter list is
   * set..
   */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  // ////////////////////////////////
  // Private data members

  mutable RCP<const Teuchos::ParameterList> validParamList_;
  RCP<Teuchos::ParameterList>  paramList_;

  RCP<const VectorSpaceBase<Scalar> > inv_g_space_;

  mutable ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  mutable ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;
  mutable bool usingObservationTargetAsParameter_;

  int obs_idx_;
  int p_idx_;


  double observationMultiplier_;
  double parameterMultiplier_; 

  bool observationTargetAsParameter_;
  
  bool observationPassThrough_;

  Teuchos::EVerbosityLevel localVerbLevel_;

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

  static const std::string ObservationPassThrough_name_;
  static const bool ObservationPassThrough_default_;

  static const std::string ParameterMultiplier_name_;
  static const double ParameterMultiplier_default_;

  static const std::string ParameterBaseVector_name_;

  // ////////////////////////////////
  // Private member functions

  void initializeDefaults();

  void initializeInArgsOutArgs() const;

  RCP<const VectorSpaceBase<Scalar> > get_obs_space() const;

};


/** \brief Non-member constructor.
 *
 * \relates DefaultInverseModelEvaluator
 */
template<class Scalar>
RCP<DefaultInverseModelEvaluator<Scalar> >
defaultInverseModelEvaluator(
  const RCP<ModelEvaluator<Scalar> > &thyraModel
  )
{
  RCP<DefaultInverseModelEvaluator<Scalar> >
    inverseModel = Teuchos::rcp(new DefaultInverseModelEvaluator<Scalar>);
  inverseModel->initialize(thyraModel);
  return inverseModel;
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
= "Parameter Subvector Ordinal";

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
DefaultInverseModelEvaluator<Scalar>::ObservationPassThrough_name_
= "Observation Pass Through";

template<class Scalar>
const bool
DefaultInverseModelEvaluator<Scalar>::ObservationPassThrough_default_
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
  :usingObservationTargetAsParameter_(false), obs_idx_(-1),p_idx_(0),
   observationTargetAsParameter_(false),
   observationPassThrough_(ObservationPassThrough_default_),
   localVerbLevel_(Teuchos::VERB_DEFAULT)
{}


template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> > &thyraModel
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  inv_g_space_= thyraModel->get_x_space()->smallVecSpcFcty()->createVecSpc(1);
  // Get ready for reinitalization
  prototypeInArgs_ = ModelEvaluatorBase::InArgs<Scalar>();
  prototypeOutArgs_ = ModelEvaluatorBase::OutArgs<Scalar>();
}


template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::uninitialize(
  RCP<ModelEvaluator<Scalar> >  *thyraModel
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
}


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{

  using Teuchos::Array;
  using Teuchos::getParameterPtr;
  using Teuchos::rcp;
  using Teuchos::sublist;

  // Validate and set the parameter list
  TEST_FOR_EXCEPT(0==paramList.get());
  paramList->validateParameters(*getValidParameters(),0);
  paramList_ = paramList;

  // Parameters for observation matching term
  obs_idx_ = paramList_->get(
    ObservationIndex_name_,ObservationIndex_default_);
  observationPassThrough_ = paramList_->get(
    ObservationPassThrough_name_, ObservationPassThrough_default_ );
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      ( obs_idx_ < 0 &&  observationPassThrough_ ), std::logic_error,
      "Error, the observation function index obs_idx = " << obs_idx_ << " is not\n"
      "allowed when the observation is simply passed through!"
      );
#endif
  observationMultiplier_ = paramList_->get(
    ObservationMultiplier_name_,ObservationMultiplier_default_);
  if (!ObservationPassThrough_default_) {
    observationTargetAsParameter_ = paramList_->get(
      ObservationTargetAsParameter_name_, ObservationTargetAsParameter_default_ );
    if(get_observationTargetIO().get()) {
      observationTargetReader_.set_vecSpc(get_obs_space());
      Teuchos::VerboseObjectTempState<ParameterDrivenMultiVectorInput<Scalar> >
        vots_observationTargetReader(
          rcp(&observationTargetReader_,false)
          ,this->getOStream(),this->getVerbLevel()
          );
      observationTargetReader_.setParameterList(
        sublist(paramList_,ObservationTargetVector_name_)
        );
      RCP<VectorBase<Scalar> >
        observationTarget;
      observationTargetReader_.readVector(
        "observation target vector",&observationTarget);
      observationTarget_ = observationTarget;
    }
  }
  else {
    observationTargetAsParameter_ = false;
    observationTarget_ = Teuchos::null;
  }
  
  // Parameters for parameter matching term
  p_idx_ = paramList_->get(
    ParameterSubvectorIndex_name_,ParameterSubvectorIndex_default_);
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
    RCP<VectorBase<Scalar> >
      parameterBase;
    parameterBaseReader_.readVector(
      "parameter base vector",&parameterBase);
    parameterBase_ = parameterBase;
  }

  // Verbosity settings
  localVerbLevel_ = this->readLocalVerbosityLevelValidatedParameter(*paramList_);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);

#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG

  // Get ready for reinitalization
  prototypeInArgs_ = ModelEvaluatorBase::InArgs<Scalar>();
  prototypeOutArgs_ = ModelEvaluatorBase::OutArgs<Scalar>();

}


template<class Scalar>
RCP<Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getValidParameters() const
{
  if(validParamList_.get()==NULL) {
    RCP<Teuchos::ParameterList>
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
    pl->set( ObservationPassThrough_name_, ObservationPassThrough_default_,
      "If true, then the observation will just be used instead of the least-squares\n"
      "function.  This allows you to add a parameter regularization term to any existing\n"
      "response function!"
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
    this->setLocalVerbosityLevelValidatedParameter(&*pl);
    Teuchos::setupVerboseObjectSublist(&*pl);
    validParamList_ = pl;
  }
  return validParamList_;
}


// Overridden from ModelEvaulator.


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultInverseModelEvaluator<Scalar>::get_p_space(int l) const
{
  if (prototypeInArgs_.Np()==0)
    initializeInArgsOutArgs();
  if ( l == prototypeInArgs_.Np()-1 && usingObservationTargetAsParameter_ )
    return get_obs_space();
  return this->getUnderlyingModel()->get_p_space(l);
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultInverseModelEvaluator<Scalar>::get_g_space(int j) const
{
  if (prototypeOutArgs_.Np()==0)
    initializeInArgsOutArgs();
  if (j==prototypeOutArgs_.Ng()-1)
    return inv_g_space_;
  return this->getUnderlyingModel()->get_g_space(j);
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultInverseModelEvaluator<Scalar>::createInArgs() const
{
  if (prototypeInArgs_.Np()==0)
    initializeInArgsOutArgs();
  return prototypeInArgs_;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultInverseModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
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


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultInverseModelEvaluator<Scalar>::createOutArgsImpl() const
{
  if (prototypeOutArgs_.Np()==0)
    initializeInArgsOutArgs();
  return prototypeOutArgs_;
}


template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using std::endl;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ModelEvaluatorBase MEB;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_LOCALVERBLEVEL_BEGIN(
    "Thyra::DefaultInverseModelEvaluator",inArgs,outArgs,localVerbLevel_
    );

  const bool trace = out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_LOW);
  const bool print_p = out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_MEDIUM);
  const bool print_x = out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_EXTREME);
  const bool print_o = print_x;

  //
  // A) See what needs to be computed
  //
  
  const RCP<VectorBase<Scalar> >
    g_inv_out = outArgs.get_g(outArgs.Ng()-1);
  const RCP<MultiVectorBase<Scalar> >
    DgDx_inv_trans_out = get_mv(
      outArgs.get_DgDx(outArgs.Ng()-1), "DgDx", MEB::DERIV_TRANS_MV_BY_ROW
      );
  const RCP<MultiVectorBase<Scalar> >
    DgDp_inv_trans_out = get_mv(
      outArgs.get_DgDp(outArgs.Ng()-1,p_idx_), "DgDp", MEB::DERIV_TRANS_MV_BY_ROW
      );
  const bool computeInverseFunction = ( nonnull(g_inv_out)
    || nonnull(DgDx_inv_trans_out) || nonnull(DgDp_inv_trans_out) );
  
  //
  // B) Compute all of the needed functions from the base model
  //

  if(trace)
    *out << "\nComputing the base point and the observation(s) ...\n";

  MEB::InArgs<Scalar>  wrappedInArgs = thyraModel->createInArgs();
  wrappedInArgs.setArgs(inArgs,true);
  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  wrappedOutArgs.setArgs(outArgs,true);
  RCP<VectorBase<Scalar> > wrapped_o;
  MEB::Derivative<Scalar> wrapped_DoDx;
  MEB::Derivative<Scalar> wrapped_DoDp_trans;
  if( obs_idx_ >= 0 && computeInverseFunction )
  {
    wrapped_o = createMember(thyraModel->get_g_space(obs_idx_));
    wrappedOutArgs.set_g(obs_idx_,wrapped_o);
    if (nonnull(DgDx_inv_trans_out)) {
      if (!observationPassThrough_)
        wrapped_DoDx = thyraModel->create_DgDx_op(obs_idx_);
      else
        wrapped_DoDx = Thyra::create_DgDx_mv(
          *thyraModel, obs_idx_, MEB::DERIV_TRANS_MV_BY_ROW );
      wrappedOutArgs.set_DgDx(obs_idx_,wrapped_DoDx);
    }
    if (nonnull(DgDp_inv_trans_out)) {
      wrapped_DoDp_trans = create_DgDp_mv(
        *thyraModel, obs_idx_, p_idx_, MEB::DERIV_TRANS_MV_BY_ROW
        );
      wrappedOutArgs.set_DgDp(obs_idx_,p_idx_,wrapped_DoDp_trans);
    }
    // 2007/07/28: rabartl: Above, we really should check if these output
    // arguments have already been set by the client.  If they are, then we
    // need to make sure that they are of the correct form or we need to throw
    // an exception!
  }

  if (!wrappedOutArgs.isEmpty()) {
    thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
  }
  else {
    if(trace)
      *out << "\nSkipping the evaluation of the underlying model since "
           << "there is nothing to compute ...\n";
  }
  
  bool failed = wrappedOutArgs.isFailed();

  //
  // C) Assemble the final observation and paramter terms
  //
  
  if ( !failed && computeInverseFunction ) {

    //
    // Compute the inverse response function and its derivatives
    //

    RCP<const VectorBase<Scalar> >
      x_in = inArgs.get_x(),
      p_in = inArgs.get_p(p_idx_);

    const MEB::InArgs<Scalar> nominalValues = this->getNominalValues();
    RCP<const VectorBase<Scalar> >
      x = ( !is_null(x_in) ? x_in : nominalValues.get_x().assert_not_null() ),
      p = ( !is_null(p_in) ? p_in : nominalValues.get_p(p_idx_).assert_not_null() );

    const RCP<const VectorSpaceBase<Scalar> >
      o_space = get_obs_space(),
      p_space = this->get_p_space(p_idx_);

    const Ordinal
      no = o_space->dim(),
      np = p_space->dim();
    
    if (trace)
      *out << "\nno = " << no
           << "\nnp = " << np
           << endl;

#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      observationPassThrough_ && no != 1, std::logic_error,
      "Error, the observation function dimension no="<<no<<" > 1 is not allowed"
      " when the observation is passed through as the observation matching term!"
      );
#endif

    // Compute diff_o if needed
    RCP<const VectorBase<Scalar> > o;
    RCP<VectorBase<Scalar> > diff_o;
    if( !observationPassThrough_ && ( nonnull(g_inv_out) || nonnull(DgDx_inv_trans_out) ) )
    {
      if (obs_idx_ < 0 ) o = x; else o = wrapped_o; // can't use ( test ? x : wrapped_o )!
      if(trace) *out << "\n||o||inf = " << norm_inf(*o) << endl;
      if (print_o) *out << "\no = " << *o;
      diff_o = createMember(o_space);
      RCP<const VectorBase<Scalar> >
        observationTarget
        = ( observationTargetAsParameter_
          ? inArgs.get_p(inArgs.Np()-1)
          : Teuchos::null
          );
      if (is_null(observationTarget) ) {
        observationTarget = observationTarget_;
        if (trace)
          *out << "\n||ot||inf = " << norm_inf(*observationTarget) << endl;
        if (print_o)
          *out << "\not = " << *observationTarget;
      }
      if (!is_null(observationTarget)) {
        V_VmV( diff_o.ptr(), *o, *observationTarget );
      }
      else {
        assign( diff_o.ptr(), *o );
      }
      if(trace) {
        *out << "\n||diff_o||inf = " << norm_inf(*diff_o) << endl;
      }
      if (print_o) {
        *out << "\ndiff_o = " << *diff_o;
      }
    }
  
    // Compute diff_p if needed
    RCP<VectorBase<Scalar> > diff_p;
    if ( nonnull(g_inv_out) || nonnull(DgDp_inv_trans_out)) {
      if(trace) *out << "\n||p||inf = " << norm_inf(*p) << endl;
      if(print_p) *out << "\np = " << Teuchos::describe(*p,Teuchos::VERB_EXTREME);
      diff_p = createMember(p_space);
      if (!is_null(parameterBase_) ) {
        if(trace) *out << "\n||pt||inf = " << norm_inf(*parameterBase_) << endl;
        if(print_p) {
          *out << "\npt = "
               << Teuchos::describe(*parameterBase_,Teuchos::VERB_EXTREME);
        }
        V_VmV( diff_p.ptr(), *p, *parameterBase_ );
      }
      else {
        assign( diff_p.ptr(), *p );
      }
      if(trace) {*out << "\n||diff_p|| = " << norm(*diff_p) << endl;}
      if(print_p) {
        *out << "\ndiff_p = "
             << Teuchos::describe(*diff_p, Teuchos::VERB_EXTREME);
      }
    }
    

    // Get and check Q_o and Q_p

    RCP<const LinearOpBase<Scalar> >
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

    // Compute Q_o * diff_o
    RCP<VectorBase<Scalar> > Q_o_diff_o;
    if ( !is_null(Q_o) && !is_null(diff_o) ) {
      Q_o_diff_o = createMember(Q_o->range()); // Should be same as domain!
      apply( *Q_o, NOTRANS, *diff_o, Q_o_diff_o.ptr() );
    }
    
    // Compute Q_p * diff_p
    RCP<VectorBase<Scalar> > Q_p_diff_p;
    if ( !is_null(Q_p)  && !is_null(diff_p)  ) {
      Q_p_diff_p = createMember(Q_p->range()); // Should be same as domain!
      apply( *Q_p, NOTRANS, *diff_p, Q_p_diff_p.ptr() );
    }

    // Compute g_inv(x,p)
    if (nonnull(g_inv_out)) {
      if(trace)
        *out << "\nComputing inverse response function ginv = g(Np-1) ...\n";
      const Scalar observationTerm
        = ( observationPassThrough_
          ? get_ele(*wrapped_o,0) // ToDo; Verify that this is already a scalar
          : ( observationMultiplier_ != ST::zero()
            ? ( !is_null(Q_o)
              ?  observationMultiplier_*0.5*dot(*diff_o,*Q_o_diff_o)
              : observationMultiplier_*(0.5/no)*dot(*diff_o,*diff_o)
              )
            : ST::zero()
            )
          );
      const Scalar parameterTerm
        = ( parameterMultiplier_ != ST::zero()
          ? ( !is_null(Q_p)
            ?  parameterMultiplier_*0.5*dot(*diff_p,*Q_p_diff_p)
            : parameterMultiplier_*(0.5/np)*dot(*diff_p,*diff_p)
            )
          : ST::zero()
          );
      const Scalar g_inv_val = observationTerm+parameterTerm;
      if(trace)
        *out
          << "\nObservation matching term of ginv = g(Np-1):"
          << "\n  observationMultiplier = " << observationMultiplier_
          << "\n  observationMultiplier*observationMatch(x,p) = " << observationTerm
          << "\nParameter regularization term of ginv = g(Np-1):"
          << "\n  parameterMultiplier = " << parameterMultiplier_
          << "\n  parameterMultiplier*parameterRegularization(p) = " << parameterTerm
          << "\nginv = " << g_inv_val
          << "\n";
      set_ele(0, observationTerm+parameterTerm, g_inv_out.ptr());
    }

    // Compute d(g_inv)/d(x)^T
    if (nonnull(DgDx_inv_trans_out)) {
      if(trace)
        *out << "\nComputing inverse response function derivative DginvDx^T:\n";
      if (!observationPassThrough_) {
        if( obs_idx_ < 0 ) {
          if (!is_null(Q_o)) {
            if (trace)
              *out << "\nDginvDx^T = observationMultiplier * Q_o * diff_o ...\n";
            V_StV(
              DgDx_inv_trans_out->col(0).ptr(),
              observationMultiplier_,
              *Q_o_diff_o
              );
          }
          else {
            if (trace)
              *out << "\nDginvDx^T = observationMultiplier * (1/no) * diff_o ...\n";
            V_StV(
              DgDx_inv_trans_out->col(0).ptr(),
              Scalar(observationMultiplier_*(1.0/no)),
              *diff_o
              );
          }
        }
        else {
          //if (trace)
          //  *out << "\n||DoDx^T||inf = " << norms_inf(*wrapped_DoDx.getMultiVector()) << endl;
          if (print_o && print_x)
            *out << "\nDoDx = " << *wrapped_DoDx.getLinearOp();
          if (!is_null(Q_o)) {
            if (trace)
              *out << "\nDginvDx^T = observationMultiplier * DoDx^T * Q_o * diff_o ...\n";
            apply(
              *wrapped_DoDx.getLinearOp(), CONJTRANS,
              *Q_o_diff_o,
              DgDx_inv_trans_out->col(0).ptr(),
              observationMultiplier_
              );
          }
          else {
            if (trace)
              *out << "\nDginvDx^T = (observationMultiplier*(1/no)) * DoDx^T * diff_o ...\n";
            apply(
              *wrapped_DoDx.getLinearOp(), CONJTRANS,
              *diff_o,
              DgDx_inv_trans_out->col(0).ptr(),
              Scalar(observationMultiplier_*(1.0/no))
              );
          }
        }
      }
      else {
        if (trace)
          *out << "\nDginvDx^T = observationMultiplier * DoDx^T ...\n";
        V_StV(
          DgDx_inv_trans_out->col(0).ptr(), observationMultiplier_,
          *wrapped_DoDx.getMultiVector()->col(0)
          );
      }
      if(trace)
        *out << "\n||DginvDx^T||inf = " << norms_inf(*DgDx_inv_trans_out) << "\n";
      if (print_x)
        *out << "\nDginvDx^T = " << *DgDx_inv_trans_out;
    }

    // Compute d(g_inv)/d(p)^T
    if (nonnull(DgDp_inv_trans_out)) {
      if(trace)
        *out << "\nComputing inverse response function derivative DginvDp^T ...\n";
      if (obs_idx_ >= 0) {
        if (trace)
          *out << "\n||DoDp^T|| = " << norms_inf(*wrapped_DoDp_trans.getMultiVector()) << endl;
        if (print_p)
          *out << "\nDoDp^T = " << Teuchos::describe(*wrapped_DoDp_trans.getMultiVector(),Teuchos::VERB_EXTREME);
      }
      if(trace)
        *out << "\nDginvDp^T = 0 ...\n";
      assign( DgDp_inv_trans_out->col(0).ptr(), ST::zero() );
      // DgDp^T += observationMultiplier * d(observationMatch)/d(p)^T
      if (!observationPassThrough_) {
        if ( obs_idx_ >= 0 ) {
          if ( !is_null(Q_o) ) {
            if(trace)
              *out << "\nDginvDp^T += observationMultiplier* * (DoDp^T) * Q_o * diff_o ...\n";
            apply(
              *wrapped_DoDp_trans.getMultiVector(), NOTRANS,
              *Q_o_diff_o,
              DgDp_inv_trans_out->col(0).ptr(),
              Scalar(observationMultiplier_*(1.0/no)),
              ST::one()
              );
          }
          else {
            if(trace)
              *out << "\nDgDp^T += observationMultiplier* * (DoDp^T) * Q_o * diff_o ...\n";
            apply(
              *wrapped_DoDp_trans.getMultiVector(), NOTRANS,
              *diff_o,
              DgDp_inv_trans_out->col(0).ptr(),
              Scalar(observationMultiplier_*(1.0/no)),
              ST::one()
              );
          }
          if(trace)
            *out << "\n||DginvDp^T||inf = " << norms_inf(*DgDp_inv_trans_out) << "\n";
          if (print_p)
            *out << "\nDginvDp^T = " << *DgDp_inv_trans_out;
        }
        else {
          // d(observationMatch)/d(p)^T = 0, nothing to do!
        }
      }
      else {
        if(trace)
          *out << "\nDginvDp^T += (observationMultiplier*(1/no)) * (DoDp^T) * diff_o ...\n";
        Vp_StV(
          DgDp_inv_trans_out->col(0).ptr(), observationMultiplier_,
          *wrapped_DoDp_trans.getMultiVector()->col(0)
          );
        
      }
      // DgDp^T += parameterMultiplier * d(parameterRegularization)/d(p)^T
      if( parameterMultiplier_ != ST::zero() ) {
        if ( !is_null(Q_p) ) {
          if(trace)
            *out << "\nDginvDp^T += parameterMultiplier * Q_p * diff_p ...\n";
          Vp_StV(
            DgDp_inv_trans_out->col(0).ptr(),
            parameterMultiplier_,
            *Q_p_diff_p
            );
        }
        else {
          if(trace)
            *out << "\nDginvDp^T += (parameterMultiplier*(1.0/np)) * diff_p ...\n";
          Vp_StV(
            DgDp_inv_trans_out->col(0).ptr(),
            Scalar(parameterMultiplier_*(1.0/np)),
            *diff_p
            );
        }
        if(trace)
          *out << "\n||DginvDp^T||inf = " << norms_inf(*DgDp_inv_trans_out) << "\n";
        if (print_p)
        *out << "\nDginvDp^T = " << *DgDp_inv_trans_out;
      }
      else {
        // This term is zero so there is nothing to do!
      }
    }

  }
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


// private


template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::initializeDefaults()
{
  obs_idx_ = ObservationIndex_default_;
  p_idx_ = ParameterSubvectorIndex_default_;
  observationMultiplier_ = ObservationMultiplier_default_;
  parameterMultiplier_ = ParameterMultiplier_default_; 
}


template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::initializeInArgsOutArgs() const
{

  typedef ModelEvaluatorBase MEB;

  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  const MEB::InArgs<Scalar> wrappedInArgs = thyraModel->createInArgs();
  const int wrapped_Np = wrappedInArgs.Np();

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  const bool supports_x = wrappedInArgs.supports(MEB::IN_ARG_x);
  usingObservationTargetAsParameter_ = ( supports_x && observationTargetAsParameter_ );
  inArgs.setSupports(
    wrappedInArgs,
    wrapped_Np + ( usingObservationTargetAsParameter_ ? 1 : 0 )
    );
  prototypeInArgs_ = inArgs;

  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int wrapped_Ng = wrappedOutArgs.Ng();

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(inArgs.modelEvalDescription());
  outArgs.set_Np_Ng( prototypeInArgs_.Np(), wrapped_Ng+1 );
  outArgs.setSupports(wrappedOutArgs);
  outArgs.setSupports(MEB::OUT_ARG_DgDx,wrapped_Ng,MEB::DERIV_TRANS_MV_BY_ROW);
  outArgs.setSupports(MEB::OUT_ARG_DgDp,wrapped_Ng,p_idx_,MEB::DERIV_TRANS_MV_BY_ROW);
  prototypeOutArgs_ = outArgs;
  
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultInverseModelEvaluator<Scalar>::get_obs_space() const
{
  return ( obs_idx_ < 0 ? this->get_x_space() : this->get_g_space(obs_idx_) );
}


} // namespace Thyra


#endif // THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP
