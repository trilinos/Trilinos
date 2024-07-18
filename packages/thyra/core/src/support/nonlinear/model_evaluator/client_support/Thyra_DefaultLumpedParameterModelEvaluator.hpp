// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAUL_LUMPED_PARAMETER_LUMPED_MODEL_EVALUATOR_HPP
#define THYRA_DEFAUL_LUMPED_PARAMETER_LUMPED_MODEL_EVALUATOR_HPP


#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

#include "sillyModifiedGramSchmidt.hpp" // This is just an example!


namespace Thyra {


/** \brief Decorator class that wraps any ModelEvaluator object and lumps
 * parameters together using a linear basis matrix.
 *
 * \section Thyra_DefaultLumpedParameterModelEvaluator_intro_sec Introduction
 *
 * The purpose of this Decorator class is to provide a linear basis reduction,
 * or "lumping", for a set of parameters.  Let <tt>p_orig</tt> be one of the
 * parameter subvectors of the underlying model
 * <tt>*getUnderlyingModel()</tt>.  This class provides an affine model for
 * the reduced parameters:
 
 \verbatim

   p_orig = B * p + p_orig_base

 \endverbatim

 * where <tt>B</tt> is some <tt>MultiVectorBase</tt> object with linearly
 * independent columns where
 * <tt>B.range()->isCompatible(*p_orig.space())==true</tt> and
 * <tt>B.domain()->dim() <= B.range()->dim()</tt>.  The basis matrix
 * <tt>B</tt> can be set by the client or can be generated automatically in
 * this object for any number of columns.  The vector <tt>p_orig_base</tt> is
 * generally selected to be the nominal values of the original parameters
 * <tt>p_orig</tt> in which case <tt>p == 0</tt> gives the original nominal
 * values and <tt>p != 0</tt> gives a perturbation from the nominal values.
 * This is the default option but a different base vector can be set.  Note
 * that storing <tt>B</tt> as a column-wise multi-vector allows the space
 * <tt>*p_orig.space()</tt> to be any arbitrary vector space, including a
 * parallel distributed vector space in an SPMD program.  It is this
 * flexibility that truly makes this decorator class an ANA class.
 *
 * The reduced parameter subvector <tt>p</tt> is what is exposed to an ANA
 * client dnd the original parameter subvector <tt>p_orig</tt> is what is
 * passed to the original underlying model.  Note that given <tt>p</tt>,
 * <tt>p_orig</tt> is uniquely determined but opposite is not true in general
 * (see Section \ref Thyra_DefaultLumpedParameterModelEvaluator_mapping_sec).
 *
 * \section Thyra_DefaultLumpedParameterModelEvaluator_derivs_sec Derivatives
 *
 * The reduced basis <tt>B</tt> of course affects the definition of all of the
 * derivatives with respect to the given parameter subvector.  The
 * relationship between the derivative of any function <tt>h</tt> with respect
 * to <tt>p_orig</tt> and <tt>p</tt> is:

 \verbatim

   d(h)/d(p) = d(h)/d(p_orig) * B

 \endverbatim

 * When <tt>d(h)/d(p)</tt> is only needed as a linear operator, both
 * <tt>d(h)/d(p_orig)</tt> and <tt>B</tt> would just need to be supplied as
 * general linear operators and then the multiplied linear operator
 * <tt>d(h)/d(p_orig)*B</tt> could be represented implicitly using a
 * <tt>Thyra::DefaultMultipliedLinearOp</tt> subclass object.
 *
 * When <tt>d(h)/d(p)</tt> is only needed as a column-wise multi-vector, then
 * <tt>d(h)/d(p_orig)</tt> is only needed as a general linear operator and
 * <tt>B</tt> must be a column-wise multi-vector.
 *
 * Lastly, when the row-wise transpose multi-vector form of
 * <tt>d(h)/d(p)^T</tt> is requested, it can be computed as:

 \verbatim

   d(h)/d(p)^T = B^T * d(h)/d(p_orig)^T

 \endverbatim

 * which requires that <tt>d(h)/d(p_orig)^T</tt> be supplied in row-wise
 * transpose multi-vector form but <tt>B</tt> could actually just be a general
 * linear operator.  This would allow for huge spaces for <tt>p_orig</tt> and
 * <tt>p</tt> which may be of some use for adjoint-based sensitivity methods.
 *
 * Since the current implementation in this class requires <tt>B</tt> be a
 * multi-vector, both forms <tt>d(h)/d(p)</tt> and <tt>d(h)/d(p)^T</tt> can be
 * supported.  In the future, it would be possible to support the
 * implementation of <tt>B</tt> as a general linear operator which would mean
 * that very large full and reduced parameter spaces could be supported.
 * However, focusing on small reduced parameter spaces is the very motivation
 * for this decorator subclass and therefore requiring that <tt>B</tt> be
 * represented as a multi-vector is not a serious limitation.
 *
 * \section Thyra_DefaultLumpedParameterModelEvaluator_mapping_sec Mapping between fully and reduced parameters
 *
 * In cases where <tt>p_orig</tt> is given and <tt>p</tt> must be determined,
 * in general there is no solution for <tt>p</tt> that will satisfy
 * <tt>p_orig=B*p+p_orig_base</tt>.
 *
 * To support the (overdetermined) mapping from <tt>p_orig</tt> to <tt>p</tt>,
 * we choose <tt>p</tt> to solve the classic linear least squares problem:
 
 \verbatim

   min   0.5 * (B*p + p_orig_base - p_orig)^T * (B*p + p_orig_base - p_orig)

 \endverbatim

 * This well known linear least squares problem has the solution:
 
 \verbatim

   p = inv(B^T * B) * (B^T * (-p_orig_base+p_oirg))

 \endverbatim
 
 * This approach has the unfortunate side effect that we can not completely
 * represent an arbitrary vector <tt>p_orig</tt> with a reduced <tt>p</tt> but
 * this is the nature of this approximation for both good and bad.
 *
 * The one case in which a unique reduced <tt>p</tt> can be determined is when
 * <tt>p_orig</tt> was computed from the afffine relationship given <tt>p</tt>
 * and therefore we expect a zero residual meaning that
 * <tt>(p_orig_base-p_oirg)</tt> lies in the null-space of <tt>B</tt>.  This
 * is handy when one wants to write out <tt>p</tt> in the form of
 * <tt>p_orig</tt> and then reconstruct <tt>p</tt> again.
 *
 * \section Thyra_DefaultLumpedParameterModelEvaluator_bounds_sec Parameter bounds
 *
 * As mentioned above, one can simply select <tt>p_orig_base</tt> to be the
 * nominal values of the original full parameters <tt>p_orig</tt> and
 * therefore allow <tt>p==0</tt> to give the exact nominal parameter values
 * which is critical for the initial guess for many ANAs.  The problem of
 * handling the bounds on the parameters is more difficult.  This
 * approximation really requires that the simple bounds on <tt>p_orig</tt> be
 * replaced with the general linear inequality constraints:

 \verbatim

   p_orig_l - p_orig_base <= B*p <= p_orig_u - p_orig_base

 \endverbatim

 * Modeling and enforcing these linear inequality constraints correctly would
 * require that <tt>B</tt> and <tt>p_orig_base</tt> be exposed to the client
 * so that the client can build these extra <em>linear</em> inequality
 * constraints into the ANA algorithm.  Certain optimization algorithms could
 * then enforce these general inequalities and therefore guarantee that the
 * original parameter bounds are never violated.  This is easy to do in both
 * active-set and interior-point methods when the size of the space
 * <tt>p_orig.space()->dim()</tt> is not too large.  When
 * <tt>p_orig.space()->dim()</tt> is large, handling these inequality
 * constraints can be very difficult to deal with.
 *
 * An alternative simpler approach would be to try to find lower and upper
 * bounds <tt>p_l</tt> and <tt>p_u</tt> that tried to match the original
 * bounds as close as possible while not restricting the feasible region
 * <tt>p_l <= p <= p_u</tt> too much. For example, one could solve the
 * inequality-constrained least squares problems:

 \verbatim

   min     0.5 * (B*p+p_orig_base-p_orig_l)^T * (B*p+p_orig_base-p_orig_l)
   s.t.    B*p+p_orig_base >= p_orig_l

 \endverbatim

 and

 \verbatim

   min     0.5 * (B*p+p_orig_base-p_orig_u)^T * (B*p+p_orig_base-p_orig_u)
   s.t.    B*p+p_orig_base <= p_orig_u

 \endverbatim

 * but in general it would not be possible to guarantee that the solution
 * <tt>p_l</tt> and <tt>p_u</tt> satisfies <tt>p_l <= p <= p_u</tt>.
 *
 * Because of the challenges of dealing with bounds on the parameters, this
 * subclass currently throws an exception if it is given an underlying model
 * with finite bounds on the parameters.  However, a parameter-list option can
 * be set that will cause the bounds to be ignored and it would be the
 * client's responsibility to deal with the implications of this choice.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultLumpedParameterModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
  , virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultLumpedParameterModelEvaluator();

  /** \brief . */
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &thyraModel
    );

  /** \brief . */
  void uninitialize(
    RCP<ModelEvaluator<Scalar> > *thyraModel
    );

  // 2007/07/30: rabartl: ToDo: Add functions to set and get the underlying
  // basis matrix!

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief .  */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    );

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

  mutable bool isInitialized_;
  mutable bool nominalValuesAndBoundsUpdated_;

  mutable RCP<const Teuchos::ParameterList> validParamList_;
  RCP<Teuchos::ParameterList> paramList_;

  // Parameters read from the parameter list
  int p_idx_;
  bool autogenerateBasisMatrix_;
  int numberOfBasisColumns_;
  bool nominalValueIsParameterBase_;
  bool ignoreParameterBounds_;
  Teuchos::EVerbosityLevel localVerbLevel_;
  bool dumpBasisMatrix_;

  // Reduced affine parameter model
  mutable RCP<const MultiVectorBase<Scalar> > B_;
  mutable RCP<const VectorBase<Scalar> > p_orig_base_;

  // Nominal values and bounds
  mutable ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  mutable ModelEvaluatorBase::InArgs<Scalar> lowerBounds_;
  mutable ModelEvaluatorBase::InArgs<Scalar> upperBounds_;

  // Static

  static const std::string ParameterSubvectorIndex_name_;
  static const int ParameterSubvectorIndex_default_;

  static const std::string AutogenerateBasisMatrix_name_;
  static const bool AutogenerateBasisMatrix_default_;

  static const std::string NumberOfBasisColumns_name_;
  static const int NumberOfBasisColumns_default_;

  static const std::string NominalValueIsParameterBase_name_;
  static const bool NominalValueIsParameterBase_default_;

  static const std::string ParameterBaseVector_name_;

  static const std::string IgnoreParameterBounds_name_;
  static const bool IgnoreParameterBounds_default_;

  static const std::string DumpBasisMatrix_name_;
  static const bool DumpBasisMatrix_default_;

  // ////////////////////////////////
  // Private member functions

  // These functions are used to implement late initialization so that the
  // need for clients to order function calls is reduced.

  // Finish enough initialization to defined spaces etc.
  void finishInitialization() const;

  // Generate the parameter basis matrix B.
  void generateParameterBasisMatrix() const;

  // Finish all of initialization needed to define nominal values, bounds, and
  // p_orig_base.  This calls finishInitialization().
  void updateNominalValuesAndBounds() const;

  // Map from p -> p_orig.
  RCP<VectorBase<Scalar> >
  map_from_p_to_p_orig( const VectorBase<Scalar> &p ) const;

  // Set up the arguments for DhDp_orig to be computed by the underlying model.
  void setupWrappedParamDerivOutArgs(
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs, // in
    ModelEvaluatorBase::OutArgs<Scalar> *wrappedOutArgs // in/out
    ) const;

  // Create DhDp_orig needed to assembled DhDp
  ModelEvaluatorBase::Derivative<Scalar>
  create_deriv_wrt_p_orig(
    const ModelEvaluatorBase::Derivative<Scalar> &DhDp,
    const ModelEvaluatorBase::EDerivativeMultiVectorOrientation requiredOrientation
    ) const;

  // After DhDp_orig has been computed, assemble DhDp or DhDp^T for all deriv
  // output arguments.
  void assembleParamDerivOutArgs(
    const ModelEvaluatorBase::OutArgs<Scalar> &wrappedOutArgs, // in
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs // in/out
    ) const;

  // Given a single DhDp_orig, assemble DhDp
  void assembleParamDeriv(
    const ModelEvaluatorBase::Derivative<Scalar> &DhDp_orig, // in
    const ModelEvaluatorBase::Derivative<Scalar> &DhDp // in/out
    ) const;

};


/** \brief Non-member constructor.
 *
 * \relates DefaultLumpedParameterModelEvaluator
 */
template<class Scalar>
RCP<DefaultLumpedParameterModelEvaluator<Scalar> >
defaultLumpedParameterModelEvaluator(
  const RCP<ModelEvaluator<Scalar> > &thyraModel
  )
{
  RCP<DefaultLumpedParameterModelEvaluator<Scalar> >
    paramLumpedModel = Teuchos::rcp(new DefaultLumpedParameterModelEvaluator<Scalar>);
  paramLumpedModel->initialize(thyraModel);
  return paramLumpedModel;
}


// /////////////////////////////////
// Implementations


// Static data members


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::ParameterSubvectorIndex_name_
= "Parameter Subvector Index";

template<class Scalar>
const int
DefaultLumpedParameterModelEvaluator<Scalar>::ParameterSubvectorIndex_default_
= 0;


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::AutogenerateBasisMatrix_name_
= "Auto-generate Basis Matrix";

template<class Scalar>
const bool
DefaultLumpedParameterModelEvaluator<Scalar>::AutogenerateBasisMatrix_default_
= true;


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::NumberOfBasisColumns_name_
= "Number of Basis Columns";

template<class Scalar>
const int
DefaultLumpedParameterModelEvaluator<Scalar>::NumberOfBasisColumns_default_
= 1;


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::NominalValueIsParameterBase_name_
= "Nominal Value is Parameter Base";

template<class Scalar>
const bool
DefaultLumpedParameterModelEvaluator<Scalar>::NominalValueIsParameterBase_default_
= true;


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::ParameterBaseVector_name_
= "Parameter Base Vector";


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::IgnoreParameterBounds_name_
= "Ignore Parameter Bounds";

template<class Scalar>
const bool
DefaultLumpedParameterModelEvaluator<Scalar>::IgnoreParameterBounds_default_
= false;


template<class Scalar>
const std::string
DefaultLumpedParameterModelEvaluator<Scalar>::DumpBasisMatrix_name_
= "Dump Basis Matrix";

template<class Scalar>
const bool
DefaultLumpedParameterModelEvaluator<Scalar>::DumpBasisMatrix_default_
= false;


// Constructors/initializers/accessors/utilities


template<class Scalar>
DefaultLumpedParameterModelEvaluator<Scalar>::DefaultLumpedParameterModelEvaluator()
  :isInitialized_(false),
   nominalValuesAndBoundsUpdated_(false),
   p_idx_(ParameterSubvectorIndex_default_),
   autogenerateBasisMatrix_(AutogenerateBasisMatrix_default_),
   numberOfBasisColumns_(NumberOfBasisColumns_default_),
   nominalValueIsParameterBase_(NominalValueIsParameterBase_default_),
   ignoreParameterBounds_(IgnoreParameterBounds_default_),
   localVerbLevel_(Teuchos::VERB_DEFAULT),
   dumpBasisMatrix_(DumpBasisMatrix_default_)
{}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> > &thyraModel
  )
{
  isInitialized_ = false;
  nominalValuesAndBoundsUpdated_ = false;
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::uninitialize(
  RCP<ModelEvaluator<Scalar> > *thyraModel
  )
{
  isInitialized_ = false;
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string
DefaultLumpedParameterModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultLumpedParameterModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{

  using Teuchos::getParameterPtr;
  using Teuchos::rcp;
  using Teuchos::sublist;

  isInitialized_ = false;
  nominalValuesAndBoundsUpdated_ = false;

  // Validate and set the parameter list
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters(),0);
  paramList_ = paramList;

  // Read in parameters
  p_idx_ = paramList_->get(
    ParameterSubvectorIndex_name_, ParameterSubvectorIndex_default_ );
  autogenerateBasisMatrix_ = paramList_->get(
    AutogenerateBasisMatrix_name_, AutogenerateBasisMatrix_default_ );
  if (autogenerateBasisMatrix_) {
    numberOfBasisColumns_ = paramList_->get(
      NumberOfBasisColumns_name_, NumberOfBasisColumns_default_ );
  }
  nominalValueIsParameterBase_ = paramList_->get(
    NominalValueIsParameterBase_name_, NominalValueIsParameterBase_default_ );
  if (!nominalValueIsParameterBase_) {
    TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement reading parameter base vector from file!");
  }
  ignoreParameterBounds_ = paramList_->get(
    IgnoreParameterBounds_name_, IgnoreParameterBounds_default_ );
  dumpBasisMatrix_ = paramList_->get(
    DumpBasisMatrix_name_, DumpBasisMatrix_default_ );

  // Verbosity settings
  localVerbLevel_ = this->readLocalVerbosityLevelValidatedParameter(*paramList_);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);

#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif

}


template<class Scalar>
RCP<Teuchos::ParameterList>
DefaultLumpedParameterModelEvaluator<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<Teuchos::ParameterList>
DefaultLumpedParameterModelEvaluator<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
DefaultLumpedParameterModelEvaluator<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
DefaultLumpedParameterModelEvaluator<Scalar>::getValidParameters() const
{
  if(validParamList_.get()==NULL) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set( ParameterSubvectorIndex_name_, ParameterSubvectorIndex_default_,
      "Determines the index of the parameter subvector in the underlying model\n"
      "for which the reduced basis representation will be determined." );
    pl->set( AutogenerateBasisMatrix_name_, AutogenerateBasisMatrix_default_,
      "If true, then a basis matrix will be auto-generated for a given number\n"
      " of basis vectors." );
    pl->set( NumberOfBasisColumns_name_, NumberOfBasisColumns_default_,
      "If a basis is auto-generated, then this parameter gives the number\n"
      "of columns in the basis matrix that will be created.  Warning!  This\n"
      "number must be less than or equal to the number of original parameters\n"
      "or an exception will be thrown!" );
    pl->set( NominalValueIsParameterBase_name_, NominalValueIsParameterBase_default_,
      "If true, then the nominal values for the full parameter subvector from the\n"
      "underlying model will be used for p_orig_base.  This allows p==0 to give\n"
      "the nominal values for the parameters." );
    /*
    if(this->get_parameterBaseIO().get())
      parameterBaseReader_.set_fileIO(this->get_parameterBaseIO());
    pl->sublist(ParameterBaseVector_name_).setParameters(
      *parameterBaseReader_.getValidParameters()
      );
    */
    pl->set( IgnoreParameterBounds_name_, IgnoreParameterBounds_default_,
      "If true, then any bounds on the parameter subvector will be ignored." );
    pl->set( DumpBasisMatrix_name_, DumpBasisMatrix_default_,
      "If true, then the basis matrix will be printed the first time it is created\n"
      "as part of the verbose output and as part of the Describable::describe(...)\n"
      "output for any verbositiy level >= \"low\"." );
    this->setLocalVerbosityLevelValidatedParameter(&*pl);
    Teuchos::setupVerboseObjectSublist(&*pl);
    validParamList_ = pl;
  }
  return validParamList_;
}


// Overridden from ModelEvaulator.


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultLumpedParameterModelEvaluator<Scalar>::get_p_space(int l) const
{
  finishInitialization();
  if (l == p_idx_)
    return B_->domain();
  return this->getUnderlyingModel()->get_p_space(l);
}


template<class Scalar>
RCP<const Array<std::string> >
DefaultLumpedParameterModelEvaluator<Scalar>::get_p_names(int l) const
{
  finishInitialization();
  if (l == p_idx_)
    return Teuchos::null; // Names for these parameters would be meaningless!
  return this->getUnderlyingModel()->get_p_names(l);
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultLumpedParameterModelEvaluator<Scalar>::getNominalValues() const
{
  updateNominalValuesAndBounds();
  return nominalValues_;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultLumpedParameterModelEvaluator<Scalar>::getLowerBounds() const
{
  updateNominalValuesAndBounds();
  return lowerBounds_;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultLumpedParameterModelEvaluator<Scalar>::getUpperBounds() const
{
  updateNominalValuesAndBounds();
  return upperBounds_;
}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
  const bool wasSolved
  )
{

  typedef ModelEvaluatorBase MEB;

  // Make sure that everything has been initialized
  updateNominalValuesAndBounds();

  const RCP<ModelEvaluator<Scalar> >
    thyraModel = this->getNonconstUnderlyingModel();

  // By default, copy all input arguments since they will all be the same
  // except for the given reduced p.  We will then replace the reduced p with
  // p_orig below.
  MEB::InArgs<Scalar> wrappedFinalPoint = thyraModel->createInArgs();
  wrappedFinalPoint.setArgs(finalPoint);

  // Replace p with p_orig.
  RCP<const VectorBase<Scalar> > p;
  if (!is_null(p=finalPoint.get_p(p_idx_))) {
    wrappedFinalPoint.set_p(p_idx_, map_from_p_to_p_orig(*p));
  }

  thyraModel->reportFinalPoint(wrappedFinalPoint,wasSolved);

}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultLumpedParameterModelEvaluator<Scalar>::createOutArgsImpl() const
{
  ModelEvaluatorBase::OutArgsSetup<Scalar>
    outArgs = this->getUnderlyingModel()->createOutArgs();
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
  // 2007/07/31: rabartl: ToDo: We need to manually set the forms of the
  // derivatives that this class object will support!  This needs to be based
  // on tests of what the forms of derivatives the underlying model supports.
}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  // This routine is pretty simple for the most part.  By default, we just
  // pass everything through to the underlying model evaluator except for
  // arguments reated to the parameter subvector with index
  // p_idx_.

  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ModelEvaluatorBase MEB;

  // Make sure that everything has been initialized
  updateNominalValuesAndBounds();

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_LOCALVERBLEVEL_BEGIN(
    "Thyra::DefaultLumpedParameterModelEvaluator",inArgs,outArgs,localVerbLevel_
    );

  //
  // A) Setup InArgs
  //

  // By default, copy all input arguments since they will all be the same
  // except for the given reduced p.  We will then replace the reduced p with
  // p_orig below.
  MEB::InArgs<Scalar> wrappedInArgs = thyraModel->createInArgs();
  wrappedInArgs.setArgs(inArgs);
  
  // Replace p with p_orig.
  RCP<const VectorBase<Scalar> > p;
  if (!is_null(p=wrappedInArgs.get_p(p_idx_))) {
    if (
      dumpBasisMatrix_
      && includesVerbLevel(localVerbLevel,Teuchos::VERB_MEDIUM)
      )
    {
      *out << "\nB = " << Teuchos::describe(*B_,Teuchos::VERB_EXTREME);
    }
    wrappedInArgs.set_p(p_idx_,map_from_p_to_p_orig(*p));
  }

  //
  // B) Setup OutArgs
  //

  // By default, copy all output arguments since they will all be the same
  // except for those derivatives w.r.t. p(p_idx).  We will then replace the
  // derivative objects w.r.t. given reduced p with the derivarive objects
  // w.r.t. p_orig below.
  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  wrappedOutArgs.setArgs(outArgs);

  // Set derivative output arguments for p_orig if derivatives for p are
  // reqeusted in outArgs
  setupWrappedParamDerivOutArgs(outArgs,&wrappedOutArgs);

  //
  // C) Evaluate the underlying model functions
  //

  if (includesVerbLevel(localVerbLevel,Teuchos::VERB_LOW))
    *out << "\nEvaluating the fully parameterized underlying model ...\n";
  // Compute the underlying functions in terms of p_orig, including
  // derivatives w.r.t. p_orig.
  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);

  //
  // D) Postprocess the output arguments
  //

  // Assemble the derivatives for p given derivatives for p_orig computed
  // above.
  assembleParamDerivOutArgs(wrappedOutArgs,outArgs);
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


// private


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::finishInitialization() const
{

  typedef ScalarTraits<Scalar> ST;
  typedef ModelEvaluatorBase MEB;

  if (isInitialized_)
    return;

  //
  // A) Get the underlying model
  //

  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  TEUCHOS_TEST_FOR_EXCEPTION(
    is_null(thyraModel), std::logic_error,
    "Error, the underlying model evaluator must be set!" );

  //
  // B) Create B for the reduced affine model for the given parameter subvector
  //

  if (autogenerateBasisMatrix_) {
    generateParameterBasisMatrix();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't handle a client-set parameter basis matrix yet!" );
  }

  isInitialized_ = true;

}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::generateParameterBasisMatrix() const
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  const RCP<const VectorSpaceBase<Scalar> >
    p_orig_space = thyraModel->get_p_space(p_idx_);

  const Ordinal p_orig_dim = p_orig_space->dim();

  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 1 <= numberOfBasisColumns_ && numberOfBasisColumns_ <= p_orig_dim ),
    std::logic_error,
    "Error, the number of basis columns = " << numberOfBasisColumns_ << " does not\n"
    "fall in the range [1,"<<p_orig_dim<<"]!" );

  //
  // Create and randomize B
  //
  // Here we make the first column all ones and then randomize columns 1
  // through numberOfBasisColumns_-1 so that the average entry is 1.0 with a
  // spread of 1.0.  This is just to give as well a scaled starting matrix as
  // possible that will hopefully yeild a well scaled orthonomal B after we
  // are finished.

  const RCP<MultiVectorBase<Scalar> >
    B = createMembers(p_orig_space,numberOfBasisColumns_);
  assign( B->col(0).ptr(), ST::one() );
  if (numberOfBasisColumns_ > 1) {
    seed_randomize<double>(0);
    Thyra::randomize( as<Scalar>(0.5*ST::one()), as<Scalar>(1.5*ST::one()),
      B.ptr() );
  }

  //
  // Create an orthogonal form of B using a modified Gram Schmidt algorithm
  //

  RCP<MultiVectorBase<double> > R;
  sillyModifiedGramSchmidt( B.ptr(), Teuchos::outArg(R) );

  // Above:
  // 1) On output, B will have orthonomal columns which makes it a good basis
  // 2) We just discard the "R" factor since we don't need it for anything 

  B_ = B;

}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::updateNominalValuesAndBounds() const
{

  typedef ScalarTraits<Scalar> ST;
  typedef ModelEvaluatorBase MEB;

  if (nominalValuesAndBoundsUpdated_)
    return;

  finishInitialization();

  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  const MEB::InArgs<Scalar> origNominalValues = thyraModel->getNominalValues();
  const MEB::InArgs<Scalar> origLowerBounds = thyraModel->getLowerBounds();
  const MEB::InArgs<Scalar> origUpperBounds = thyraModel->getUpperBounds();

  // p_orig_base

  if (nominalValueIsParameterBase_) {
    const RCP<const VectorBase<Scalar> >
      p_orig_init = origNominalValues.get_p(p_idx_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      is_null(p_orig_init), std::logic_error,
      "Error, if the user requested that the nominal values be used\n"
      "as the base vector p_orig_base then that vector has to exist!" );
    p_orig_base_ = p_orig_init->clone_v();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't handle reading in the parameter base vector yet!" );
  }

  // Nominal values

  nominalValues_ = origNominalValues;
  
  if (nominalValueIsParameterBase_) {
    // A value of p==0 will give p_orig = p_orig_init!
    const RCP<VectorBase<Scalar> >
      p_init = createMember(B_->domain());
    assign( p_init.ptr(), ST::zero() );
    nominalValues_.set_p(p_idx_, p_init);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't handle creating p_init when p_orig_base != p_orig_init yet!" );
  }

  // Bounds

  lowerBounds_ = origLowerBounds;
  upperBounds_ = origUpperBounds;

  lowerBounds_.set_p(p_idx_,Teuchos::null);
  upperBounds_.set_p(p_idx_,Teuchos::null);

  if (!ignoreParameterBounds_) {
    const RCP<const VectorBase<Scalar> >
      p_orig_l = origLowerBounds.get_p(p_idx_),
      p_orig_u = origUpperBounds.get_p(p_idx_);
    if ( !is_null(p_orig_l) || !is_null(p_orig_u) ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error, we don't handle bounds on p_orig yet!" );
    }
  }

  nominalValuesAndBoundsUpdated_ = true;

}


template<class Scalar>
RCP<VectorBase<Scalar> >
DefaultLumpedParameterModelEvaluator<Scalar>::map_from_p_to_p_orig(
  const VectorBase<Scalar> &p
  ) const
{
  // p_orig = B*p + p_orig_base
  const RCP<VectorBase<Scalar> > p_orig = createMember(B_->range());
  apply( *B_, NOTRANS, p, p_orig.ptr() );
  Vp_V( p_orig.ptr(), *p_orig_base_ );
  return p_orig;
}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::setupWrappedParamDerivOutArgs(
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs, // in
  ModelEvaluatorBase::OutArgs<Scalar> *wrappedOutArgs_inout // in/out
  ) const
{

  typedef ModelEvaluatorBase MEB;
  typedef MEB::Derivative<Scalar> Deriv;

  TEUCHOS_TEST_FOR_EXCEPT(wrappedOutArgs_inout==0);
  MEB::OutArgs<Scalar> &wrappedOutArgs = *wrappedOutArgs_inout;
    
  Deriv DfDp;
  if ( !(DfDp=outArgs.get_DfDp(p_idx_)).isEmpty() ) {
    wrappedOutArgs.set_DfDp(p_idx_,create_deriv_wrt_p_orig(DfDp,MEB::DERIV_MV_BY_COL));
  }

  const int Ng = outArgs.Ng();
  for ( int j = 0; j < Ng; ++j ) {
    Deriv DgDp;
    if ( !(DgDp=outArgs.get_DgDp(j,p_idx_)).isEmpty() ) {
      wrappedOutArgs.set_DgDp(
        j, p_idx_,
        create_deriv_wrt_p_orig(DgDp,DgDp.getMultiVectorOrientation())
        );
    }
  }

}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
DefaultLumpedParameterModelEvaluator<Scalar>::create_deriv_wrt_p_orig(
  const ModelEvaluatorBase::Derivative<Scalar> &DhDp,
  const ModelEvaluatorBase::EDerivativeMultiVectorOrientation requiredOrientation
  ) const
{

  typedef ModelEvaluatorBase MEB;

  const RCP<const MultiVectorBase<Scalar> >
    DhDp_mv = DhDp.getMultiVector();
  TEUCHOS_TEST_FOR_EXCEPTION(
    is_null(DhDp_mv) || (DhDp.getMultiVectorOrientation() != requiredOrientation),
    std::logic_error,
    "Error, we currently can't handle non-multi-vector derivatives!" );

  RCP<MultiVectorBase<Scalar> > DhDp_orig_mv;
  switch (requiredOrientation) {
    case MEB::DERIV_MV_BY_COL:
      // DhDp = DhDp_orig * B
      DhDp_orig_mv = createMembers(DhDp_mv->range(),B_->range()->dim());
      // Above, we could just request DhDp_orig as a LinearOpBase object since
      // we just need to apply it!
      break;
    case MEB::DERIV_TRANS_MV_BY_ROW:
      // (DhDp^T) = B^T * (DhDp_orig^T)  [DhDp_orig_mv is transposed!]
      DhDp_orig_mv = createMembers(B_->range(),DhDp_mv->domain()->dim());
      // Above, we really do need DhDp_orig as the gradient form multi-vector
      // since it must be the RHS for a linear operator apply!
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
  
  return MEB::Derivative<Scalar>(DhDp_orig_mv,requiredOrientation);
  
}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::assembleParamDerivOutArgs(
  const ModelEvaluatorBase::OutArgs<Scalar> &wrappedOutArgs, // in
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs // in/out
  ) const
{

  typedef ModelEvaluatorBase MEB;
  typedef MEB::Derivative<Scalar> Deriv;
    
  Deriv DfDp;
  if ( !(DfDp=outArgs.get_DfDp(p_idx_)).isEmpty() ) {
    assembleParamDeriv( wrappedOutArgs.get_DfDp(p_idx_), DfDp );
  }

  const int Ng = outArgs.Ng();
  for ( int j = 0; j < Ng; ++j ) {
    Deriv DgDp;
    if ( !(DgDp=outArgs.get_DgDp(j,p_idx_)).isEmpty() ) {
      assembleParamDeriv( wrappedOutArgs.get_DgDp(j,p_idx_), DgDp );
    }
  }

}


template<class Scalar>
void DefaultLumpedParameterModelEvaluator<Scalar>::assembleParamDeriv(
  const ModelEvaluatorBase::Derivative<Scalar> &DhDp_orig, // in
  const ModelEvaluatorBase::Derivative<Scalar> &DhDp // in/out
  ) const
{

  typedef ModelEvaluatorBase MEB;

  const RCP<const MultiVectorBase<Scalar> >
    DhDp_orig_mv = DhDp_orig.getMultiVector();
  TEUCHOS_TEST_FOR_EXCEPTION(
    is_null(DhDp_orig_mv), std::logic_error,
    "Error, we currently can't handle non-multi-vector derivatives!" );

  const RCP<MultiVectorBase<Scalar> >
    DhDp_mv = DhDp.getMultiVector();
  TEUCHOS_TEST_FOR_EXCEPTION(
    is_null(DhDp_mv), std::logic_error,
    "Error, we currently can't handle non-multi-vector derivatives!" );

  switch( DhDp_orig.getMultiVectorOrientation() ) {
    case MEB::DERIV_MV_BY_COL:
      // DhDp = DhDp_orig * B
#ifdef TEUCHSO_DEBUG
      TEUCHOS_ASSERT(
        DhDp.getMultiVectorOrientation() == MEB::DERIV_MV_BY_COL );
#endif
      apply( *DhDp_orig_mv, NOTRANS, *B_, DhDp_mv.ptr() );
      // Above, we could generalize DhDp_oirg to just be a general linear
      // operator.
      break;
    case MEB::DERIV_TRANS_MV_BY_ROW:
      // (DhDp^T) = B^T * (DhDp_orig^T)  [DhDp_orig_mv is transposed!]
#ifdef TEUCHSO_DEBUG
      TEUCHOS_ASSERT(
        DhDp.getMultiVectorOrientation() == MEB::DERIV_TRANS_MV_BY_ROW );
#endif
      apply( *B_, CONJTRANS, *DhDp_orig_mv, DhDp_mv.ptr() );
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }

}


} // namespace Thyra


#endif // THYRA_DEFAUL_LUMPED_PARAMETER_LUMPED_MODEL_EVALUATOR_HPP
