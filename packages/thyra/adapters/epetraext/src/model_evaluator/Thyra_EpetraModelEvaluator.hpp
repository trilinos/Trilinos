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

#ifndef THYRA_EPETRA_MODEL_EVALUATOR_HPP
#define THYRA_EPETRA_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Teuchos_Array.hpp"


namespace Thyra {


/** \brief Concrete Adapter subclass that takes an
 * <tt>EpetraExt::ModelEvaluator</tt> object and wraps it as a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * This class takes care of the basic details of wrapping and unwrapping
 * Epetra from Thyra objects.  This class is highly configurable and will be
 * maintained and modified in the future as the basic link between the Epetra
 * world and the Thyra world for nonlinear models and nonlinear algorithms.
 *
 * \section Thyra_EpetraModelEvaluator_Scaling_sec Scaling
 *
 * This class can handle scaling of the state function f(...) and of the state
 * variables x and all of the affected derivatives.
 *
 * The scaling for the state function can be set manually using
 * <tt>setStateFunctionScalingVec()</tt> or can be computed automatically
 * using the parameter <tt>"State Function Scaling"</tt> (see documentation
 * output from <tt>this->getValidParameters()->print(...)</tt>) in the input
 * parameter list set by <tt>setParameterList()</tt>.  Reguardless of how the
 * state function scaling is computed, it will compute a positive vector
 * <tt>s_f</tt> that defines a diagonal matrix <tt>S_f = diag(s_f)</tt> that
 * transforms the state function:
 
 \verbatim

    f(...) = S_f * f_orig(...)

 \endverbatim

 * where <tt>f_orig(...)</tt> is the original state function as computed by the
 * underlying <tt>EpetraExt::ModelEvaluator</tt> object and <tt>f(...)</tt> is
 * the state function as computed by <tt>evalModel()</tt>.
 *
 * The scaling for the state variables must be set manually using
 * <tt>Thyra::setStateVariableScalingVec()</tt>.  The vector that is set
 * <tt>s_x>/tt> defines a diagonal scaling matrix <tt>S_x = diag(s_x)</tt>
 * that transforms the variables as:
 
 \verbatim

    x = S_x * x_orig

 \endverbatim

 * where <tt>x_orig</tt> is the original unscaled state variable vector as
 * defined by the underlying <tt>EpetraExt::ModelEvaluator</tt> object and
 * <tt>x</tt> is the scaled state varaible vector as returned from
 * <tt>getNominalValues()</tt> and as accepted by <tt>evalModel()</tt>.  Note
 * that when the scaled variables <tt>x</tt> are passed into
 * <tt>evalModel</tt> that they are unscaled as:
 
 \verbatim

    x_orig = inv(S_x) * x

 \endverbatim

 * where <tt>inv(S_x)</tt> is the inverse of the diagonals of <tt>S_x</tt>
 * which is stored as a vector <tt>inv_s_x</tt>.
 *
 * Note how these scalings affect the state function:
 
 \verbatim

    f(x,...) = S_f * f_orig( inv(S_x)*x...)

 \endverbatim

 * which as the state/state Jacobian:
 
 \verbatim

    W = d(f)/d(x) = S_f * d(f_orig)/d(x_orig) * inv(S_x)

 \endverbatim

 * Currently, this class does not handle scalings of the parameters
 * <tt>p(l)</tt> or of the auxilary response functions <tt>g(j)(...)</tt>.
 *
 * The state varaible and state function scaling gives the following scaled
 * quantities:
 
 \verbatim

    f = S_f * f_orig

    W = S_f * W_orig * inv(S_x)

    DfDp(l) = S_f * DfDp_orig(l)

    g(j) = g_orig(j)

    DgDx_dot(j) = DgDx_dot_orig(j) * inv(S_x)

    DgDx(j) = DgDx_orig(j) * inv(S_x)
    
    DgDp(j,l) = DgDp_orig(j,l)

 \endverbatim

 * Since the scaling is done explicitly, the client never even sees the
 * orginal scaling and the linear solver (and contained preconditioner) are
 * computed from the scaled W shown above.
 *
 * ToDo: Describe how scaling affects the Hessian-vector products an how you just
 * need to scale the Lagrange mutipliers as:
 
 \verbatim

  u^T * f(...) = u^T * (S_f * f_orig(...)) = u_f^T * f_orig(...)

 \endverbatim

 * where <tt>u_f = S_f * u</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup EpetraExt_Thyra_Op_Vec_adapters_grp
 */
class EpetraModelEvaluator
  : public ModelEvaluatorDefaultBase<double>,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  EpetraModelEvaluator();

  /** \brief . */
  EpetraModelEvaluator(
    const RCP<const EpetraExt::ModelEvaluator> &epetraModel,
    const RCP<LinearOpWithSolveFactoryBase<double> > &W_factory
    );

  /** \brief . */
  void initialize(
    const RCP<const EpetraExt::ModelEvaluator> &epetraModel,
    const RCP<LinearOpWithSolveFactoryBase<double> > &W_factory
    );

  /** \brief . */
  RCP<const EpetraExt::ModelEvaluator> getEpetraModel() const;

  /** \brief Set the nominal values.
   *
   * Warning, if scaling is being used, these must be according to the scaled
   * values, not the original unscaled values.
   */
  void setNominalValues( const ModelEvaluatorBase::InArgs<double>& nominalValues );
  
  /** \brief Set the state variable scaling vector <tt>s_x</tt> (see above).
   *
   * This function must be called after <tt>intialize()</tt> or the
   * constructur in order to set the scaling vector correctly!
   *
   * ToDo: Move this into an external strategy class object!
   */
  void setStateVariableScalingVec(
    const RCP<const Epetra_Vector> &stateVariableScalingVec
    );
  
  /** \brief Get the state variable scaling vector <tt>s_x</tt> (see
   * above). */
  RCP<const Epetra_Vector>
  getStateVariableInvScalingVec() const;
  
  /** \brief Get the inverse state variable scaling vector <tt>inv_s_x</tt>
   * (see above). */
  RCP<const Epetra_Vector>
  getStateVariableScalingVec() const;
  
  /** \brief Set the state function scaling vector <tt>s_f</tt> (see
   * above). */
  void setStateFunctionScalingVec(
    const RCP<const Epetra_Vector> &stateFunctionScalingVec
    );
  
  /** \brief Get the state function scaling vector <tt>s_f</tt> (see
   * above). */
  RCP<const Epetra_Vector>
  getStateFunctionScalingVec() const;

  /** \brief . */
  void uninitialize(
    RCP<const EpetraExt::ModelEvaluator> *epetraModel = NULL,
    RCP<LinearOpWithSolveFactoryBase<double> > *W_factory = NULL
    );
  
  /** \brief . */
  const ModelEvaluatorBase::InArgs<double>& getFinalPoint() const;

  /** \brief . */
  bool finalPointWasSolved() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
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
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  RCP<const VectorSpaceBase<double> > get_x_space() const;
  /** \brief . */
  RCP<const VectorSpaceBase<double> > get_f_space() const;
  /** \brief . */
  RCP<const VectorSpaceBase<double> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  RCP<const VectorSpaceBase<double> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getUpperBounds() const;
  /** \brief . */
  RCP<LinearOpBase<double> > create_W_op() const;
  /** \brief Returns null currently. */
  RCP<PreconditionerBase<double> > create_W_prec() const;
  /** \breif . */
  RCP<const LinearOpWithSolveFactoryBase<double> > get_W_factory() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> createInArgs() const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<double>      &finalPoint
    ,const bool                                   wasSolved
    );

  //@}

  // Made public to simplify implementation but this is harmless to be public.
  // Clients should not deal with this type.
  enum EStateFunctionScaling { STATE_FUNC_SCALING_NONE, STATE_FUNC_SCALING_ROW_SUM };

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  RCP<LinearOpBase<double> > create_DfDp_op_impl(int l) const;
  /** \brief . */
  RCP<LinearOpBase<double> > create_DgDx_dot_op_impl(int j) const;
  /** \brief . */
  RCP<LinearOpBase<double> > create_DgDx_op_impl(int j) const;
  /** \brief . */
  RCP<LinearOpBase<double> > create_DgDp_op_impl(int j, int l) const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<double> &inArgs,
    const ModelEvaluatorBase::OutArgs<double> &outArgs
    ) const;

  //@}

private:

  // ////////////////////
  // Private types

  typedef Teuchos::Array<RCP<const Epetra_Map> > p_map_t;
  typedef Teuchos::Array<RCP<const Epetra_Map> > g_map_t;
  typedef std::vector<bool> p_map_is_local_t;
  typedef std::vector<bool> g_map_is_local_t;

  typedef Teuchos::Array<RCP<const VectorSpaceBase<double> > >
  p_space_t;
  typedef Teuchos::Array<RCP<const VectorSpaceBase<double> > >
  g_space_t;

  // /////////////////////
  // Private data members

  RCP<const EpetraExt::ModelEvaluator> epetraModel_;

  RCP<Teuchos::ParameterList> paramList_;

  RCP<LinearOpWithSolveFactoryBase<double> > W_factory_;

  RCP<const Epetra_Map> x_map_;
  p_map_t p_map_;
  g_map_t g_map_;
  p_map_is_local_t p_map_is_local_;
  p_map_is_local_t g_map_is_local_;
  RCP<const Epetra_Map> f_map_;

  RCP<const VectorSpaceBase<double> > x_space_;
  p_space_t p_space_;
  RCP<const VectorSpaceBase<double> > f_space_;
  g_space_t g_space_;

  mutable ModelEvaluatorBase::InArgs<double> nominalValues_;
  mutable ModelEvaluatorBase::InArgs<double> lowerBounds_;
  mutable ModelEvaluatorBase::InArgs<double> upperBounds_;
  mutable bool nominalValuesAndBoundsAreUpdated_;

  ModelEvaluatorBase::InArgs<double> finalPoint_;

  EStateFunctionScaling stateFunctionScaling_;
  mutable RCP<const Epetra_Vector> stateFunctionScalingVec_;

  RCP<const Epetra_Vector> stateVariableScalingVec_; // S_x
  mutable RCP<const Epetra_Vector> invStateVariableScalingVec_; // inv(S_x)
  mutable EpetraExt::ModelEvaluator::InArgs epetraInArgsScaling_;
  mutable EpetraExt::ModelEvaluator::OutArgs epetraOutArgsScaling_;
  
  mutable RCP<Epetra_Vector> x_unscaled_;
  mutable RCP<Epetra_Vector> x_dot_unscaled_;

  mutable ModelEvaluatorBase::InArgs<double> prototypeInArgs_;
  mutable ModelEvaluatorBase::OutArgs<double> prototypeOutArgs_;
  mutable bool currentInArgsOutArgs_;

  bool finalPointWasSolved_;

  // //////////////////////////
  // Private member functions

  /** \brief . */
  void convertInArgsFromEpetraToThyra(
    const EpetraExt::ModelEvaluator::InArgs &epetraInArgs,
    ModelEvaluatorBase::InArgs<double> *inArgs
    ) const;

  /** \brief . */
  void convertInArgsFromThyraToEpetra(
    const ModelEvaluatorBase::InArgs<double> &inArgs,
    EpetraExt::ModelEvaluator::InArgs *epetraInArgs
    ) const;

  /** \brief . */
  void convertOutArgsFromThyraToEpetra(
    // Thyra form of the outArgs
    const ModelEvaluatorBase::OutArgs<double> &outArgs,
    // Epetra form of the unscaled output arguments 
    EpetraExt::ModelEvaluator::OutArgs *epetraUnscaledOutArgs,
    // The passed-in form of W
    RCP<LinearOpBase<double> > *W_op,
    RCP<EpetraLinearOp> *efwdW,
    // The actual Epetra object passed to the underylying EpetraExt::ModelEvaluator
    RCP<Epetra_Operator> *eW
    ) const;

  /** \brief . */
  void preEvalScalingSetup(
    EpetraExt::ModelEvaluator::InArgs *epetraInArgs,
    EpetraExt::ModelEvaluator::OutArgs *epetraUnscaledOutArgs,
    const RCP<Teuchos::FancyOStream> &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  /** \brief . */
  void postEvalScalingSetup(
    const EpetraExt::ModelEvaluator::OutArgs &epetraUnscaledOutArgs,
    const RCP<Teuchos::FancyOStream> &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  /** \brief . */
  void finishConvertingOutArgsFromEpetraToThyra(
    const EpetraExt::ModelEvaluator::OutArgs &epetraOutArgs,
    RCP<LinearOpBase<double> > &W_op,
    RCP<EpetraLinearOp> &efwdW,
    RCP<Epetra_Operator> &eW,
    const ModelEvaluatorBase::OutArgs<double> &outArgs // Output!
    ) const;
  // 2007/08/03: rabartl: Above, I pass many of the RCP objects by non-const
  // reference since I don't want the compiler to perform any implicit
  // conversions on this RCP objects.

  /** \brief . */
  void updateNominalValuesAndBounds() const;

  /** \brief . */
  void updateInArgsOutArgs() const;

  /** \brief . */
  RCP<EpetraLinearOp> create_epetra_W_op() const;
  
};


//
// Utility functions
//


/** \brief .
 * \relates EpetraModelEvaluator
 */
RCP<EpetraModelEvaluator>
epetraModelEvaluator(
  const RCP<const EpetraExt::ModelEvaluator> &epetraModel,
  const RCP<LinearOpWithSolveFactoryBase<double> > &W_factory
  );


/** \brief .
 * \relates EpetraModelEvaluator
 */
ModelEvaluatorBase::EDerivativeMultiVectorOrientation
convert( const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation &mvOrientation );


/** \brief .
 * \relates EpetraModelEvaluator
 */
EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation
convert( const ModelEvaluatorBase::EDerivativeMultiVectorOrientation &mvOrientation );


/** \brief .
 * \relates EpetraModelEvaluator
 */
ModelEvaluatorBase::DerivativeProperties
convert( const EpetraExt::ModelEvaluator::DerivativeProperties &derivativeProperties );


/** \brief .
 * \relates EpetraModelEvaluator
 */
ModelEvaluatorBase::DerivativeSupport
convert( const EpetraExt::ModelEvaluator::DerivativeSupport &derivativeSupport );


/** \brief .
 * \relates EpetraModelEvaluator
 */
EpetraExt::ModelEvaluator::Derivative
convert(
  const ModelEvaluatorBase::Derivative<double> &derivative,
  const RCP<const Epetra_Map> &fnc_map,
  const RCP<const Epetra_Map> &var_map
  );


} // namespace Thyra


#endif // THYRA_EPETRA_MODEL_EVALUATOR_HPP
