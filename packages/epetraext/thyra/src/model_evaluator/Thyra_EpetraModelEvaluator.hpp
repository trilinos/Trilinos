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

#ifndef THYRA_EPETRA_MODEL_EVALUATOR_HPP
#define THYRA_EPETRA_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"


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

    f(...) = S_f * f_hat(...)

 \endverbatim

 * where <tt>f_hat(...)</tt> is the original state function as computed by the
 * underlying <tt>EpetraExt::ModelEvaluator</tt> object and <tt>f(...)</tt> is
 * the state function as computed by <tt>evalModel()</tt>.
 *
 * The scaling for the state variables must be set manually using
 * <tt>Thyra::setStateVariableScalingVec()</tt>.  The vector that is set
 * <tt>s_x>/tt> defines a diagonal scaling matrix <tt>S_x = diag(s_x)</tt>
 * that transforms the variables as:
 
 \verbatim

    x = S_x * x_hat

 \endverbatim

 * where <tt>x_hat</tt> is the original unscaled state variable vector as
 * defined by the underlying <tt>EpetraExt::ModelEvaluator</tt> object and
 * <tt>x</tt> is the scaled state varaible vector as returned from
 * <tt>getNominalValues()</tt> and as accepted by <tt>evalModel()</tt>.  Note
 * that when the scaled variables <tt>x</tt> are passed into
 * <tt>evalModel</tt> that they are unscaled as:
 
 \verbatim

    x_hat = inv(S_x) * x

 \endverbatim

 * where <tt>inv(S_x)</tt> is the inverse of the diagonals of <tt>S_x</tt>
 * which is stored as a vector <tt>inv_s_x</tt>.
 *
 * Note how these scalings affect the state function:
 
 \verbatim

    f(x,...) = S_f * f_hat( inv(S_x)*x...)

 \endverbatim

 * which as the state/state Jacobian:
 
 \verbatim

    W = d(f)/d(x) = S_f * d(f_hat)/d(x_hat) * inv(S_x)

 \endverbatim

 * Currently, this class does not handle scalings of the parameters
 * <tt>p(l)</tt> or of the auxilary response functions <tt>g(j)(...)</tt>.
 *
 * The state varaible and state function scaling gives the following scaled
 * quantities:
 
 \verbatim

    f = S_f * f_hat

    W = S_f * W_hat * inv(S_x)

    DfDp(l) = S_f * DfDp_hat(l)

    g(j) = g_hat(j)

    DgDx(j) = DgDx_hat(j) * inv(S_x)
    
    DgDp(j,l) = DgDp_hat(j,l)

 \endverbatim

 * Since the scaling is done explicitly, the client never even sees the
 * orginal scaling and the linear solver (and contained preconditioner) are
 * computed from the scaled W shown above.
 *
 * ToDo: Describe how scaling affects the Hessian-vector products an how you just
 * need to scale the Lagrange mutipliers as:
 
 \verbatim

  u^T * f(...) = u^T * (S_f * f_hat(...)) = u_f^T * f_hat(...)

 \endverbatim

 * where <tt>u_f = S_f * u</tt>.
 *
 * ToDo: Finish documentation!
 */
class EpetraModelEvaluator
  : public ModelEvaluator<double>,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  EpetraModelEvaluator();

  /** \brief . */
  EpetraModelEvaluator(
    const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> &epetraModel,
    const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > &W_factory
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> &epetraModel,
    const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > &W_factory
    );

  /** \brief . */
  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> getEpetraModel() const;

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
    const Teuchos::RefCountPtr<const Epetra_Vector> &stateVariableScalingVec
    );
  
  /** \brief Get the state variable scaling vector <tt>s_x</tt> (see above). */
  Teuchos::RefCountPtr<const Epetra_Vector>
  getStateVariableScalingVec() const;
  
  /** \brief Set the state function scaling vector <tt>s_f</tt> (see above). */
  void setStateFunctionScalingVec(
    const Teuchos::RefCountPtr<const Epetra_Vector> &stateFunctionScalingVec
    );
  
  /** \brief Get the state function scaling vector <tt>s_f</tt> (see above). */
  Teuchos::RefCountPtr<const Epetra_Vector>
  getStateFunctionScalingVec() const;

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> *epetraModel = NULL,
    Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > *W_factory = NULL
    );
  
  /** \brief . */
  const ModelEvaluatorBase::InArgs<double>& getFinalPoint() const;

  /** \brief . */
  bool finalPointWasSolved() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_x_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_f_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getUpperBounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > create_W() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_W_op() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_DfDp_op(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_DgDx_op(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_DgDp_op( int j, int l ) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> createInArgs() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<double> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<double>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<double>  &outArgs
    ) const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<double>      &finalPoint
    ,const bool                                   wasSolved
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  // Made public to simplify implementation but this is harmless to be public.
  // Clients should not deal with this type.
  enum EStateFunctionScaling { STATE_FUNC_SCALING_NONE, STATE_FUNC_SCALING_ROW_SUM };

private:

  // ////////////////////
  // Private types

  typedef std::vector<Teuchos::RefCountPtr<const Epetra_Map> > p_map_t;
  typedef std::vector<Teuchos::RefCountPtr<const Epetra_Map> > g_map_t;
  typedef std::vector<bool> p_map_is_local_t;
  typedef std::vector<bool> g_map_is_local_t;

  typedef std::vector<Teuchos::RefCountPtr<const SpmdVectorSpaceDefaultBase<double> > >
  p_space_t;
  typedef std::vector<Teuchos::RefCountPtr<const SpmdVectorSpaceDefaultBase<double> > >
  g_space_t;

  // /////////////////////
  // Private data members

  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> epetraModel_;

  Teuchos::RefCountPtr<Teuchos::ParameterList> paramList_;

  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > W_factory_;

  Teuchos::RefCountPtr<const Epetra_Map> x_map_;
  p_map_t p_map_;
  g_map_t g_map_;
  p_map_is_local_t p_map_is_local_;
  p_map_is_local_t g_map_is_local_;
  Teuchos::RefCountPtr<const Epetra_Map> f_map_;

  Teuchos::RefCountPtr<const SpmdVectorSpaceDefaultBase<double> > x_space_;
  p_space_t p_space_;
  Teuchos::RefCountPtr<const SpmdVectorSpaceDefaultBase<double> > f_space_;
  g_space_t g_space_;

  mutable ModelEvaluatorBase::InArgs<double> nominalValues_;
  mutable ModelEvaluatorBase::InArgs<double> lowerBounds_;
  mutable ModelEvaluatorBase::InArgs<double> upperBounds_;
  mutable bool nominalValuesAndBoundsAreUpdated_;

  ModelEvaluatorBase::InArgs<double> finalPoint_;

  EStateFunctionScaling stateFunctionScaling_;
  mutable Teuchos::RefCountPtr<const Epetra_Vector> stateFunctionScalingVec_;

  Teuchos::RefCountPtr<const Epetra_Vector> stateVariableScalingVec_; // S_x
  mutable Teuchos::RefCountPtr<const Epetra_Vector> invStateVariableScalingVec_; // inv(S_x)
  mutable EpetraExt::ModelEvaluator::InArgs epetraInArgsScaling_;
  mutable EpetraExt::ModelEvaluator::OutArgs epetraOutArgsScaling_;
  
  mutable Teuchos::RefCountPtr<Epetra_Vector> x_unscaled_;
  mutable Teuchos::RefCountPtr<Epetra_Vector> x_dot_unscaled_;

  bool finalPointWasSolved_;

  // //////////////////////////
  // Private member functions

  void convertInArgsFromEpetraToThyra(
    const EpetraExt::ModelEvaluator::InArgs &epetraInArgs,
    ModelEvaluatorBase::InArgs<double> *inArgs
    ) const;

  void convertInArgsFromThyraToEpetra(
    const ModelEvaluatorBase::InArgs<double> &inArgs,
    EpetraExt::ModelEvaluator::InArgs *epetraInArgs
    ) const;

  void convertOutArgsFromThyraToEpetra(
    const ModelEvaluatorBase::OutArgs<double> &outArgs,
    EpetraExt::ModelEvaluator::OutArgs *epetraOutArgs
    ) const;

  void convertOutArgsFromEpetraToThyra(
    const EpetraExt::ModelEvaluator::OutArgs &epetraOutArgs,
    ModelEvaluatorBase::OutArgs<double> *outArgs
    ) const;

  void updateNominalValuesAndBounds() const;
  
};


//
// Utility functions
//


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
  const Teuchos::RefCountPtr<const Epetra_Map> &fnc_map,
  const Teuchos::RefCountPtr<const Epetra_Map> &var_map
  );


} // namespace Thyra


#endif // THYRA_EPETRA_MODEL_EVALUATOR_HPP
