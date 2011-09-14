//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRA_EXT_MODEL_EVALUATOR_SCALING_TOOLS_H
#define EPETRA_EXT_MODEL_EVALUATOR_SCALING_TOOLS_H


#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_Utils.hpp"


namespace EpetraExt {


/** \defgroup EpetraExt_ModelEvaluator_ScalingTools_grp Scaling Tools for EpetraExt::ModelEvaluator.
 *
 * These scaling functions implement scaling of input variables and output
 * functions and their derivatives.
 *
 * The scaling vectors are stored in
 * <tt>EpetraExt::ModelEvaluator::InArgs</tt> and
 * <tt>EpetraExt::ModelEvaluator::OutArgs</tt> objects in order to enhance
 * maintainability and to avoid programming errors.  This will result in some
 * wasted space but it should not be excessive if used carefully.
 *
 * First, consider scaling of the state function.  Reguardless of how the
 * state function scaling is computed, it will be represented as a positive
 * vector <tt>s_f</tt> that defines a diagonal matrix <tt>S_f = diag(s_f)</tt>
 * that transforms the state function:
 
 \verbatim

    f(...) = S_f * f_hat(...)

 \endverbatim

 * where <tt>f_hat(...)</tt> is the original unscaled state function as
 * computed by the underlying <tt>EpetraExt::ModelEvaluator</tt> object and
 * <tt>f(...)</tt> is the scaled state function.
 *
 * Next, consider the scaling of the state varaibles.  The scaling for the
 * state variables is defined by a positive vector <tt>s_x>/tt> defines a
 * diagonal scaling matrix <tt>S_x = diag(s_x)</tt> that transforms the
 * variables as:
 
 \verbatim

    x = S_x * x_hat

 \endverbatim

 * where <tt>x_hat</tt> is the original unscaled state variable vector as
 * defined by the underlying <tt>EpetraExt::ModelEvaluator</tt> object and
 * <tt>x</tt> is the scaled state varaible vector.  Note that when the scaled
 * variables <tt>x</tt> are passed into <tt>evalModel</tt> that they must be
 * unscaled as:
 
 \verbatim

    x_hat = inv(S_x) * x

 \endverbatim

 * where <tt>inv(S_x)</tt> is the inverse of the diagonals of <tt>S_x</tt>
 * which is stored as a positive vector <tt>inv_s_x</tt>.  Since unscaling the
 * variables as shown above is more common than scaling the original
 * variables, the scaling vector will be stored as <tt>inv_s_x</tt> and not as
 * <tt>s_x</tt>.
 *
 * Note how these scalings affect the state function:
 
 \verbatim

    f( x_dot, x, ... ) = S_f * f_hat( inv(S_x)*x_dot, inv(S_x)*x, ... )

 \endverbatim

 * which has the state/state Jacobian:
 
 \verbatim

    W = alpha * d(f)/d(x_dot) + beta * d(f)/d(x)
      = S_f * ( alpha * d(f_hat)/d(x_hat) + beta * d(f_hat)/d(x) ) * inv(S_x)

 \endverbatim

 * Currently, these functions do not handle scalings of the parameters
 * <tt>p(l)</tt> or of the auxilary response functions <tt>g(j)(...)</tt>.
 *
 * The state varaible and state function scaling gives the following scaled
 * quantities:
 
 \verbatim

    f = S_f * f_hat

    W = S_f * W_hat * inv(S_x)

    DfDp(l) = S_f * DfDp_hat(l), for l=0...Np-1

    g(j) = g_hat(j), for j=0...Ng-1

    DgDx_dot(j) = DgDx_dot_hat(j) * inv(S_x), for j=0...Ng-1

    DgDx(j) = DgDx_hat(j) * inv(S_x), for j=0...Ng-1
    
    DgDp(j,l) = DgDp_hat(j,l), for j=0...Ng-1, l=0...Np-1

 \endverbatim

 * ToDo: Describe how scaling of the state function <tt>S_f</tt> affects the
 * Hessian-vector products an how you just need to scale the Lagrange
 * mutipliers as:
 
 \verbatim

  u^T * f(...) = u^T * (S_f * f_hat(...)) = u_f^T * f_hat(...)

 \endverbatim

 * where <tt>u_f = S_f * u</tt>.
 *
 * ToDo: Also describe how scaling of the state varaibles <tt>S_x</tt> affects
 * Hessian-vector products and other related quantities.
 *
 * \section EpetraExt_ModelEvaluator_ScalingTools_Maintenance_sec Maintenance of these tools
 *
 * These scaling tools must be updated whenever the <tt>InArgs</tt> or
 * <tt>OutArgs</tt> classes are augmented.  However, not every use case with
 * the model evaluator requires scaling so scaling with respect to some inputs
 * and some outputs may never be needed and therefore never need to be seen by
 * these tools.  However, there is some danger in ignoring inputs and outputs
 * in these scaling tools since some objects may be silently unscaled and
 * could cause hard to track down bugs.
 *
 * ToDo: Finish documentation!
 *
 */
//@{


/** \brief Gather the nominal values from a model evaluator.
 *
 * ToDo: Finish documentation!
 *
 * ToDo: Perhaps refactor the EpetraExt::ModelEvaluator interface to return
 * nominal values as a single InArgs object?
 */
void gatherModelNominalValues(
  const ModelEvaluator &model,
  ModelEvaluator::InArgs *nominalValues
  );


/** \brief Gather the lower and upper bounds from a model evaluator.
 *
 * ToDo: Finish documentation!
 *
 * ToDo: Perhaps refactor the EpetraExt::ModelEvaluator interface to return
 * lower and upper bounds as two different InArgs objects?
 */
void gatherModelBounds(
  const ModelEvaluator &model,
  ModelEvaluator::InArgs *lowerBounds,
  ModelEvaluator::InArgs *upperBounds
  );


/** \brief Scale the original unscaled variables into the scaled variables.
 *
 * \param origVars
 *          [in] The orginal unscaled variables.  The input data pointed to in
 *          this object will not be changed by this function call.
 *
 * \param varScalings
 *          [in] The variable scaling vectors.  These scaling vectors must be
 *          stored as the inverse scaling vector, such as <tt>inv_s_x</tt> as
 *          described in \ref EpetraExt_ModelEvaluator_ScalingTools_grp.
 *
 * \param scaledVars
 *          [in/out] On output, contains copies of the scaled varaibles.  On
 *          first call, <tt>*scaledVars</tt> may be empty.  Any storage that
 *          does not exist will be created.  On subsequent calls the storage
 *          will be reused.  Warning! const casting will be used to allow the
 *          modification of the vector objects pointed to by
 *          <tt>*scaledVars</tt> so don't bank on those vectors not being
 *          modified.  Any vectors pointed to by <tt>*scaledVars</tt> is fair
 *          game to be modified in this function. <it>Precondition:</it>
 *          <tt>scaledVars!=0</tt>.
 *
 * \param out
 *          [out] If <tt>out != 0</tt> then output will be sent to
 *          <tt>*out</tt>.
 *
 * \param verbLevel
 *          [in] Determines the verbosity level for output sent to <tt>*out</tt>.
 *
 * The assumpition, of course, is that the InArgs objects <tt>origVars</tt>,
 * <tt>varScalings</tt>, and <tt>*scaledVars</tt> will all have been created
 * by the same <tt>EpetraExt::ModelEvaluator::createOutArgs()</tt> function
 * call.
 */
void scaleModelVars(
  const ModelEvaluator::InArgs &origVars,
  const ModelEvaluator::InArgs &varScalings,
  ModelEvaluator::InArgs *scaledVars,
  Teuchos::FancyOStream *out = 0,
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_LOW
  );

/** \brief Scale the lower and upper model variable bounds.
 *
 * ToDo: Finish documentation!
 */
void scaleModelBounds(
  const ModelEvaluator::InArgs &origLowerBounds,
  const ModelEvaluator::InArgs &origUpperBounds,
  const double infBnd,
  const ModelEvaluator::InArgs &varScalings,
  ModelEvaluator::InArgs *scaledLowerBounds,
  ModelEvaluator::InArgs *scaledUpperBounds,
  Teuchos::FancyOStream *out = 0,
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_LOW
  );


/** \brief Unscale the scaled variables.
 *
 * \param scaledVars
 *          [in] The scaled varaibles.  The input data pointed to in this
 *          object will not be modified by this function call.
 *
 * \param varScalings
 *          [in] The variable scaling vectors.  These scaling vectors must be
 *          stored as the inverse scaling vector, such as <tt>inv_s_x</tt> as
 *          described in \ref EpetraExt_ModelEvaluator_ScalingTools_grp.
 *
 * \param origVars
 *          [in/out] On output, contains copies of the unscaled varaibles.  On
 *          first call, <tt>*origVars</tt> may be empty.  Any storage that
 *          does not exist will be created.  On subsequent calls the storage
 *          will be reused.  Warning! const casting will be used to allow the
 *          modification of the vector objects pointed to by
 *          <tt>*origVars</tt> so don't bank on those vectors not being
 *          modified.  Any vectors pointed to by <tt>*scaledVars</tt> is fair
 *          game to be modified in this function. <it>Precondition:</it>
 *          <tt>origVars!=0</tt>.
 *
 * \param out
 *          [out] If <tt>out != 0</tt> then output will be sent to
 *          <tt>*out</tt>.
 *
 * \param verbLevel
 *          [in] Determines the verbosity level for output sent to <tt>*out</tt>.
 * 
 */
void unscaleModelVars(
  const ModelEvaluator::InArgs &scaledVars,
  const ModelEvaluator::InArgs &varScalings,
  ModelEvaluator::InArgs *origVars,
  Teuchos::FancyOStream *out = 0,
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_LOW
  );


/** \brief Scale the output functions and their derivative objects.
 *
 * \param origFuncs
 *          [in/out] On input, contains the unscaled functions and derivative
 *          objects.  On output, many to most of the objects pointed to by
 *          this object will be scaled in place.  See details below.
 *
 * \param varScalings
 *          [in] The variable scaling vectors.  These scaling vectors must be
 *          stored as the inverse scaling vector, such as <tt>inv_s_x</tt> as
 *          described in \ref EpetraExt_ModelEvaluator_ScalingTools_grp.
 *
 * \param funcScalings
 *          [in] The function scaling vectors.  These scaling vectors must be
 *          stored as the forward scaling vector, such as <tt>s_f</tt> as
 *          described in \ref EpetraExt_ModelEvaluator_ScalingTools_grp.
 *
 * \param scaledFuncs
 *          [out] On output, points to in-place scaled functions and
 *          derivatives.  No new storage is created in this object.  Any
 *          functions or derivative objects in <tt>origFuncs</tt> that could
 *          not be scaled will not be presented in this object.  An output
 *          object may not be scaled if the object does not allow scaling.
 *          For example, if a derivative object is defined only as an
 *          <tt>Epetra_Operator</tt> object and can not be dynamic cased to a
 *          <tt>Epetra_RowMatrix</tt>, then no explicit scaling will be
 *          possible.  <it>Precondition:</it> <tt>scaledFuncs!=0</tt>.
 *
 * \param allFuncsWhereScaled
 *          [out] On output, determines if all of the functions and
 *          derivatives in <tt>origFuncs</tt> where successfully scaled.  If
 *          <tt>*allFuncsWhereScaled==true</tt> on output, then all of the
 *          functions in <tt>origFuncs</tt> where scaled and are represented
 *          in <tt>*scaledFuncs</tt>.  If <tt>*allFuncsWhereScaled==false</tt>
 *          on output, then at least one of the functions or derivative
 *          objects present in <tt>origFuncs</tt> was not scaled and is not
 *          present in <tt>*scaledFuncs</tt>.  It is up to the client to
 *          search <tt>*scaledFuncs</tt> and compare to <tt>origFuncs</tt> to
 *          see what is missing; Sorry :-(.
 *
 * \param out
 *          [out] If <tt>out != 0</tt> then output will be sent to
 *          <tt>*out</tt>.
 *
 * \param verbLevel
 *          [in] Determines the verbosity level for output sent to <tt>*out</tt>.
 *
 * In general, any output objects that are <tt>Epetra_MultiVector</tt> (or
 * <tt>Epetra_Vector</tt>) or dynamic castable to <tt>Epetra_RowMatrix</tt>
 * can be explicitly scaled by this function.  Objects that are simply
 * <tt>Epetra_Operator</tt> objects can not and will not be scaled by this
 * function and the client is on thier own.
 *
 * ToDo: Consider using some composite Epetra_Operator classes to create
 * implicitly scaled Epetra_Operator objects and put in an option for doing
 * this or not.
 */
void scaleModelFuncs(
  const ModelEvaluator::OutArgs &origFuncs,
  const ModelEvaluator::InArgs &varScalings,
  const ModelEvaluator::OutArgs &funcScalings,
  ModelEvaluator::OutArgs *scaledFuncs,
  bool *allFuncsWhereScaled,
  Teuchos::FancyOStream *out = 0,
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_LOW
  );


/** \brief Create an inverse scaling vector.
 *
 * This function may actually result in <tt>scalingVector</tt> being saved for
 * use later embedded within <tt>returnValue</tt>.
 */
Teuchos::RefCountPtr<const Epetra_Vector>
createInverseModelScalingVector(
  Teuchos::RefCountPtr<const Epetra_Vector> const& scalingVector
  );


/** \brief Scale a vector given its inverse scaling vector.
 *
 * \param origVars
 *          [in] The vector of original unscaled varaibles.
 *
 * \param invVarScaling
 *          [in] The inverse scaling vector.
 *
 * \param scaledVars
 *          [out] On output, will contain the scaled varaibles:
 *          <tt>scaledVars[i] = origVars[i] / invScaleVector[i]</tt>, for
 *          <tt>i=0...n-1</tt>.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>scaledVars != 0</tt>
 * <li><tt>origVars.Map().SameAs(invVarScaling.Map()) == true</tt>
 * <li><tt>origVars.Map().SameAs(scaledVars->Map()) == true</tt>
 * </ul>
 *
 * This function is used by the <tt>scaleModelVars()</tt> function to scale
 * each of the varaible vectors.
 */
void scaleModelVarsGivenInverseScaling(
  const Epetra_Vector &origVars,
  const Epetra_Vector &invVarScaling,
  Epetra_Vector *scaledVars
  );


/** \brief Scale the model variable bounds. */
void scaleModelVarBoundsGivenInverseScaling(
  const Epetra_Vector &origLowerBounds,
  const Epetra_Vector &origUpperBounds,
  const double infBnd,
  const Epetra_Vector &invVarScaling,
  Epetra_Vector *scaledLowerBounds,
  Epetra_Vector *scaledUpperBounds
  );


/** \brief Unscale a vector given its inverse scaling vector.
 *
 * \param scaledVars
 *          [in] The vector of original unscaled varaibles.
 *
 * \param invVarScaling
 *          [in] The inverse scaling vector.
 *
 * \param origVars
 *          [out] On output, will contain the unscaled varaibles:
 *          <tt>origVars[i] = invScaleVector[i] * scaledVars[i]</tt>, for
 *          <tt>i=0...n-1</tt>.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>origVars != 0</tt>
 * <li><tt>scaledVars.Map().SameAs(invVarScaling.Map()) == true</tt>
 * <li><tt>scaledVars.Map().SameAs(origVars->Map()) == true</tt>
 * </ul>
 *
 * This function is used by the function <tt>unscaleModelVars()</tt> function
 * to unscale each of the varaible vectors.
 */
void unscaleModelVarsGivenInverseScaling(
  const Epetra_Vector &origVars,
  const Epetra_Vector &invVarScaling,
  Epetra_Vector *scaledVars
  );


/** \brief Scale (in place) an output function vector given its forward
 * scaling vector.
 *
 * \param fwdFuncScaling
 *          [in] The forward scaling vector.
 *
 * \param funcs
 *          [in/out] On input, contains the vector of unscaled functions.  On
 *          output, contains the scaled functions: <tt>scaledFuncs[i] *=
 *          fwdFuncScaling[i]</tt>.
 *
 * <b>Preconditions:</b><ul>
 * <li> ???
 * </ul>
 *
 * This function is used by the <tt>scaleModelFuncs()</tt> function to scale
 * each of the otuput function vectors.
 */
void scaleModelFuncGivenForwardScaling(
  const Epetra_Vector &fwdFuncScaling,
  Epetra_Vector *funcs
  );


/** \brief Scale (in place) an output first-order function derivative object
 * represented as an Epetra_Operator given its function and variable scaling.
 *
 * \param invVarScaling
 *          [in] If <tt>invVarScaling !=0</tt>, then this represents the
 *          inverse varaible scaling (e.g. <tt>inv_s_x</tt>).  If
 *          <tt>invVarScaling==0</tt>, the identity scaling is assumed.
 *
 * \param fwdFuncScaling
 *          [in] If <tt>fwdFuncScaling !=0</tt>, then this represents the
 *          forward function scaling (e.g. <tt>s_f</tt>).  If
 *          <tt>fwdFuncScaling==0</tt>, the identity scaling is assumed.
 *
 * \param funcDerivOp
 *          [in/out] If scaling could be performed, then on output, this object
 *          will be scaled.  Otherwise, it will not be scaled (see <tt>didScaling</tt>)..
 *
 * \param didScaling
 *          [out] On output <tt>*didScaling==true</tt> if the scaling
 *          could be performed.
 *
 * <b>Preconditions:</b><ul>
 * <li> ???
 * </ul>
 *
 * This function is used by the <tt>scaleModelFuncs()</tt> function to scale
 * each of the otuput function first derivative objects.
 */
void scaleModelFuncFirstDerivOp(
  const Epetra_Vector *invVarScaling,
  const Epetra_Vector *fwdFuncScaling,
  Epetra_Operator *funcDerivOp,
  bool *didScaling
  );


/** \brief Scale (in place) an output first-order function derivative object
 * given its function and variable scaling.
 *
 * \param origFuncDeriv
 *          [in/out] The vector of original unscaled function derivative.  If
 *          this object can be scaled, then on output it will be scaled in
 *          place and <tt>scaledFuncDeriv</tt> will also point to the scaled
 *          derivative object.
 *
 * \param invVarScaling
 *          [in] If <tt>invVarScaling !=0</tt>, then this represents the
 *          inverse varaible scaling (e.g. <tt>inv_s_x</tt>).  If
 *          <tt>invVarScaling==0</tt>, the identity scaling is assumed.
 *
 * \param fwdFuncScaling
 *          [in] If <tt>fwdFuncScaling !=0</tt>, then this represents the
 *          forward function scaling (e.g. <tt>s_f</tt>).  If
 *          <tt>fwdFuncScaling==0</tt>, the identity scaling is assumed.
 *
 * \param scaledFuncDeriv
 *          [out] If scaling could be performed, then on output, this object
 *          will point to the scaled function derivative.  If not, then
 *          <tt>scaledFuncDeriv.isEmpty() == true</tt> on output.
 *
 * \param didScaling
 *          [out] On output <tt>*didScaling==true</tt> if the scaling
 *          could be performed.
 *
 * <b>Preconditions:</b><ul>
 * <li> ???
 * </ul>
 *
 * This function is used by the <tt>scaleModelFuncs()</tt> function to scale
 * each of the otuput function first derivative objects.
 */
void scaleModelFuncFirstDeriv(
  const ModelEvaluator::Derivative &origFuncDeriv,
  const Epetra_Vector *invVarScaling,
  const Epetra_Vector *fwdFuncScaling,
  ModelEvaluator::Derivative *scaledFuncDeriv,
  bool *didScaling
  );


/** \brief Class that gets and sets x_dot in an InArgs object. */
class InArgsGetterSetter_x_dot {
public:

  std::string getName() const { return "x_dot"; }

  Teuchos::RefCountPtr<const Epetra_Vector>
  getVector( const ModelEvaluator::InArgs &inArgs ) const
  {
    return inArgs.get_x_dot();
  }

  void setVector(
    const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot,
    ModelEvaluator::InArgs *inArgs
    ) const
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(!inArgs);
#endif
    inArgs->set_x_dot(x_dot);
  }

};


/** \brief Class that gets and sets x in an InArgs object. */
class InArgsGetterSetter_x {
public:

  std::string getName() const { return "x"; }

  Teuchos::RefCountPtr<const Epetra_Vector>
  getVector( const ModelEvaluator::InArgs &inArgs ) const
  {
    return inArgs.get_x();
  }

  void setVector(
    const Teuchos::RefCountPtr<const Epetra_Vector> &x,
    ModelEvaluator::InArgs *inArgs
    ) const
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(!inArgs);
#endif
    inArgs->set_x(x);
  }

};


/** \brief Class that gets and sets p(l) in an InArgs object. */
class InArgsGetterSetter_p {
public:

  InArgsGetterSetter_p( int l ) : l_(l) {}

  std::string getName() const
  { return "p["+Teuchos::Utils::toString(l_)+"]"; }

  Teuchos::RefCountPtr<const Epetra_Vector>
  getVector( const ModelEvaluator::InArgs &inArgs ) const
  {
    return inArgs.get_p(l_);
  }

  void setVector(
    const Teuchos::RefCountPtr<const Epetra_Vector> &p_l,
    ModelEvaluator::InArgs *inArgs
    ) const
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(!inArgs);
#endif
    inArgs->set_p(l_,p_l);
  }

private:

  int l_;

  InArgsGetterSetter_p(); // Not defined!

};


/** \brief Class that gets and sets f in an OutArgs object. */
class OutArgsGetterSetter_f {
public:

  Teuchos::RefCountPtr<Epetra_Vector>
  getVector( const ModelEvaluator::OutArgs &outArgs ) const
  {
    return outArgs.get_f();
  }

  void setVector(
    const Teuchos::RefCountPtr<Epetra_Vector> &f,
    ModelEvaluator::OutArgs *outArgs
    ) const
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(!outArgs);
#endif
    outArgs->set_f(f);
  }

};


/** \brief Class that gets and sets g(j) in an OutArgs object. */
class OutArgsGetterSetter_g {
public:

  OutArgsGetterSetter_g( int j ) : j_(j) {}

  Teuchos::RefCountPtr<Epetra_Vector>
  getVector( const ModelEvaluator::OutArgs &outArgs ) const
  {
    return outArgs.get_g(j_);
  }

  void setVector(
    const Teuchos::RefCountPtr<Epetra_Vector> &g_j,
    ModelEvaluator::OutArgs *outArgs
    ) const
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(!outArgs);
#endif
    outArgs->set_g(j_,g_j);
  }

private:

  int j_;

  OutArgsGetterSetter_g(); // Not defined!

};


//@}


} // namespace EpetraExt


#endif // EPETRA_EXT_MODEL_EVALUATOR_SCALING_TOOLS_H
