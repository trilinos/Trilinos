// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SINCOSMODEL_H
#define SINCOSMODEL_H

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"


/** \brief Sine-Cosine model problem from Tempus.
  * This is a canonical Sine-Cosine differential equation
  *   \f[
  *   \mathbf{\ddot{x}}=-\mathbf{x}
  *   \f]
  * with a few enhancements. We start with the exact solution to the
  * differential equation
  *   \f{eqnarray*}{
  *     x_{0}(t) & = & a+b*\sin((f/L)*t+\phi)\\
  *     x_{1}(t) & = & b*(f/L)*\cos((f/L)*t+\phi)
  *   \f}
  * then the form of the model is
  * \f{eqnarray*}{
  *   \frac{d}{dt}x_{0}(t) & = & x_{1}(t)\\
  *   \frac{d}{dt}x_{1}(t) & = & \left(\frac{f}{L}\right)^{2}(a-x_{0}(t))
  * \f}
  * where the default parameter values are \f$a=0\f$, \f$f=1\f$, and \f$L=1\f$,
  * and the initial conditions
  * \f{eqnarray*}{
  *   x_{0}(t_{0}=0) & = & \gamma_{0}[=0]\\
  *   x_{1}(t_{0}=0) & = & \gamma_{1}[=1]
  * \f}
  * determine the remaining coefficients
  * \f{eqnarray*}{
  *   \phi & = & \arctan(((f/L)/\gamma_{1})*(\gamma_{0}-a))-(f/L)*t_{0}[=0]\\
  *   b & = & \gamma_{1}/((f/L)*cos((f/L)*t_{0}+\phi))[=1]
  * \f}

  * Therefore this model has three model parameters and two initial conditions
  * which effect the exact solution as above.
  * \f[
  *   \mathbf{p}=(a,f,L)
  * \f]
  * \f[
  *   \dot{\mathbf{x}}=\mathbf{F}(\mathbf{x},t,\mathbf{p})
  * \f]
  * where
  * \f{eqnarray*}{
  *   F_{0} & = & x_{1}\\
  *   F_{1} & = & \left(\frac{f}{L}\right)^{2}(a-x_{0})
  * \f}

  * The exact sensitivities, \f$\mathbf{s}=\partial\mathbf{x}/\partial\mathbf{p}\f$,
  * for the problem are specified as
  * \f[
  *   \mathbf{s}(t)=\left[\begin{array}{cc}
                                                   *   1 & 0\\
  *   \left(\frac{b}{L}\right)t\,\cos\left(\left(\frac{f}{L}\right)t+\phi\right) & \left(\frac{b}{L}\right)\cos\left(\left(\frac{f}{L}\right)t+\phi\right)-\frac{b\, f\, t}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)\\
  *   -\frac{b\, f\, t}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right) & -\frac{b\, f}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{2}\, t}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)
  *   \end{array}\right]
  * \f]
  * and for the default initial conditions, \f$\phi=0\f$ and \f$b=1\f$
  * \f[
  *   \mathbf{s}(t=0)=\left[\begin{array}{cc}
  *   1 & 0\\
  *   0 & \frac{b}{L}\\
  *   0 & -\frac{f}{L^{2}}
  *   \end{array}\right]
  * \f]
  * The time differentiated sensitivities, \f$\dot{\mathbf{s}}=\partial\mathbf{s}/\partial t=\partial/\partial t(\partial\mathbf{x}/\partial\mathbf{p})=\partial/\partial\mathbf{p}(\partial\mathbf{x}/\partial t)\f$
  * are
  * \f[
  *   \dot{\mathbf{s}}(t)=\left[\begin{array}{cc}
  *   0 & 0\\
  *   \left(\frac{b}{L}\right)\cos\left(\left(\frac{f}{L}\right)t+\phi\right)-\frac{b\, f\, t}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right) & -\frac{2b\, f}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)\left(\frac{b}{L}\right)-\frac{b\, f^{2}\, t}{L^{3}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)\\
 *   -\frac{b\, f}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{2}\, t}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right) & \frac{2b\, f^{2}}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{3}\, t}{L^{4}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)
  *   \end{array}\right]
  * \f]
  */


using LO = Tpetra::Map<>::local_ordinal_type;
using GO = Tpetra::Map<>::global_ordinal_type;
using Scalar = double; 
typedef Tpetra::Map<LO,GO>  Tpetra_Map;
typedef Tpetra::Vector<Scalar,LO,GO>  Tpetra_Vector;
typedef Tpetra::MultiVector<Scalar,LO,GO>  Tpetra_MultiVector;
typedef Tpetra::Operator<Scalar,LO,GO>  Tpetra_Operator;
typedef Tpetra::CrsGraph<LO,GO>  Tpetra_CrsGraph;
typedef Tpetra::CrsMatrix<Scalar,LO,GO>  Tpetra_CrsMatrix;
typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LO, GO> ConverterT;



class SinCosModel
    : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
      public Teuchos::ParameterListAcceptorDefaultBase

{
  public:

  // Constructor
  SinCosModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  // Exact sensitivity solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSensSolution(int j, double t) const;

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  //@}

  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

private:

  void setupInOutArgs_() const;

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
    ) const;
  //@}

  void calculateCoeffFromIC_();

protected:
  int dim_;         ///< Number of state unknowns (2)
  int Np_;          ///< Number of parameter vectors (1)
  int np_;          ///< Number of parameters in this vector (2)
  int Ng_;          ///< Number of observation functions (1)
  int ng_;          ///< Number of elements in this observation function (1)
  bool haveIC_;     ///< false => no nominal values are provided (default=true)
  bool acceptModelParams_; ///< Changes inArgs to require parameters
  bool useDfDpAsTangent_; ///< Treat DfDp OutArg as tangent (df/dx*dx/dp+df/dp)
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > DxDp_space_;

  // Parameters for the model:  x_0(t) = a + b*sin(f*t+phi)
  //                            x_1(t) = b*f*cos(f*t+phi)
  double a_;     ///< Model parameter
  double f_;     ///< Model parameter
  double L_;     ///< Model parameter
  double phi_;   ///< Parameter determined from the IC
  double b_;     ///< Parameter determined from the IC
  double t0_ic_; ///< Time value where the initial condition is specified
  double x0_ic_; ///< Initial condition for x0
  double x1_ic_; ///< Initial condition for x1
};


/// Non-member constructor
//Teuchos::RCP<SinCosModel> sineCosineModel(
//  Teuchos::RCP<Teuchos::ParameterList> pList_)
//{
//  Teuchos::RCP<SinCosModel> model = rcp(new SinCosModel(pList_));
//  return(model);
//}

//! Adjoint for SinCosModel
/*!
 * This model evaluator modifies evalModel() to compute the adjoint operator
 * instead of the forward operator.
 */
class SinCosModelAdjoint
  : public SinCosModel
{
  public:

  // Constructor
  SinCosModelAdjoint(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null) : SinCosModel(pList) {}

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
    ) const;
  //@}
};


#endif // SINCOSMODEL
