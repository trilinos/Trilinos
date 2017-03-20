// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef TEMPUS_TEST_SINCOSXDOTDOT_MODEL_DECL_HPP
#define TEMPUS_TEST_SINCOSXDOTDOT_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief This problem is analogous  to the SinCosModel problem, except written in its usual xdotdot form. 
 * The governing ODE is:  
  *   \f[
  *   \ddot{x}=-x
  *   \f]
  * The initial conditions are:
  *   \f{eqnarray*}{
  *     x(0) & = & b*sin(\phi)\\
  *     \dot{x}(0) & = & b*cos(\phi) 
  *   \f}
  * for some constants \f$\phi\f$ and \f$b\f$.  It can be shown that the exact solution to this problem is:
  *    \f[
  *    x(t) = b*sin(t+\phi) 
  *    \f]
  * We consider the problem for \f$t\in [0,1]\f$, and we take \f$b = 1\f$, \f$\phi = 0\f$.  
  */

template<class Scalar>
class SinCosXdotdotModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  SinCosXdotdotModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

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

private:
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_vec_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_vec_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_dot_vec_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > p_init_;
  int vecLength_; //Number of state unknowns (1)
  int numResponses_; //Number of responses (1)
  double b_; //parameter defining problem
  double phi_; //parameter defining problem
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  mutable bool isInitialized_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;

};



} // namespace Tempus_Test
#endif // TEMPUS_TEST_SINCOSXDOTDOT_MODEL_DECL_HPP
