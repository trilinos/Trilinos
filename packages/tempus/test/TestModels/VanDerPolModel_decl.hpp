// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef TEMPUS_TEST_VANDERPOL_MODEL_DECL_HPP
#define TEMPUS_TEST_VANDERPOL_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief van der Pol model problem for nonlinear electrical circuit.
  * This is a canonical equation of a nonlinear oscillator (Hairer, Norsett,
  * and Wanner, pp. 111-115, and Hairer and Wanner, pp. 4-5) for an electrical
  * circuit.  The scaled form of this problem can be written as
  * \f{eqnarray*}{
  *   f_0 & = & \dot{x}_0(t) - x_1(t) \\
  *   f_1 & = & \dot{x}_1(t) - [(1-x_0^2)x_1-x_0]/\epsilon
  * \f}
  * where \f$\epsilon = 10^{-6}\f$ and the initial conditions are
  * \f{eqnarray*}{
  *   x_0(t_0=0) & = & 2 \\
  *   x_1(t_0=0) & = & 0
  * \f}
  * and the initial time derivatives are
  * \f{eqnarray*}{
  *   \dot{x}_0(t_0=0) & = & x_1(t_0=0) = 0 \\
  *   \dot{x}_1(t_0=0) & = & [(1-x_0^2)x_1-x_0]/\epsilon = -2/\epsilon
  * \f}
  * Hairer and Wanner suggest the output times of \f$t = 1,2,3,4,...,11\f$
  * For \f$\epsilon = 0\f$, the solution becomes
  * \f{eqnarray*}{
  *   \ln \left|x_0\right| - \frac{x_0^2}{2} & = & t + C \\
  *   x_1 & = & \frac{x_0}{1-x_0^2}
  * \f}
  * where \f$C =\ln \left|x_0(t=0)\right| - \frac{x_0^2(t=0)}{2} =-1.306853.\f$
  */

template<class Scalar>
class VanDerPolModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  VanDerPolModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

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

  int dim_;         ///< Number of state unknowns (2)
  int Np_;          ///< Number of parameter vectors (1)
  int np_;          ///< Number of parameters in this vector (1)
  int Ng_;          ///< Number of observation functions (0)
  int ng_;          ///< Number of elements in this observation function (0)
  bool haveIC_;     ///< false => no nominal values are provided (default=true)
  bool acceptModelParams_; ///< Changes inArgs to require parameters
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  // Parameters for the model:
  Scalar epsilon_; ///< This is a model parameter
  Scalar t0_ic_;   ///< initial time
  Scalar x0_ic_;   ///< initial condition for x0
  Scalar x1_ic_;   ///< initial condition for x1
};


/// Non-member constructor
//Teuchos::RCP<VanDerPolModel> sineCosineModel(
//  Teuchos::RCP<Teuchos::ParameterList> pList_)
//{
//  Teuchos::RCP<VanDerPolModel> model = rcp(new VanDerPolModel(pList_));
//  return(model);
//}


} // namespace Tempus_Test
#endif // TEMPUS_TEST_VANDERPOL_MODEL_DECL_HPP
