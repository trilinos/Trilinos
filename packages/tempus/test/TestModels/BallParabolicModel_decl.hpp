// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef TEMPUS_TEST_BALLPARABOLIC_MODEL_DECL_HPP
#define TEMPUS_TEST_BALLPARABOLIC_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Tempus_Test {

/** \brief This is the "parabolic ball" model problem from Piro.
  * This is a canonical differential equation model of a ball thrown up 
  * in the air, and taking on a parabolic trajectory:
  *   \f[
  *   \ddot{x}=-1
  *   \f]
  * The initial conditions are:
  *   \f{eqnarray*}{
  *     x(0) & = & 0\\
  *     \dot{x}(0) & = & 1
  *   \f}
  * It can be shown that the exact solution to this problem is:
  *    \f[ 
  *    x(t) = t(1-0.5t)
  *    \f]
  * We consider the problem for \f$t\in [0,2]\f$ .  An EpetraExt version of this test is implemented in 
  * Piro::MockModelEval_B (see Trilinos/packages/piro/test), where it is used to test the Piro (EpetraExt)
  * Newmark-Beta scheme (see input_Solver_NB.xml input file).  
  */

template<class Scalar>
class BallParabolicModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  BallParabolicModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

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
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
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
  int numParameters_; //Number of parameters (1) 
  mutable ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  mutable bool isInitialized_;
  
};



} // namespace Tempus_Test
#endif // TEMPUS_TEST_BALLPARABOLIC_MODEL_DECL_HPP
