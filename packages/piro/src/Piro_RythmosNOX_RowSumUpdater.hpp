//@HEADER
//@HEADER

#ifndef PIRO_RYTHMOS_NOX_ROWSUM_UPDATER_OBSERVER_HPP
#define PIRO_RYTHMOS_NOX_ROWSUM_UPDATER_OBSERVER_HPP

#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_RCP.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Thyra_NonlinearSolver_NOX.hpp"

namespace Piro {


/** \brief For a Rythmos/NOX solve, this object updates the row sum scaling
 */
template<class Scalar>
class RythmosNOXRowSumUpdaterObserver : virtual public Rythmos::IntegrationObserverBase<Scalar>
{
public:

  RythmosNOXRowSumUpdaterObserver();  

  void set_f_scaling(const Teuchos::RCP< ::Thyra::VectorBase<Scalar> >& f_scaling);

  /** \name Overridden from IntegrationObserverBase */
  //@{

  Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > cloneIntegrationObserver() const;
  
  void 
  resetIntegrationObserver(const Rythmos::TimeRange<Scalar> &integrationTimeDomain);

  void observeStartTimeIntegration(const Rythmos::StepperBase<Scalar> &stepper);

  void observeEndTimeIntegration(const Rythmos::StepperBase<Scalar> &stepper);

  void observeStartTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  void observeCompletedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  void observeFailedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  //@}

private:

  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > inv_row_sum_vec_;

};


/** \brief Nonmember constructor.
 *
 * \relates RythmosNOXRowSumUpdaterObserver
 */
template<class Scalar>
Teuchos::RCP<RythmosNOXRowSumUpdaterObserver<Scalar> >
createRythmosNOXRowSumUpdaterObserver()
{
  const Teuchos::RCP<RythmosNOXRowSumUpdaterObserver<Scalar> > o = 
    Teuchos::rcp(new RythmosNOXRowSumUpdaterObserver<Scalar>);

  return o;
}


// //////////////////////////////////////////////////////
// Implementations

template<typename Scalar>
RythmosNOXRowSumUpdaterObserver<Scalar>::RythmosNOXRowSumUpdaterObserver()
{ 
}

template<typename Scalar>
void 
RythmosNOXRowSumUpdaterObserver<Scalar>::
set_f_scaling(const Teuchos::RCP< ::Thyra::VectorBase<Scalar> >& f_scaling)
{
  inv_row_sum_vec_ = f_scaling;
}

template<typename Scalar>
Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > 
RythmosNOXRowSumUpdaterObserver<Scalar>::cloneIntegrationObserver() const
{
  Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer = 
    Teuchos::rcp(new RythmosNOXRowSumUpdaterObserver<Scalar>(*this));
  return observer;
}
  
template<typename Scalar>
void 
RythmosNOXRowSumUpdaterObserver<Scalar>::
resetIntegrationObserver(const Rythmos::TimeRange<Scalar> &integrationTimeDomain)
{
}

template<typename Scalar>
void RythmosNOXRowSumUpdaterObserver<Scalar>::
observeStartTimeIntegration(const Rythmos::StepperBase<Scalar> &stepper)
{
}				

template<typename Scalar>
void RythmosNOXRowSumUpdaterObserver<Scalar>::
observeEndTimeIntegration(const Rythmos::StepperBase<Scalar> &stepper)
{
}				

template<typename Scalar>
void RythmosNOXRowSumUpdaterObserver<Scalar>::observeStartTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  //std::cout << "ROGER - computing a new RYTHMOS scaling" << std::endl;

  // First time through, the nox solver is not constructed (ctor on
  // demand) so don't use scaling on first time step.
  if (timeStepIter == 0)
    return;

  // Create rcp to do safe casting of object chain
  const RCP<const Rythmos::StepperBase<Scalar> > rcp_stepper =
    rcp(&stepper, false);

  const RCP<const Rythmos::SolverAcceptingStepperBase<Scalar> > stepper_with_solver =
    Teuchos::rcp_dynamic_cast<const Rythmos::SolverAcceptingStepperBase<Scalar> >(rcp_stepper);
  
  TEUCHOS_ASSERT(nonnull(stepper_with_solver));

  const RCP<const ::Thyra::NonlinearSolverBase<Scalar> > thyra_solver = 
    stepper_with_solver->getSolver();

  // NOTE:: NOX only supports double!!!
  RCP<const ::Thyra::NOXNonlinearSolver > nox_thyra_solver =
    Teuchos::rcp_dynamic_cast<const ::Thyra::NOXNonlinearSolver>(thyra_solver);

  TEUCHOS_ASSERT(nonnull(nox_thyra_solver));

  RCP<const NOX::Solver::Generic> solver = nox_thyra_solver->getNOXSolver();

  RCP<const NOX::Abstract::Group> group = solver->getSolutionGroupPtr();
  
  RCP<const NOX::Thyra::Group> thyra_group = 
    rcp_dynamic_cast<const NOX::Thyra::Group>(group);
  
  if (!thyra_group->isJacobian()) {
    RCP<NOX::Thyra::Group> tmp_nox_thyra_group = 
      Teuchos::rcp_const_cast<NOX::Thyra::Group>(thyra_group);
    TEUCHOS_ASSERT( !tmp_nox_thyra_group.is_null() );
    tmp_nox_thyra_group->computeJacobian();
  }
  
  RCP< const ::Thyra::LinearOpBase< double > > jac = 
    thyra_group->getJacobianOperator(); 	

  RCP< const ::Thyra::RowStatLinearOpBase< double > > row_stat_jac = 
    Teuchos::rcp_dynamic_cast< const ::Thyra::RowStatLinearOpBase< double > >(jac);

  TEUCHOS_ASSERT( !row_stat_jac.is_null() );

  if (inv_row_sum_vec_.is_null())
    inv_row_sum_vec_ = ::Thyra::createMember(jac->range());

  row_stat_jac->getRowStat( ::Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
			    inv_row_sum_vec_.ptr());

  //inv_row_sum_vec_->describe(*Teuchos::fancyOStream(Teuchos::RCP<std::ostream>(&std::cout, false)), Teuchos::VERB_EXTREME);

  //jac->describe(*Teuchos::fancyOStream(Teuchos::RCP<std::ostream>(&std::cout, false)), Teuchos::VERB_EXTREME);
}				

template<typename Scalar>
void RythmosNOXRowSumUpdaterObserver<Scalar>::observeCompletedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
}				

template<typename Scalar>
void RythmosNOXRowSumUpdaterObserver<Scalar>::observeFailedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
}				

}


#endif
