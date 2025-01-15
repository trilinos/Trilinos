// @HEADER
// @HEADER

#ifndef USER_APP_TRANSIENT_REDUCED_OBJECTIVE_HPP
#define USER_APP_TRANSIENT_REDUCED_OBJECTIVE_HPP

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientForwardSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientAdjointSensitivity.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"

#include "user_app_Utilities.hpp"

namespace ROL {

template <typename Real>
class TransientReducedObjective : public virtual ROL::Objective<Real> {
public:

  TransientReducedObjective(
    const Teuchos::RCP<Teuchos::ParameterList>& input_params,
    const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
    const Teuchos::RCP<Teuchos::ParameterList>& objective_params,
    const Teuchos::RCP<std::ostream>& os);

  virtual ~TransientReducedObjective() {}

  //! Compute value of objective
  Real value( const ROL::Vector<Real> &x, Real &tol );

  //! Compute gradient of objective
  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol );

  //! Set response target for computing objective
  void set_target(const Teuchos::RCP<Thyra::VectorBase<Real> >& target);
  void set_target(const Teuchos::RCP<ROL::Vector<Real> >& target);

  //! Helper function to create optimization vector
  Teuchos::RCP<ROL::Vector<Real> > create_design_vector() const;

  //! Helper function to create a response vector
  Teuchos::RCP<ROL::Vector<Real> > create_response_vector() const;

  //! Helper function to run tempus, computing responses and derivatives
  void run_tempus(ROL::Vector<Real>& r, const ROL::Vector<Real>& p);
  void run_tempus(const Thyra::ModelEvaluatorBase::InArgs<Real>&  inArgs,
                  const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs);

private:

  // Teuchos::RCP<Thyra::ModelEvaluator<Real> > model_;
  Teuchos::RCP<Teuchos::ParameterList> tempus_params_;
  Teuchos::RCP<Thyra::VectorBase<Real> > target_;
  std::string objective_type_;
  std::string sensitivity_method_;
  int param_index_;
  int response_index_;
  bool use_fd_gradient_;

  // Objects neeeded to rebuild the time integrator since it seems to
  // not be able to reset itself
  Teuchos::RCP<Teuchos::ParameterList> input_params_;
  Teuchos::RCP<const Teuchos::Comm<int>> comm_;
  Teuchos::RCP<Thyra::ModelEvaluator<Real>> model_;
  Teuchos::RCP<panzer::GlobalData> global_data_;
  Teuchos::RCP<panzer_stk::STK_Interface> mesh_;
  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> response_library_;
  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> stk_io_response_library_;
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>> lin_obj_factory_;
  Teuchos::RCP<panzer::GlobalIndexer> global_indexer_;
  Teuchos::RCP<std::ostream> os_;
  bool print_debug_;
};

}

#endif
