//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorPFD_impl_hpp
#define Tempus_PhiEvaluatorPFD_impl_hpp

#include "Tempus_PhiEvaluatorPFD.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
PhiEvaluatorPFD<Scalar>::PhiEvaluatorPFD(std::string name)
  : PhiEvaluator<Scalar>(name),
    pfdMethod_("CN"),
    pfdMethods_()
{
  // Intialize PFD methods map
  pfdMethods_["IE"] = PFDMethod{
      {RealPF{1.0, -1.0}},
      {}};
  pfdMethods_["CN"] = PFDMethod{
      {RealPF{2.0, -4.0}},
      {-1.0}};
  pfdMethods_["EDIRK3"] = PFDMethod{
      {RealPF{3.0, -27.0 / 2.0}},
      {-7.0 / 2.0, -1.0 / 2.0}};

  // setup timers, if available
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  std::stringstream ss;
  ss << "Tempus::" << name;

  // set up timers for overall Phi-evaluation, Linop-vector applications, and solves
  const std::string phiLabel   = ss.str() + ": PhiEval";
  timerPhi_                    = Teuchos::TimeMonitor::getNewCounter(phiLabel);
  const std::string linOpLabel = ss.str() + ": LinOp";
  timerLinOp_                  = Teuchos::TimeMonitor::getNewCounter(linOpLabel);
  const std::string solveLabel = ss.str() + ": Solve";
  timerSolve_                  = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorPFD<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set(
      "PhiEvaluator Type", "PFD",
      "Method to approximate the phi-function evaluation.");

  pl->set(
      "PFD Method", "CN",
      "Partial fraction decomposition method:\n"
      "'IE': Implicit-Euler compatible,\n"
      "'CN': Crank-Nicolson compatible,\n"
      "'EDIRK3': EDIRK3 compatible.\n");

  return pl;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorPFD<Scalar>::computePhi(
    const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
    const int phi_order, const Scalar cdt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      phi_order < 0, std::logic_error,
      "Error - PhiEvaluator::computePhi() phi_order must be non-negative!\n");

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor phiTimer(*timerPhi_);
#endif

  auto pfdMethod = pfdMethods_.find(pfdMethod_);
  TEUCHOS_TEST_FOR_EXCEPTION(
      pfdMethod == pfdMethods_.end(), std::logic_error,
      "Error - pfdMethod not found!\n");

  typename PhiEvaluatorPFD<Scalar>::PFList pfs        = pfdMethod->second.pfs;
  typename PhiEvaluatorPFD<Scalar>::PFPolynomial poly = pfdMethod->second.poly;

  // check that model and PhiLinearSolver have been initialized properly
  this->checkInitialized();

  // compute phiv as a sum, start from zero
  Thyra::assign(phiv, Scalar(0.0));

  // allocate two temp vectors
  Teuchos::RCP<Thyra::VectorBase<Scalar>> Mtemp = Mrhs_b->clone_v();
  Teuchos::RCP<Thyra::VectorBase<Scalar>> temp  = Mrhs_b->clone_v();

  // polynomial part:
  // initialize the vector iteration
  Thyra::copy(*Mrhs_b, Mtemp.ptr());
  for (typename PhiEvaluatorPFD<Scalar>::PFPolynomial::size_type k = phi_order; k < poly.size(); k++) {
    const double poly_coeff = poly[k];
    if (std::cmp_greater(k, phi_order)) {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor linOpTimer(*timerLinOp_);
#endif
      // apply (-cdt*M*J) to temp, next step in the vector iteration
      this->applyJacobian(Mtemp.ptr(), temp);
      Thyra::scale(-cdt, Mtemp.ptr());
    }
    // invert the mass matrix and add the scaled result
    Thyra::assign(temp.ptr(), Scalar(0.0));
    this->solveMass(temp.ptr(), Mtemp);
    Thyra::Vp_StV(phiv, Scalar(poly_coeff), *temp);
  }

  // partial fraction part
  Thyra::SolveStatus<Scalar> sStatus;
  for (auto pf = pfs.begin(); pf != pfs.end(); ++pf) {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor solveTimer(*timerSolve_);
#endif

    const Scalar alpha = Scalar(pf->pole / cdt);
    const Scalar beta  = Scalar(1.0);

    Thyra::assign(temp.ptr(), Scalar(0.0));

    // we call solveMpJ without inArgs, it will use the existing Mass (lumped or not)
    // and Jacobian (without preconditioning).
    sStatus = this->phiLinSolv_->solveMpJ(temp.ptr(), Mrhs_b, alpha, beta);

    // this version will reassemble MpJ (never lumped), but has preconditoning
    // sStatus = this->phiLinSolv_->assembleAndsolveMpJ(this->inArgs_lin_, temp.ptr(), Mrhs_b, alpha, beta);

    // the factor for each pole is coeff / pole^{phi_order} / (-cdt)
    const double coeff_k = -pf->coeff / std::pow(pf->pole, phi_order) / cdt;
    Thyra::Vp_StV(phiv, Scalar(coeff_k), *temp);
  }

  // return the status of the final solve
  return sStatus;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorPFD<Scalar>::computePhis(
    const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
    const Scalar cdt,
    const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B)
{
  // This method gives the identical result to the slow reference implementation,
  // but at a fraction of the cost if Mrhs_B contains more than one non-null entry
  //
  // Thyra::assign(phiv, Scalar(0.0));
  // Teuchos::RCP<Thyra::VectorBase<Scalar>> temp = x->clone_v();
  // for (Thyra::Ordinal k = 0; k < Mrhs_B.size(); k++) {
  //   if (Mrhs_B[k] != Teuchos::null) {
  //     Thyra::assign(temp.ptr(), Scalar(0.0));
  //     sStatus = this->computePhi(temp.ptr(), k, cdt, Mrhs_B[k]);
  //     Thyra::Vp_V(phiv, *temp);
  //   }
  // }

  TEUCHOS_TEST_FOR_EXCEPTION(
      Mrhs_B.size() < 1, std::logic_error,
      "Error - PhiEvaluator::computePhi() Mrhs_B list must not be empty!\n");

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor phiTimer(*timerPhi_);
#endif

  auto pfdMethod = pfdMethods_.find(pfdMethod_);
  TEUCHOS_TEST_FOR_EXCEPTION(
      pfdMethod == pfdMethods_.end(), std::logic_error,
      "Error - pfdMethod not found!\n");

  typename PhiEvaluatorPFD<Scalar>::PFList pfs        = pfdMethod->second.pfs;
  typename PhiEvaluatorPFD<Scalar>::PFPolynomial poly = pfdMethod->second.poly;

  // check that model and PhiLinearSolver have been initialized properly
  this->checkInitialized();

  // compute phiv as a sum, start from zero
  Thyra::assign(phiv, Scalar(0.0));

  // allocate two zero initialized temp vectors
  Teuchos::RCP<Thyra::VectorBase<Scalar>> Mtemp = phiv->clone_v();
  Teuchos::RCP<Thyra::VectorBase<Scalar>> temp  = phiv->clone_v();

  // cannot underflow due to assertion above
  const Thyra::Ordinal max_phi_order = Mrhs_B.size() - 1;

  // polynomial part:
  Thyra::Ordinal first_nonnull_rhs = max_phi_order + 1;
  for (typename PhiEvaluatorPFD<Scalar>::PFPolynomial::size_type k = 0; k < poly.size(); k++) {
    if (std::cmp_less_equal(k, max_phi_order) && Mrhs_B[k] != Teuchos::null) {
      // if RHS is provided, insert it at the appropriate time into the iteration
      Thyra::Vp_V(Mtemp.ptr(), *Mrhs_B[k]);
      first_nonnull_rhs = std::min<Thyra::Ordinal>(first_nonnull_rhs, k);
    }
    if (std::cmp_greater_equal(k, first_nonnull_rhs)) {
      // invert the mass matrix and add the scaled result
      const double poly_coeff = poly[k];
      Thyra::assign(temp.ptr(), Scalar(0.0));
      this->solveMass(temp.ptr(), Mtemp);
      Thyra::Vp_StV(phiv, Scalar(poly_coeff), *temp);

      if (k + 1 < poly.size()) {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor linOpTimer(*timerLinOp_);
#endif
        // apply (-cdt*M*J) to temp, next step in the vector iteration
        this->applyJacobian(Mtemp.ptr(), temp);
        Thyra::scale(-cdt, Mtemp.ptr());
      }
    }
  }

  // partial fraction part
  Thyra::SolveStatus<Scalar> sStatus;
  for (auto pf = pfs.begin(); pf != pfs.end(); ++pf) {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor solveTimer(*timerSolve_);
#endif

    const Scalar alpha = Scalar(pf->pole / cdt);
    const Scalar beta  = Scalar(1.0);

    // assemble the combined RHS for partial fractions
    Thyra::assign(Mtemp.ptr(), Scalar(0.));
    for (Thyra::Ordinal k = 0; k <= max_phi_order; k++) {
      if (Mrhs_B[k] != Teuchos::null) {
        // the factor for each pole is coeff / pole^(phi_order) / (-cdt)
        const double coeff_k = -pf->coeff / std::pow(pf->pole, k) / cdt;
        Thyra::Vp_StV(Mtemp.ptr(), Scalar(coeff_k), *Mrhs_B[k]);
      }
    }

    Thyra::assign(temp.ptr(), Scalar(0.0));
    // we call solveMpJ without inArgs, it will use the existing Mass (lumped or not)
    // and Jacobian (without preconditioning).
    sStatus = this->phiLinSolv_->solveMpJ(temp.ptr(), Mtemp, alpha, beta);

    // this version will reassemble MpJ (never lumped), but has preconditoning
    // sStatus = this->phiLinSolv_->assembleAndsolveMpJ(this->inArgs_lin_, temp.ptr(), Mrhs_b, alpha, beta);

    Thyra::Vp_V(phiv, *temp);
  }

  // return the status of the final solve
  return sStatus;
}

template <class Scalar>
void PhiEvaluatorPFD<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  // pl->validateParametersAndSetDefaults(*getValidParameters());

  pfdMethod_ = pl->get<std::string>("PFD Method", pfdMethod_);
  if (pfdMethods_.find(pfdMethod_) == pfdMethods_.end()) {
    Teuchos::RCP<Teuchos::FancyOStream> l_out = this->getOStream();
    l_out->setOutputToRootOnly(0);
    *l_out << "PFD Method '" << pfdMethod_ << "'\n"
           << "is not yet implemented, continuing with 'CN'!\n";
    pfdMethod_ = "CN";
  }
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorPFD<Scalar>> createPhiEvaluatorPFD(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto phi = rcp(new PhiEvaluatorPFD<Scalar>());
  phi->setName("From createPhiEvaluatorPFD");

  if (pl != Teuchos::null)
    phi->setPhiEvaluatorValues(pl);

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorPFD_impl_hpp
