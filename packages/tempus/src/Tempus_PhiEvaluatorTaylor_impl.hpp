//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_impl_hpp
#define Tempus_PhiEvaluatorTaylor_impl_hpp

#include "Tempus_PhiEvaluatorTaylor.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_RCPDecl.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorBase.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorTaylor<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set(
      "PhiEvaluator Type", "Taylor",
      "Method to approximate the phi-function evaluation.");

  pl->set<int>(
      "Expansion Order", 10,
      "Taylor degree N in sum_{j=0}^N L^j v/(j+p)!.\n"
      "\n"
      "The default is 10.");

  return pl;
}

template <class Scalar>
std::string PhiEvaluatorTaylor<Scalar>::description() const
{
  return ("Tempus::PhiEvaluatorTaylor - '" + this->name_ + "'");
}

template <class Scalar>
void PhiEvaluatorTaylor<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  PhiEvaluator<Scalar>::describe(out, verbLevel);

  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  if ((Teuchos::as<int>(verbLevel) ==
       Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))) {
    *l_out << "  Expansion Order  = " << getExpansionOrder() << std::endl;
  }

  *l_out << std::string(this->description().length() + 8, '-') << std::endl;
}


template <class Scalar>
Thyra::SolveStatus<Scalar>
PhiEvaluatorTaylor<Scalar>::computeLinOpPhi(const int phi_order,
                                            const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
                                            const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v, const Scalar cdt)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      phi_order < 0,
      std::invalid_argument,
      "LinOpPhi: phi_order must be nonnegative.");

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor phitimer(*timerPhi_);
#endif

  const int expansionOrder = getExpansionOrder();

  // phi_k(L) * v is in range(L)
  const auto rangeSpace = L->range();

  // build 1. / (phi_order!)
  int k;
  Scalar inv_factorial = Scalar(1.);
  for (k = 1; k <= phi_order; ++k)
  {
    inv_factorial /= Scalar(k);
  }

  // Iteration vector d_0 = v / (phi_order!)
  Teuchos::RCP<Thyra::VectorBase<Scalar>> d_k = Thyra::createMember(rangeSpace);

  if (phi_order > 0)
  {
    Thyra::V_StV(d_k.ptr(), inv_factorial, *v);
    // v := d_0
    Thyra::assign(v, *d_k);
  }
  else
  {
    Thyra::assign(d_k.ptr(), *v);
  }

  // allocate temporary vector
  Teuchos::RCP<Thyra::VectorBase<Scalar>> next = Thyra::createMember(rangeSpace);

  Scalar norm_d_k = Thyra::norm_inf(*d_k);
  Scalar overflow = Thyra::norm_inf(*v);

  Thyra::SolveStatus<Scalar> sStatus;
  sStatus.achievedTol = norm_d_k;
  sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;

  // Iteratively compute d_k = (L^(k-phi_order) d_{k-1}) / (k!) and add to result
  for (k = phi_order + 1; k <= expansionOrder + phi_order; ++k)
  {
    {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor linoptimer(*timerLinOp_);
#endif
      // next <- K * d_k
      // TODO: do we need the temp vector?
      Thyra::apply(*L, Thyra::NOTRANS, *d_k, next.ptr());

      // multiply the update by 1/k and store in d_k
      Thyra::V_StV(d_k.ptr(), Scalar(1.) / Scalar(k), *next);
    }
    // add d_k to the final result
    Thyra::Vp_V(v, *d_k);

    norm_d_k = Thyra::norm_inf(*d_k);

    // overflow is an upper bound on Thyra::norm_inf(*v);
    // it tracks how large the update matExpTemp could have gotten in intermediate iterations
    // any number larger than overflow * machine_eps should not be affected much by roundoff
    overflow += norm_d_k;

    // TODO: refine this and make dependent on Scalar type
    const Scalar cutoff = 1e22;
    if (overflow > cutoff)
    {
      sStatus.achievedTol = norm_d_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_UNCONVERGED;
      break;
    }

    // terminate if the update drops below likely significance
    if (norm_d_k < overflow / cutoff)
    {
      sStatus.achievedTol = norm_d_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
      break;
    }
  }

  std::stringstream ss;
  ss << "Taylor: Norm of solution=" << Thyra::norm_inf(*v)
     << " overflow=" << overflow
     << " final update=" << norm_d_k
     << " achieved in it. " << k << ".";
  sStatus.message = ss.str();

  // std::cout << sStatus.message << std::endl;

  return sStatus;
}

template <class Scalar>
void PhiEvaluatorTaylor<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  setExpansionOrder(pl->get<int>("Expansion Order", 10));

  // TODO: make this configurable?
  this->useAtildeForSingleRHS_ = false;
}


// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> phi = Teuchos::rcp(new PhiEvaluatorTaylor<Scalar>());
  phi->setName("From createPhiEvaluatorTaylor");

  if (pl != Teuchos::null)
    phi->setPhiEvaluatorValues(pl);

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorTaylor_impl_hpp
