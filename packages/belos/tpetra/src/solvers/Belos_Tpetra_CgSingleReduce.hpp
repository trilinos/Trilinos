// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_CG_SINGLE_REDUCE_HPP
#define BELOS_TPETRA_CG_SINGLE_REDUCE_HPP

#include "Belos_Tpetra_Krylov.hpp"
#include "Tpetra_idot.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::MultiVector<>::scalar_type,
         class MV = Tpetra::MultiVector<SC>,
         class OP = Tpetra::Operator<SC>>
class CgSingleReduce: public Krylov<SC, MV, OP> {
private:
  using base_type = Krylov<SC, MV, OP>;

public:
  CgSingleReduce () = default;

  CgSingleReduce (const Teuchos::RCP<const OP>& A) :
    base_type::Krylov (A)
  {}

  virtual ~CgSingleReduce () = default;

protected:
  using vec_type = typename base_type::vec_type;
  
  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in X0/out X
               vec_type& B, // in B/out R
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input) override
  {
    using std::endl;
    using dev_type = typename MV::device_type;
    using ATS = Kokkos::ArithTraits<SC>;
    using magnitude_type = typename ATS::mag_type;
    using ATM = Kokkos::ArithTraits<magnitude_type>;
    using dot_type = typename MV::dot_type;
    
    const SC ONE = ATS::one ();
    SolverOutput<SC> output {};

    // compute initial residual
    vec_type MR (B.getMap ());
    magnitude_type beta_old = ATM::zero ();
    if (input.precoSide == "none") {
      beta_old = ATS::real (B.dot (B));
    }
    else {
      M.apply (B, MR);
      beta_old = ATS::real (B.dot (MR));
    }
    magnitude_type r_norm = std::sqrt (beta_old);
    magnitude_type r_norm_orig = r_norm;

    // quick return
    magnitude_type metric = this->getConvergenceMetric (r_norm, r_norm_orig, input);
    if (metric <= input.tol) {
      output.absResid = r_norm;
      output.relResid = ATM::one ();
      output.numIters = 0;
      output.converged = true;
      // R doesn't exist yet, so we don't have to copy R back to B
      // here, as we do below.
      return output;
    }

    // local vectors
    // matrix containing R and AR
    MV R_AR (B.getMap (), 2);
    vec_type  R =  * (R_AR.getVectorNonConst (0));
    vec_type AR =  * (R_AR.getVectorNonConst (1));
    Tpetra::deep_copy (R, B);
    // results of [R AR]'*R
    magnitude_type RAR;
    Kokkos::View<dot_type*, dev_type> RR_RAR ("results[numVecs]", 2);
    auto RR_RAR_host = Kokkos::create_mirror_view (RR_RAR);
    vec_type P (R, Teuchos::Copy);
    vec_type AP (P.getMap ());

    // Initial step
    // AR = A*R
    magnitude_type PAP;
    if (input.precoSide == "none") {
      A.apply (R, AR);
      PAP = ATS::real (R.dot (AR));
    }
    else {
      M.apply (R, MR);
      A.apply (MR, AR);

      // [beta_old, PAP] = [MR,AR]'*MR
      // TODO: idot is used for now.
      //beta_old = MR.dot(  R );
      //PAP      = MR.dot( AR );
      auto req = Tpetra::idot (RR_RAR, R_AR, MR);
      req->wait ();

      Kokkos::deep_copy (RR_RAR_host, RR_RAR);
      beta_old = ATS::real (RR_RAR_host(0));
      PAP = ATS::real (RR_RAR_host(1));

      r_norm = std::sqrt (beta_old);
    }
    magnitude_type alpha      = beta_old / PAP;
    magnitude_type beta       = ATM::zero ();
    magnitude_type beta_new   = ATM::zero ();
    // main loop
    for (int iter = 0; iter < input.maxNumIters; ++iter) {
      if (outPtr != nullptr) {
        *outPtr << "Iteration " << (iter+1) << " of " << input.maxNumIters
                << ": r_norm: " << r_norm;
      }

      // * search direction *
      // P = R + beta*P
      if (input.precoSide == "none") {
        P.update (ONE, R, static_cast<SC> (beta));
      } else {
        P.update (ONE, MR, static_cast<SC> (beta));
      }
      // AP = AR + beta*AP
      AP.update (ONE, AR, static_cast<SC> (beta));

      // * solution update *
      // X = X + alpha*P
      X.update (static_cast<SC> (alpha), P, ONE);
      // R = R - alpha*AP
      R.update (static_cast<SC> (-alpha), AP, ONE);

      // * matrix op *
      // AR = A*R
      if (input.precoSide == "none") {
        A.apply (R, AR);
        // [RR,RAR] = [R,AR]'*R
        // TODO: idot is used for now.
        auto req = Tpetra::idot (RR_RAR, R_AR, R);
        req->wait ();
      }
      else {
        M.apply (R, MR);
        A.apply (MR, AR);
        // [RR,RAR] = [R,AR]'*MR
        // TODO: idot is used for now.
        // TODO: need to compute R'*R for convergence check.
        auto req = Tpetra::idot (RR_RAR, R_AR, MR);
        req->wait ();
      }
      // * all-reduce *
      Kokkos::deep_copy (RR_RAR_host, RR_RAR);
      beta_new = ATS::real (RR_RAR_host(0));
      RAR = ATS::real (RR_RAR_host(1));

      // convergence check
      r_norm = std::sqrt( beta_new );
      metric = this->getConvergenceMetric (r_norm, r_norm_orig, input);
      if (outPtr != nullptr) {
        *outPtr << ", r_norm: " << r_norm << ", RAR: " << RAR << ", metric: " << metric;
      }
      if (metric <= input.tol) {
        if (outPtr != nullptr) {
          *outPtr << endl;
        }
        output.absResid = r_norm;
        output.relResid = r_norm / r_norm_orig;
        output.numIters = iter + 1;
        output.converged = true;

        Tpetra::deep_copy (B, R);
        return output;
      }
      else if (iter + 1 < input.maxNumIters) { // not last iteration
        // beta
        beta = beta_new / beta_old;
        if (outPtr != nullptr) {
          *outPtr << ", beta: " << beta;
        }
        // PAP
        PAP = RAR - beta_new * (beta /alpha);
        TEUCHOS_TEST_FOR_EXCEPTION
          (RAR <= ATM::zero (), std::runtime_error, "At iteration " << (iter+1)
           << " out of " << input.maxNumIters << ", R.dot(AR) = " << RAR <<
           " <= 0.  This usually means that the matrix A is not symmetric "
           "(Hermitian) positive definite.");

        // alpha = 
        alpha = beta_new / PAP;
        if (outPtr != nullptr) {
          *outPtr << ", alpha: " << alpha << endl;
        }
        // beta_old
        beta_old = beta_new;
      }
    }

    // Reached max iteration count without converging
    output.absResid = r_norm;
    output.relResid = r_norm / r_norm_orig;
    output.numIters = input.maxNumIters;
    output.converged = false;

    Tpetra::deep_copy (B, R);
    return output;
  }
};

template<class SC, class MV, class OP,
         template<class, class, class> class KrylovSubclassType>
class SolverManager;

// This is the Belos::SolverManager subclass that gets registered with
// Belos::SolverFactory.
template<class SC, class MV, class OP>
using CgSingleReduceSolverManager = SolverManager<SC, MV, OP, CgSingleReduce>;

/// \brief Register CgSingleReduceSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_CgSingleReduce (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_CG_SINGLE_REDUCE_HPP
