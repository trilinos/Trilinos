// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_CG_HPP
#define BELOS_TPETRA_CG_HPP

#include "Belos_Tpetra_Krylov.hpp"
#include "Tpetra_idot.hpp"

namespace BelosTpetra {
namespace Impl {  

template<class SC = Tpetra::MultiVector<>::scalar_type,
	 class MV = Tpetra::MultiVector<SC>,
	 class OP = Tpetra::Operator<SC>>
class Cg: public Krylov <SC, MV, OP> {
private:
  using base_type = Krylov<SC, MV, OP>;

public:
  Cg () = default;

  Cg (const Teuchos::RCP<const OP>& A) :
    base_type::Krylov (A)
  {}

  virtual ~Cg () = default;

protected:
  using vec_type = typename base_type::vec_type;
  
  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in X/out X
               vec_type& B, // in B/out R
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input)
  {
    using std::endl;
    using dev_type = typename MV::device_type;
    using val_type = typename MV::dot_type;
    using ATS = Kokkos::ArithTraits<val_type>;
    using magnitude_type = typename ATS::mag_type;
    using ATM = Kokkos::ArithTraits<magnitude_type>;

    const auto ONE  = ATS::one ();

    SolverOutput<SC> output {};

    // scalars
    magnitude_type PAP;
    magnitude_type alpha;

    vec_type P  (B.getMap ());
    vec_type AP (B.getMap ());
    MV R_MR (B.getMap (), 2);
    vec_type R = * (R_MR.getVectorNonConst (0));
    vec_type MR = * (R_MR.getVectorNonConst (1));

    // compute initial residual
    Kokkos::View<val_type*, dev_type> r_beta ("results[numVecs]", 2);
    magnitude_type beta_old = ATM::zero ();
    magnitude_type r_norm = ATM::zero ();
    magnitude_type r_norm_orig = ATM::zero ();

    A.apply (X, R);
    R.update (ONE, B, -ONE);

    if (input.precoSide == "none") { // no preconditioner
      Tpetra::deep_copy (P, R);
      beta_old = ATS::real (R.dot (R));
      r_norm = beta_old;
    }
    else {
      M.apply (R, MR);
      Tpetra::deep_copy (P, MR);

      // Compute [MR, R]'*[R],
      // R'*R is used for convergence check
      //TODO: idot is used for now.
      auto req = Tpetra::idot (r_beta, R_MR, R);
      req->wait ();
      r_norm = ATS::real (r_beta(0));
      beta_old = ATS::real (r_beta(1));
    }
    r_norm = std::sqrt (r_norm);
    r_norm_orig = r_norm;

    // quick-return
    magnitude_type metric =
      this->getConvergenceMetric (r_norm, r_norm_orig, input);
    if (metric <= input.tol) {
      if (outPtr != nullptr) {
        *outPtr << "Initial guess' residual norm " << r_norm
                << " meets tolerance " << input.tol << endl;
      }
      output.absResid = r_norm;
      output.relResid = ATM::one ();
      output.numIters = 0;
      output.converged = true;
      return output;
    }

    // main loop
    magnitude_type beta_new = ATM::zero ();
    for (int iter = 0; iter < input.maxNumIters; ++iter) {
      if (outPtr != nullptr) {
        *outPtr << "Iteration " << (iter+1) << " of " << input.maxNumIters << ":" << endl;
        outPtr->pushTab ();
        *outPtr << "r_norm: " << r_norm << endl;
      }

      A.apply (P, AP); // AP = A*P
      PAP = ATS::real (P.dot (AP));

      if (outPtr != nullptr) {
        *outPtr << "PAP: " << PAP << endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION
        (ATS::real (PAP) <= ATM::zero (), std::runtime_error,
	 "At iteration " << (iter+1) << " out of " << input.maxNumIters
         << ", P.dot(AP) = " << PAP << " <= 0.  This usually means that "
	 "the matrix A is not symmetric (Hermitian) positive definite.");

      alpha = beta_old / PAP;
      if (outPtr != nullptr) {
        *outPtr << "alpha: " << alpha << endl;
      }

      X.update (static_cast<SC> (alpha), P, ONE); // X = X + alpha*P
      R.update (static_cast<SC> (-alpha), AP, ONE); // R = R - alpha*AP

      if (input.precoSide == "none") { // no preconditioner
        beta_new = ATS::real (R.dot (R));
        r_norm = beta_new;
      }
      else {
        M.apply (R, MR);
        //TODO: idot is used to compute [MR, R]'*[R] for now.
        auto req = Tpetra::idot (r_beta, R_MR, R);
        req->wait ();
        r_norm = ATS::real (r_beta(0));
        beta_new = ATS::real (r_beta(1));
      }
      r_norm = std::sqrt (r_norm);
      metric = this->getConvergenceMetric (r_norm, r_norm_orig, input);
      if (outPtr != nullptr) {
        *outPtr << "r_norm: " << r_norm << endl;
        *outPtr << "metric: " << metric << endl;
      }

      // convergence check
      if (metric <= input.tol) {
        output.absResid = r_norm;
        output.relResid = r_norm / r_norm_orig;
        output.numIters = iter + 1;
        output.converged = true;
        return output;
      }
      else if (iter + 1 < input.maxNumIters) { // not last iteration
        const magnitude_type beta = beta_new / beta_old;
        if (input.precoSide == "none") {
          P.update (ONE, R, static_cast<SC> (beta)); // P += beta*R
        }
	else {
          P.update (ONE, MR, static_cast<SC> (beta));
        }
        beta_old = beta_new;
      }

      if (outPtr != nullptr) {
        outPtr->popTab ();
      }
    }

    // Reached max iteration count without converging
    output.absResid = r_norm;
    output.relResid = r_norm / r_norm_orig;
    output.numIters = input.maxNumIters;
    output.converged = false;
    return output;
  }
};

template<class SC, class MV, class OP,
	 template<class, class, class> class KrylovSubclassType>
class SolverManager;

// This is the Belos::SolverManager subclass that gets registered with
// Belos::SolverFactory.
template<class SC, class MV, class OP>
using CgSolverManager = SolverManager<SC, MV, OP, Cg>;

/// \brief Register CgSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_Cg (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_CG_HPP
