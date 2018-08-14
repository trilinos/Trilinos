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
  CgSingleReduce () :
    base_type::Krylov ()
  {}

  CgSingleReduce (const Teuchos::RCP<const OP>& A) :
    base_type::Krylov (A)
  {}

  virtual ~CgSingleReduce () {}

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
    using device_type = typename MV::device_type;
    using STS = Teuchos::ScalarTraits<SC>;
    using mag_type = typename STS::magnitudeType;
    using STM = Teuchos::ScalarTraits<mag_type>;
    using dot_type = typename MV::dot_type;
    
    const SC ONE = STS::one ();
    SolverOutput<SC> output {};

    // compute initial residual
    vec_type MR (B.getMap ());
    mag_type beta_old = STM::zero ();
    if (input.precoSide == "none") {
      beta_old = STS::real (B.dot (B));
    }
    else {
      M.apply (B, MR);
      beta_old = STS::real (B.dot (MR));
    }
    mag_type r_norm = std::sqrt (beta_old);
    mag_type r_norm_orig = r_norm;

    // quick return
    mag_type metric = this->getConvergenceMetric (r_norm, r_norm_orig, input);
    if (metric <= input.tol) {
      output.absResid = r_norm;
      output.relResid = STM::one ();
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
    mag_type RAR;
    Kokkos::View<dot_type*, device_type> RR_RAR ("results[numVecs]", 2);
    vec_type P (R, Teuchos::Copy);
    vec_type AP (P.getMap ());

    // Initial step
    // AR = A*R
    mag_type PAP;
    if (input.precoSide == "none") {
      A.apply (R, AR);
      PAP = STS::real (R.dot (AR));
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

      beta_old = STS::real (RR_RAR(0));
      PAP = STS::real (RR_RAR(1));

      r_norm = std::sqrt (beta_old);
    }
    mag_type alpha      = beta_old / PAP;
    mag_type beta       = STM::zero ();
    mag_type beta_new   = STM::zero ();
    // main loop
    for (int iter = 0; iter < input.maxNumIters; ++iter) {
      if (outPtr != nullptr) {
        *outPtr << "Iteration " << (iter+1) << " of " << input.maxNumIters
		<< ":" << endl;
        outPtr->pushTab ();
        *outPtr << "r_norm: " << r_norm << endl;
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
      beta_new = STS::real (RR_RAR(0));
      RAR = STS::real (RR_RAR(1));

      // convergence check
      r_norm = std::sqrt( beta_new );
      metric = this->getConvergenceMetric (r_norm, r_norm_orig, input);
      if (outPtr != nullptr) {
        *outPtr << "r_norm: " << r_norm << endl;
        *outPtr << "RAR: " << RAR << endl;
        *outPtr << "metric: " << metric << endl;
      }
      if (metric <= input.tol) {
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
          *outPtr << "beta: " << beta << endl;
        }
        // PAP
        PAP = RAR - beta_new * (beta /alpha);
        TEUCHOS_TEST_FOR_EXCEPTION
          (RAR <= STM::zero (), std::runtime_error, "At iteration " << (iter+1)
           << " out of " << input.maxNumIters << ", R.dot(AR) = " << RAR <<
           " <= 0.  This usually means that the matrix A is not symmetric "
           "(Hermitian) positive definite.");

        // alpha = 
        alpha = beta_new / PAP;
        if (outPtr != nullptr) {
          *outPtr << "alpha: " << alpha << endl;
        }
        // beta_old
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
