#ifndef BELOS_TPETRA_CGPIPELINE_HPP
#define BELOS_TPETRA_CGPIPELINE_HPP

#include "Belos_Tpetra_Krylov.hpp"
#include "Tpetra_idot.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::MultiVector<>::scalar_type,
	 class MV = Tpetra::MultiVector<SC>,
	 class OP = Tpetra::Operator<SC>>
class CgPipeline : public Krylov<SC, MV, OP> {
private:
  using base_type = Krylov<SC, MV, OP>;

public:
  CgPipeline () :
    base_type::Krylov ()
  {}

  CgPipeline (const Teuchos::RCP<const OP>& A) :
    base_type::Krylov (A)
  {}

  virtual ~CgPipeline () {}

protected:
  using vec_type = typename base_type::vec_type;

  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in/out
               vec_type& R_in, // in/out
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input) override
  {
    using std::endl;
    using device_type = typename MV::device_type;
    using val_type = typename MV::dot_type;
    using STS = Kokkos::ArithTraits<val_type>;
    using mag_type = typename STS::mag_type;
    using STM = Kokkos::ArithTraits<mag_type>;
    
    const auto ONE = STS::one ();

    SolverOutput<SC> output {};

    MV R_AR (R_in.getMap (), 2); // [R, A*R]
    vec_type R = * (R_AR.getVectorNonConst (0));
    vec_type AR = * (R_AR.getVectorNonConst (1));
    Tpetra::deep_copy (R, R_in);

    // aux matrix containing R and R
    MV R_R (R_in.getMap (), 2);
    vec_type R1 =  * (R_R.getVectorNonConst (0));
    vec_type R2 =  * (R_R.getVectorNonConst (1));

    // results of [R R]'*[R AR]
    Kokkos::View<val_type*, device_type> RR_RAR ("results[numVecs]", 2);
    vec_type P (R, Teuchos::Copy);
    vec_type AP (P.getMap ());
    vec_type AAR(R.getMap ());
    vec_type AW (R.getMap ());
    vec_type MR (R.getMap ());
    vec_type U (R.getMap ());
    vec_type Q (R.getMap ());

    // local vars
    val_type RAR, PAP;
    mag_type alpha = STM::zero ();
    mag_type beta = STM::zero ();
    mag_type r_norm = STM::zero ();
    mag_type r_norm_orig = STM::zero ();
    mag_type beta_new = STM::zero ();
    mag_type beta_old = STM::zero ();

    // Initial step
    if (input.precoSide == "none") {
      // W = A*R
      A.apply (R, AR);
    }
    else {
      // preconditioner
      M.apply (R, U);
      // W = A*R
      A.apply (U, AR);
    }

    // main loop
    for (int iter = 0; iter < input.maxNumIters; ++iter) {
      if (outPtr != nullptr) {
        *outPtr << "Iteration " << (iter+1) << " of " << input.maxNumIters << ":" << endl;
        outPtr->pushTab ();
        *outPtr << "r_norm: " << r_norm << endl;
      }

      // * all-reduce *
      // [RR,RAR] = [R,RAR]'*R
      using request_ptr = decltype (Tpetra::idot (RR_RAR, R_AR, R));
      request_ptr req;
      {
        if (input.precoSide == "none") {
          req = Tpetra::idot (RR_RAR, R_AR, R);
        } else {
          req = Tpetra::idot (RR_RAR, R_AR, U);
        }
      }

      // * matrix op (moved up to overlap with all-reduce) *
      if (input.precoSide == "none") {
        // n = A*w (AR is w, and AAR is n)
        A.apply (AR, AAR);
      }
      else {
        // preconditioner
        // m = M^{-1}*w (AR is w, and MR is m)
        M.apply (AR, MR);
        // n = A*m (MR is m, and AAR is n)
        A.apply (MR, AAR);
      }

      // * check for convergence *
      req->wait ();
      RAR = RR_RAR(1);
      beta_new = STS::real (RR_RAR(0));
      r_norm = std::sqrt( beta_new );
      if (iter == 0) {
        r_norm_orig = r_norm;
      }
      const mag_type metric =
	this->getConvergenceMetric (r_norm, r_norm_orig, input);
      if (outPtr != nullptr) {
        *outPtr << "RAR: " << RAR << endl;
        *outPtr << "r_norm: " << r_norm << endl;
        *outPtr << "metric: " << metric << endl;
      }
      if (metric <= input.tol) {
        output.absResid = r_norm;
        output.relResid = r_norm / r_norm_orig;
        output.numIters = iter + 1;
        output.converged = true;

        Tpetra::deep_copy (R_in, R);
        return output;
      }

      if (iter == 0) {
        alpha = beta_new / STS::real (RAR);
        beta  = STM::zero ();
      }
      else {
        // beta
        beta = beta_new / beta_old;
        if (outPtr != nullptr) {
          *outPtr << "beta: " << beta << endl;
        }
        // PAP
        PAP = RAR - beta_new * (beta / alpha);
        TEUCHOS_TEST_FOR_EXCEPTION
          (STS::real (RAR) <= STM::zero (), std::runtime_error,
	   "At iteration " << (iter+1)
           << " out of " << input.maxNumIters << ", R.dot(AR) = " << RAR <<
           " <= 0.  This usually means that the matrix A is not symmetric "
           "(Hermitian) positive definite.");

        // alpha
        alpha = beta_new / STS::real (PAP);
        if (outPtr != nullptr) {
          *outPtr << "alpha: " << alpha << endl;
        }
      }
      // beta_old
      beta_old = beta_new;

      // * search direction *
      if (input.precoSide == "none") {
        // P = R + beta*P
        P.update (ONE, R, static_cast<SC> (beta));
      }
      else {
        // p = u + beta*p (P is p, and U is u)
        P.update (ONE, U, static_cast<SC> (beta));
        // q = m + beta*q (Q is q, and MR is m)
        Q.update (ONE, MR, static_cast<SC> (beta));
        // u = u - alpha*q (U is u, and Q is q)
        U.update (static_cast<SC> (-alpha), Q, ONE);
      }
      // s = w + beta*s (AP is s, and AR is w)
      AP.update (ONE, AR, static_cast<SC> (beta));

      // * solution update *
      // x = x + alpha*p (X is x, and P is p)
      X.update (static_cast<SC> (alpha), P, ONE);
      // r = r - alpha*s (R is r, and AP is s)
      R.update (static_cast<SC> (-alpha), AP, ONE);

      // * next search directions update *
      // z = n + beta*z (AW is z, and AAR is n)
      AW.update (ONE, AAR, static_cast<SC> (beta));
      // w = w - alpha*n (AR is w, and AW is n)
      AR.update (static_cast<SC> (-alpha), AW, ONE);
    }
    if (outPtr != nullptr) {
      outPtr->popTab ();
    }

    // Reached max iteration count without converging
    Tpetra::deep_copy (R_in, R);
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
using CgPipelineSolverManager = SolverManager<SC, MV, OP, CgPipeline>;

/// \brief Register CgPipelineSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_CgPipeline (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_CGPIPELINE_HPP
