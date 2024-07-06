// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_GMRES_PIPELINE_HPP
#define BELOS_TPETRA_GMRES_PIPELINE_HPP

#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_UpdateNewton.hpp"
#include "Tpetra_idot.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::Operator<>::scalar_type,
         class MV = Tpetra::MultiVector<SC>,
         class OP = Tpetra::Operator<SC>>
class GmresPipeline : public Gmres<SC, MV, OP> {
private:
  using base_type = Gmres<SC, MV, OP>;
  using MVT = Belos::MultiVecTraits<SC, MV>;
  using LO = typename MV::local_ordinal_type;
  using STS = Teuchos::ScalarTraits<SC>;
  using mag_type = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<mag_type>;
  using complex_type = std::complex<mag_type>;
  using dense_matrix_type = Teuchos::SerialDenseMatrix<LO, SC>;
  using dense_vector_type = Teuchos::SerialDenseVector<LO, SC>;
  using vec_type = typename Krylov<SC, MV, OP>::vec_type;
  using device_type = typename MV::device_type;
  using dot_type = typename MV::dot_type;
  using dot_view_type = Kokkos::View<dot_type*, device_type>;

public:
  GmresPipeline () :
    base_type::Gmres ()
  {
    this->input_.computeRitzValues = true;
  }

  GmresPipeline (const Teuchos::RCP<const OP>& A) :
    base_type::Gmres (A)
  {
    this->input_.computeRitzValues = true;
  }

  virtual ~GmresPipeline () = default;

private:
  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in X/out X
               vec_type& B, // in B
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input)
  {
    using std::endl;

    int restart = input.resCycle;
    int ell = 1;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();
    const mag_type eps = STS::eps ();
    const mag_type tolOrtho = mag_type (10.0) * STM::squareroot (eps);
    const bool computeRitzValues = input.computeRitzValues;

    // timers
    Teuchos::RCP< Teuchos::Time > spmvTimer = Teuchos::TimeMonitor::getNewCounter ("GmresPipeline::matrix-apply");

    // initialize output parameters
    SolverOutput<SC> output {};
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;

    Teuchos::BLAS<LO ,SC> blas;

    mag_type b_norm; // initial residual norm
    mag_type b0_norm; // initial residual norm, not left-preconditioned
    mag_type r_norm;
    mag_type r_norm_imp = -STM::one ();
    dense_matrix_type  G (restart+1, restart+1, true);
    dense_matrix_type  H (restart+1, restart, true);
    dense_vector_type  y (restart+1, true);
    dense_vector_type  h (restart+1, true);
    std::vector<mag_type> cs (restart);
    std::vector<SC> sn (restart);

    bool zeroOut = false; // Kokkos::View:init can take a long time on GPU?
    MV  Q (B.getMap (), restart+1, zeroOut);
    MV  V (B.getMap (), restart+1, zeroOut);
    vec_type R (B.getMap (), zeroOut);
    vec_type Y (B.getMap (), zeroOut);
    vec_type MZ (B.getMap (), zeroOut);
    vec_type Z0 = * (V.getVectorNonConst (0));

    // initial residual (making sure R = B - Ax)
    {
      Teuchos::TimeMonitor LocalTimer (*spmvTimer);
      A.apply (X, R);
    }
    R.update (one, B, -one);
    // TODO: this should be idot?
    b0_norm = STM::squareroot (STS::real (R.dot (R))); // initial residual norm, no preconditioned
    if (input.precoSide == "left") {
      M.apply (R, Z0);
      // TODO: this should be idot?
      b_norm = STS::real (Z0.dot( Z0 )); //Z.norm2 (); // initial residual norm, preconditioned
    }
    else {
      b_norm = b0_norm;
    }
    r_norm = b_norm;

    if (computeRitzValues) {
      // Invoke standard Gmres for the first restart cycle, to compute
      // Ritz values as Newton shifts
      SolverInput<SC> input_gmres = input;
      input_gmres.maxNumIters = input.resCycle;
      input_gmres.computeRitzValues = true;

      Tpetra::deep_copy (R, B);
      output = Gmres<SC, MV, OP>::solveOneVec (outPtr, X, R, A, M,
                                               input_gmres);
      if (output.converged) {
        return output; // standard GMRES converged
      }
      if (input.precoSide == "left") {
        M.apply (R, Z0);
        r_norm = Z0.norm2 (); // residual norm
      }
      else {
        r_norm = output.absResid;
      }
      output.numRests++;
    }
    if (input.precoSide != "left") {
      Tpetra::deep_copy (Z0, R);
    }

    // for idot
    std::shared_ptr<Tpetra::Details::CommRequest> req;
    dot_view_type vals ("results[numVecs]", restart+1);
    auto vals_h = Kokkos::create_mirror_view (vals);

    // Initialize starting vector
    //Z.scale (one / b_norm);
    G(0, 0) = r_norm*r_norm;
    y[0] = r_norm;

    // Main loop
    int iter = 0;
    mag_type metric = 2*input.tol; // to make sure to hit the first synch
    while (output.numIters < input.maxNumIters && ! output.converged) {
      if (iter == 0) {
        if (input.maxNumIters < output.numIters+restart) {
          restart = input.maxNumIters-output.numIters;
        }

        // Normalize initial vector
        MVT::MvScale (Z0, one/std::sqrt(G(0, 0)));

        // Copy initial vector
        vec_type AP = * (Q.getVectorNonConst (0));
        Tpetra::deep_copy (AP, Z0);
      }

      // Restart cycle
      for (; iter < restart+ell && metric > input.tol; ++iter) {
        if (iter < restart) {
          // W = A*Z
          vec_type Z = * (V.getVectorNonConst (iter));
          vec_type W = * (V.getVectorNonConst (iter+1));
          if (input.precoSide == "none") {
            Teuchos::TimeMonitor LocalTimer (*spmvTimer);
            A.apply (Z, W);
          }
          else if (input.precoSide == "right") {
            M.apply (Z, MZ);
            {
              Teuchos::TimeMonitor LocalTimer (*spmvTimer);
              A.apply (MZ, W);
            }
          }
          else {
            {
              Teuchos::TimeMonitor LocalTimer (*spmvTimer);
              A.apply (Z, MZ);
            }
            M.apply (MZ, W);
          }
          // Shift for Newton basis, explicitly for the first iter
          // (rest is done through change-of-basis)
          if (computeRitzValues && iter == 0) {
            //W.update (-output.ritzValues(iter%ell),  Z, one);
            const complex_type theta = output.ritzValues[iter%ell];
            UpdateNewton<SC, MV>::updateNewtonV (iter, V, theta);
          }
          output.numIters ++;
        }
        int k = iter+1 - ell; // we synch idot from k-th iteration
        if (outPtr != nullptr && k > 0) {
          *outPtr << "Current pipeline iteration: iter=" << iter
                  << ", restart=" << restart
                  << ", metric=" << metric << endl;
          Indent indent3 (outPtr);
        }

        // Compute G and H
        if (k >= 0) {
          if (k > 0) {
            req->wait (); // wait for idot
            auto v_iter = Kokkos::subview(vals, std::pair<int, int>(0, iter+1));
            auto h_iter = Kokkos::subview(vals_h, std::pair<int, int>(0, iter+1));
            Kokkos::deep_copy (h_iter, v_iter);

            for (int i = 0; i <= iter; i++) {
              G(i, k) = vals_h[i];
            }
            for (int i = 0; i <= k; ++i) {
              H(i, k-1) = G(i, k);
            }
            // Integrate shift for Newton basis (applied through
            // change-of-basis)
            if (computeRitzValues) {
              //H(k-1, k-1) += output.getRitzValue((k-1)%ell);
              const complex_type theta = output.ritzValues[(k-1)%ell];
              UpdateNewton<SC, MV>::updateNewtonH (k-1, H, theta);
            }

            // Fix H
            for (int i = 0; i < k; ++i) {
              H(k, k-1) -= (G(i, k)*G(i, k));
            }
            TEUCHOS_TEST_FOR_EXCEPTION
              (STS::real (H(k, k-1)) < STM::zero (), std::runtime_error,
               "At iteration " << iter << ", H(" << k << ", "
               << k-1 << ") = " << H(k, k-1) << " < 0.");
            H(k, k-1) = std::sqrt( H(k, k-1) );
          }

          if (k > 0) {
            // Orthogonalize V(:, k), k = iter+1-ell
            vec_type AP = * (Q.getVectorNonConst (k));
            Teuchos::Range1D index_prev (0, k-1);
            const MV Qprev = * (Q.subView (index_prev));
            dense_matrix_type g_prev (Teuchos::View, G, k, 1, 0, k);

            MVT::MvTimesMatAddMv (-one, Qprev, g_prev, one, AP);
            MVT::MvScale (AP, one/H(k, k-1));
          }
        }

        if (iter < restart) {
          // Apply change-of-basis to W
          vec_type W = * (V.getVectorNonConst (iter+1));
          if (k > 0) {
            Teuchos::Range1D index_prev (ell, iter);
            const MV Zprev = * (V.subView (index_prev));

            dense_matrix_type h_prev (Teuchos::View, H, k, 1, 0, k-1);
            MVT::MvTimesMatAddMv (-one, Zprev, h_prev, one, W);

            MVT::MvScale (W, one/H(k, k-1));
          }
        }

        if (k > 0) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (STS::real (H(k, k-1)) < STM::zero (), std::runtime_error,
             "At iteration " << k << ", H(" << k << ", " << k-1 << ") = "
             << H(k, k-1) << " < 0.");
          // NOTE (mfh 16 Sep 2018) It's not entirely clear to me
          // whether the code as given to me was correct for complex
          // arithmetic.  I'll do my best to make it compile.
          if (STS::real (H(k, k-1)) > STS::real (tolOrtho*G(k, k))) {
            // Apply Givens rotations to new column of H and y
            this->reduceHessenburgToTriangular (k-1, H, cs, sn, y.values());
            // Convergence check
            metric = this->getConvergenceMetric (STS::magnitude (y(k)),
                                                 b_norm, input);
          }
          else { // breakdown
            H(k, k-1) = zero;
            metric = STM::zero ();
          }
        }

        if (iter < restart && metric > input.tol) {
          // Copy the new vector
          vec_type AP = * (Q.getVectorNonConst (iter+1));
          Tpetra::deep_copy (AP, * (V.getVectorNonConst (iter+1)));

          // Start all-reduce to compute G(:, iter+1)
          // [Q(:,1:k-1), V(:,k:iter+1)]'*W
          Teuchos::Range1D index_prev(0, iter+1);
          const MV Qprev  = * (Q.subView(index_prev));

          dot_view_type v_iter = Kokkos::subview(vals, std::pair<int, int>(0, iter+2));
          vec_type W = * (V.getVectorNonConst (iter+1));
          req = Tpetra::idot (v_iter, Qprev, static_cast<const MV&> (W));
        }
      } // End of restart cycle

      if (iter < restart+ell) {
        // save the old solution, just in case explicit residual norm failed the convergence test
        Tpetra::deep_copy (Y, X);
        blas.COPY (1+iter, y.values(), 1, h.values(), 1);
      }
      if (iter >= ell) {
        r_norm_imp = STS::magnitude (y (iter - ell)); // save implicit residual norm

        // Update solution
        blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
                   Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
                   iter-ell, 1, one,
                   H.values(), H.stride(), y.values(), y.stride());
        Teuchos::Range1D cols(0, (iter-ell)-1);
        Teuchos::RCP<const MV> Qj = Q.subView(cols);
        y.resize (iter);
        if (input.precoSide == "right") {
          dense_vector_type y_iter (Teuchos::View, y.values (), iter-ell);
          MVT::MvTimesMatAddMv (one, *Qj, y_iter, zero, R);
          M.apply (R, MZ);
          X.update (one, MZ, one);
        }
        else {
          dense_vector_type y_iter (Teuchos::View, y.values (), iter-ell);
          MVT::MvTimesMatAddMv (one, *Qj, y_iter, one, X);
        }
        y.resize (restart+1);
      }
      // Compute real residual
      {
        Teuchos::TimeMonitor LocalTimer (*spmvTimer);
        A.apply (X, R);
      }
      R.update (one, B, -one);
      // TODO: compute residual norm with all-reduce, should be idot?
      r_norm = R.norm2 ();
      output.absResid = r_norm;
      output.relResid = r_norm / b_norm;
      if (outPtr != nullptr) {
        *outPtr << "Implicit and explicit residual norms at restart: " << r_norm_imp << ", " << r_norm << endl;
      }

      // Convergence check (with explicitly computed residual norm)
      metric = this->getConvergenceMetric (r_norm, b_norm, input);
      if (metric <= input.tol) {
        output.converged = true;
      }
      else if (output.numIters < input.maxNumIters) {
        // Restart, only if max inner-iteration was reached.
        // Otherwise continue the inner-iteration.
        if (iter >= restart+ell) {
          // Restart: Initialize starting vector for restart
          iter = 0;
          Z0 = * (V.getVectorNonConst (0));
          if (input.precoSide == "left") {
            M.apply (R, Z0);
            r_norm = STS::real (Z0.dot (Z0)); //norm2 (); // residual norm
          }
          else {
            // set the starting vector
            Tpetra::deep_copy (Z0, R);
          }
          G(0, 0) = r_norm;
          r_norm = STM::squareroot (r_norm);
          //Z.scale (one / r_norm);
          y[0] = r_norm;
          for (int i=1; i < restart+1; ++i) {
            y[i] = zero;
          }
          // Restart
          output.numRests ++;
        }
        else {
          // reset to the old solution
          Tpetra::deep_copy (X, Y);
          blas.COPY (1+iter, h.values(), 1, y.values(), 1);
          {
            // Copy the new vector
            vec_type AP = * (Q.getVectorNonConst (iter));
            Tpetra::deep_copy (AP, * (V.getVectorNonConst (iter)));

            // Start all-reduce to compute G(:, iter)
            // [Q(:,1:k-1), V(:,k:iter)]'*W
            Teuchos::Range1D index_prev(0, iter);
            const MV Qprev  = * (Q.subView(index_prev));

            vec_type W = * (V.getVectorNonConst (iter));
            req = Tpetra::idot (vals, Qprev, static_cast<const MV&> (W));
          }
        }
      }
    }

    return output;
  }
};

template<class SC, class MV, class OP,
         template<class, class, class> class KrylovSubclassType>
class SolverManager;

// This is the Belos::SolverManager subclass that gets registered with
// Belos::SolverFactory.
template<class SC, class MV, class OP>
using GmresPipelineSolverManager = SolverManager<SC, MV, OP, GmresPipeline>;

/// \brief Register GmresPipelineSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_GmresPipeline (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_GMRES_PIPELINE_HPP
