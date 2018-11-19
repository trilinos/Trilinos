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

  virtual ~GmresPipeline ()
  {}

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

    // initialize output parameters
    SolverOutput<SC> output {};
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;

    Teuchos::BLAS<LO ,SC> blas;

    mag_type b_norm; // initial residual norm
    mag_type b0_norm; // initial residual norm, not left-preconditioned
    mag_type r_norm;
    dense_matrix_type  G (restart+1, restart+1, true);
    dense_matrix_type  H (restart+1, restart, true);
    dense_vector_type  y (restart+1, true);
    std::vector<mag_type> cs (restart);
    std::vector<SC> sn (restart);
    MV  Q (B.getMap (), restart+1);
    MV  V (B.getMap (), restart+1);
    vec_type Z = * (V.getVectorNonConst (0));
    vec_type R (B.getMap ());
    vec_type MZ (B.getMap ());

    // initial residual (making sure R = B - Ax)
    A.apply (X, R);
    R.update (one, B, -one);
    // TODO: this should be idot?
    b0_norm = STM::squareroot (STS::real (R.dot (R))); // initial residual norm, no preconditioned
    if (input.precoSide == "left") {
      M.apply (R, Z);
      // TODO: this should be idot?
      b_norm = STS::real (Z.dot( Z )); //Z.norm2 (); // initial residual norm, preconditioned
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
        M.apply (R, Z);
        r_norm = Z.norm2 (); // residual norm
      }
      else {
        r_norm = output.absResid;
      }
      output.numRests++;
    }
    if (input.precoSide != "left") {
      Tpetra::deep_copy (Z, R);
    }

    // for idot
    std::shared_ptr<Tpetra::Details::CommRequest> req;
    Kokkos::View<dot_type*, device_type> vals ("results[numVecs]",
                                               restart+1);
    auto vals_h = Kokkos::create_mirror_view (vals);

    // Initialize starting vector
    //Z.scale (one / b_norm);
    G(0, 0) = r_norm*r_norm;
    y[0] = r_norm;

    // Main loop
    mag_type metric = 2*input.tol; // to make sure to hit the first synch
    while (output.numIters < input.maxNumIters && ! output.converged) {
      int iter = 0;
      if (input.maxNumIters < output.numIters+restart) {
        restart = input.maxNumIters-output.numIters;
      }

      // Normalize initial vector
      MVT::MvScale (Z, one/std::sqrt(G(0, 0)));

      // Copy initial vector
      vec_type AP = * (Q.getVectorNonConst (0));
      Tpetra::deep_copy (AP, Z);

      // Restart cycle
      for (iter = 0; iter < restart+ell && metric > input.tol; ++iter) {
        if (iter < restart) {
          // W = A*Z
          vec_type Z = * (V.getVectorNonConst (iter));
          vec_type W = * (V.getVectorNonConst (iter+1));
          if (input.precoSide == "none") {
            A.apply (Z, W);
          }
          else if (input.precoSide == "right") {
            M.apply (Z, MZ);
            A.apply (MZ, W);
          }
          else {
            A.apply (Z, MZ);
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

        // Compute G and H
        if (k >= 0) {
          if (k > 0) {
            req->wait (); // wait for idot
            Kokkos::deep_copy (vals_h, vals);

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
              UpdateNewton<SC, MV>::updateNewtonH(k-1, H, theta);
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
            Teuchos::Range1D index_prev(0, k-1);
            const MV Qprev = * (Q.subView(index_prev));
            dense_matrix_type g_prev (Teuchos::View, G, k, 1, 0, k);

            MVT::MvTimesMatAddMv (-one, Qprev, g_prev, one, AP);
            MVT::MvScale (AP, one/H(k, k-1));
          }
        }

        if (iter < restart) {
          // Apply change-of-basis to W
          vec_type W = * (V.getVectorNonConst (iter+1));
          if (k > 0) {
            Teuchos::Range1D index_prev(ell, iter);
            const MV Zprev = * (V.subView(index_prev));

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

          vec_type W = * (V.getVectorNonConst (iter+1));
          req = Tpetra::idot (vals, Qprev, W);
        }
      } // End of restart cycle
      if (iter > 0) {
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

          //MVT::MvTimesMatAddMv (one, *Qj, y, zero, R);
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
      Z = * (V.getVectorNonConst (0));
      A.apply (X, Z);
      Z.update (one, B, -one);
      r_norm = Z.norm2 (); // residual norm
      output.absResid = r_norm;
      output.relResid = r_norm / b_norm;
      // Convergence check (with explicitly computed residual norm)
      metric = this->getConvergenceMetric (r_norm, b_norm, input);
      if (metric <= input.tol) {
        output.converged = true;
      }
      else if (output.numIters < input.maxNumIters) {
        // Initialize starting vector for restart
        if (input.precoSide == "left") {
          Tpetra::deep_copy (R, Z);
          M.apply (R, Z);
        }
        // TODO: recomputing all-reduce, should be idot?
        r_norm = STS::real (Z.dot (Z)); //norm2 (); // residual norm
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
