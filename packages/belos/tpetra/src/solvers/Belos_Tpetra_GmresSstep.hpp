// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_GMRES_SSTEP_HPP
#define BELOS_TPETRA_GMRES_SSTEP_HPP

#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_UpdateNewton.hpp"

#include "KokkosBlas3_trsm.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::Operator<>::scalar_type,
         class MV = Tpetra::MultiVector<>,
         class OP = Tpetra::Operator<> >
class CholQR {
private:
  using LO = typename MV::local_ordinal_type;
  using blas_type = Teuchos::BLAS<LO, SC>;
  using lapack_type = Teuchos::LAPACK<LO, SC>;

public:
  /// \typedef FactorOutput
  /// \brief Return value of \c factor().
  ///
  /// Here, FactorOutput is just a minimal object whose value is
  /// irrelevant, so that this class' interface looks like that of
  /// \c CholQR.
  using FactorOutput = int;
  using STS = Teuchos::ScalarTraits<SC>;
  using MVT = Belos::MultiVecTraits<SC, MV>;
  using mag_type = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<mag_type>;
  using dense_matrix_type = Teuchos::SerialDenseMatrix<LO, SC>;

  /// \brief Constructor
  ///
  /// \param theCacheSizeHint [in] Cache size hint in bytes.  If 0,
  ///   the implementation will pick a reasonable size, which may be
  ///   queried by calling cache_size_hint().
  CholQR () = default;

  /// \brief Compute the QR factorization of the matrix A.
  ///
  /// Compute the QR factorization of the nrows by ncols matrix A,
  /// with nrows >= ncols, stored either in column-major order (the
  /// default) or as contiguous column-major cache blocks, with
  /// leading dimension lda >= nrows.
  FactorOutput
  factor (Teuchos::FancyOStream* outPtr, MV& A, dense_matrix_type& R)
  {
    Teuchos::RCP< Teuchos::Time > factorTimer = Teuchos::TimeMonitor::getNewCounter ("CholQR::factor");
    Teuchos::TimeMonitor LocalTimer (*factorTimer);

    blas_type blas;
    lapack_type lapack;

    const SC zero = STS::zero ();
    const SC one  = STS::one ();

    // quick return
    size_t ncols = A.getNumVectors ();
    size_t nrows = A.getLocalLength ();
    if (ncols == 0 || nrows == 0) {
      return 0;
    }

    // Compute R := A^T * A, using a single BLAS call.
    // MV with "static" memory (e.g., Tpetra manages the static GPU memory pool)
    MV R_mv = impl::makeStaticLocalMultiVector (A, ncols, ncols);
    //R_mv.putScalar (STS::zero ());

    // compute R := A^T * A
    R_mv.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, one, A, A, zero);

    // Compute the Cholesky factorization of R in place, so that
    // A^T * A = R^T * R, where R is ncols by ncols upper
    // triangular.
    int info = 0;
    {
      auto R_h = R_mv.getLocalViewHost (Tpetra::Access::ReadWrite);
      int ldr = int (R_h.extent (0));
      SC *Rdata = reinterpret_cast<SC*> (R_h.data ());
      lapack.POTRF ('U', ncols, Rdata, ldr, &info);
      if (info > 0) {
        // FIXME (mfh 17 Sep 2018) Don't throw; report an error code.
        //ncols = info;
        //throw std::runtime_error("Cholesky factorization failed");
        *outPtr << "  >  POTRF( " << ncols << " ) failed with info = " << info << std::endl;
        for (size_t i=info-1; i<ncols; i++) {
          R_h(i, i) = one;
          for (size_t j=i+1; j<ncols; j++) {
            R_h(i, j) = zero;
          }
        }
      }
    }
    // Copy to the output R
    Tpetra::deep_copy (R, R_mv);
    // TODO: replace with a routine to zero out lower-triangular
    //     : not needed, but done for testing
    for (size_t i=0; i<ncols; i++) {
      for (size_t j=0; j<i; j++) {
        R(i, j) = zero;
      }
    }

    // Compute A := A * R^{-1}.  We do this in place in A, using
    // BLAS' TRSM with the R factor (form POTRF) stored in the upper
    // triangle of R.

    // Compute A_cur / R (Matlab notation for A_cur * R^{-1}) in place.
    {
      auto A_d = A.getLocalViewDevice (Tpetra::Access::ReadWrite);
      auto R_d = R_mv.getLocalViewDevice (Tpetra::Access::ReadOnly);
      KokkosBlas::trsm ("R", "U", "N", "N",
                        one, R_d, A_d);
    }
    return (info > 0 ? info-1 : ncols);
  }

  // recursive call to factor
  FactorOutput
  reFactor (Teuchos::FancyOStream* outPtr, MV& A, dense_matrix_type& R)
  {
    int ncols = int (A.getNumVectors ());
    int rank = 0;
    int old_rank = -1;

    // recursively call factor while cols remaining and has made progress
    while (rank < ncols && old_rank != rank) {
      Teuchos::Range1D next_index(rank, ncols-1);
      MV nextA = * (A.subView(next_index));

      dense_matrix_type nextR (Teuchos::View, R, ncols-rank, ncols-rank, rank, rank);
      old_rank = rank;
      auto new_rank = factor (outPtr, nextA, nextR);
      if (outPtr != nullptr) {
        if (rank > 0) {
          *outPtr << "  ++ reCholQR(";
        } else {
          *outPtr << "  >>   CholQR(";
        }
        *outPtr << rank << ":" << ncols-1 << "), new_rank = " << new_rank << std::endl;
      }
      rank += new_rank;
    }
    return rank;
  }
};


template<class SC = Tpetra::Operator<>::scalar_type,
         class MV = Tpetra::MultiVector<SC>,
         class OP = Tpetra::Operator<SC>>
class GmresSstep : public Gmres<SC, MV, OP>  {
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
  GmresSstep () :
    base_type::Gmres (),
    useCholQR2_ (false),
    cholqr_ (Teuchos::null)
  {}

  GmresSstep (const Teuchos::RCP<const OP>& A) :
    base_type::Gmres (A),
    cholqr_ (Teuchos::null)
  {}

  virtual ~GmresSstep () = default;

  virtual void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const
  {
    base_type::getParameters (params, defaultValues);

    const int stepSize = defaultValues ? 5 : this->input_.stepSize;
    params.set ("Step Size", stepSize );
  }

  virtual void
  setParameters (Teuchos::ParameterList& params) {
    base_type::setParameters (params);
    int stepSize = params.get<int> ("Step Size", this->input_.stepSize);
    this->input_.stepSize = stepSize;

    bool computeRitzValuesOnFly 
      = params.get<bool> ("Compute Ritz Values on Fly", this->input_.computeRitzValuesOnFly);
    this->input_.computeRitzValuesOnFly = computeRitzValuesOnFly;

    constexpr bool useCholQR_default = true;
    bool useCholQR = params.get<bool> ("CholeskyQR", useCholQR_default);

    bool useCholQR2 = params.get<bool> ("CholeskyQR2", useCholQR2_);
    useCholQR2_ = useCholQR2;

    if ((!useCholQR && !useCholQR2) && !cholqr_.is_null ()) {
      cholqr_ = Teuchos::null;
    } else if ((useCholQR || useCholQR2) && cholqr_.is_null ()) {
      cholqr_ = Teuchos::rcp (new CholQR<SC, MV, OP> ());
    }
  }

private:
  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in X/out X
               vec_type& B, // in B/out R (not left-preconditioned)
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input)
  {
    using std::endl;
    int stepSize = input.stepSize;
    int restart = input.resCycle;
    int step = stepSize;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();

    // timers
    Teuchos::RCP< Teuchos::Time > spmvTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSstep::matrix-apply");
    Teuchos::RCP< Teuchos::Time > bortTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSstep::BOrtho");
    Teuchos::RCP< Teuchos::Time > tsqrTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSstep::TSQR");

    // initialize output parameters
    SolverOutput<SC> output {};
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;

    if (outPtr != nullptr) {
      *outPtr << "GmresSstep" << endl;
      Indent indent1 (outPtr);
      *outPtr << "Solver input:" << endl;
      Indent indentInner (outPtr);
      *outPtr << input;
    }

    Teuchos::BLAS<LO, SC> blas;
    Teuchos::LAPACK<LO, SC> lapack;
    dense_matrix_type  H (restart+1, restart, true); // Hessenburg matrix
    dense_matrix_type  T (restart+1, restart, true); // H reduced to upper-triangular matrix
    dense_matrix_type  G (restart+1, step+1, true);  // Upper-triangular matrix from ortho process
    dense_matrix_type  G2(restart+1, restart, true); // a copy of Hessenburg matrix for computing Ritz values
    dense_vector_type  y (restart+1, true);
    dense_matrix_type  h (restart+1, 1, true); // used for reorthogonalization
    std::vector<mag_type> cs (restart);
    std::vector<SC> sn (restart);

    mag_type b_norm;  // initial residual norm
    mag_type b0_norm; // initial residual norm, not left preconditioned
    mag_type r_norm;
    mag_type r_norm_imp;
    mag_type metric;

    bool zeroOut = false; // Kokkos::View:init can take a long time on GPU?
    vec_type R (B.getMap (), zeroOut);
    vec_type Y (B.getMap (), zeroOut);
    vec_type MP (B.getMap (), zeroOut);
    MV  Q (B.getMap (), restart+1, zeroOut);
    vec_type P0 = * (Q.getVectorNonConst (0));

    // Compute initial residual (making sure R = B - Ax)
    {
      Teuchos::TimeMonitor LocalTimer (*spmvTimer);
      A.apply (X, R);
    }
    R.update (one, B, -one);
    b0_norm = R.norm2 (); // initial residual norm, not preconditioned
    if (input.precoSide == "left") {
      M.apply (R, P0);
      r_norm = P0.norm2 (); // initial residual norm, left-preconditioned
    } else {
      r_norm = b0_norm;
    }
    b_norm = r_norm;

    metric = this->getConvergenceMetric (b0_norm, b0_norm, input);
    if (metric <= input.tol) {
      if (outPtr != nullptr) {
        *outPtr << "Initial guess' residual norm " << b_norm
                << " meets tolerance " << input.tol << endl;
      }
      output.absResid = r_norm;
      output.relResid = r_norm / b0_norm;
      output.converged = true;
      // Return residual norm as B
      Tpetra::deep_copy (B, P0);
      return output;
    } else if (STM::isnaninf (metric)) {
      if (outPtr != nullptr) {
        *outPtr << "Initial guess' residual norm " << b_norm
                << " is nan " << endl;
      }
      output.absResid = r_norm;
      output.relResid = r_norm / b0_norm;
      output.converged = false;
      // Return residual norm as B
      Tpetra::deep_copy (B, P0);
      return output;
    } else if (input.computeRitzValues && !input.computeRitzValuesOnFly) {
      // Invoke standard Gmres for the first restart cycle, to compute
      // Ritz values for use as Newton shifts
      if (outPtr != nullptr) {
        *outPtr << "Run standard GMRES for first restart cycle" << endl;
      }
      SolverInput<SC> input_gmres = input;
      input_gmres.maxNumIters = input.resCycle;
      input_gmres.maxNumIters = std::min(input.resCycle, input.maxNumIters);
      input_gmres.computeRitzValues = true;

      Tpetra::deep_copy (R, B);
      output = Gmres<SC, MV, OP>::solveOneVec (outPtr, X, R, A, M,
                                               input_gmres);
      if (output.converged) {
        return output; // standard GMRES converged
      }

      if (input.precoSide == "left") {
        M.apply (R, P0);
        r_norm = P0.norm2 (); // residual norm
      }
      else {
        r_norm = output.absResid;
      }
      output.numRests++;
    }

    // Initialize starting vector
    if (input.precoSide != "left") {
      Tpetra::deep_copy (P0, R);
    }
    P0.scale (one / r_norm);
    y[0] = r_norm;

    // Main loop
    int iter = 0;
    while (output.numIters < input.maxNumIters && ! output.converged) {
      if (outPtr != nullptr) {
        *outPtr << "Restart cycle " << output.numRests << ":" << endl;
        Indent indent2 (outPtr);
        *outPtr << output;
      }

      if (input.maxNumIters < output.numIters+restart) {
        restart = input.maxNumIters-output.numIters;
      }

      // Restart cycle
      for (; iter < restart && metric > input.tol; iter+=step) {
        if (outPtr != nullptr) {
          *outPtr << "Current s-step iteration: iter=" << iter
                  << ", restart=" << restart
                  << ", step=" << step
                  << ", metric=" << metric << endl;
          Indent indent3 (outPtr);
        }

        // Compute matrix powers
        if (input.computeRitzValuesOnFly && output.numIters < input.stepSize) {
          stepSize = 1;
        } else {
          stepSize = input.stepSize;
        }
        for (step=0; step < stepSize && iter+step < restart; step++) {
          // AP = A*P
          vec_type P  = * (Q.getVectorNonConst (iter+step));
          vec_type AP = * (Q.getVectorNonConst (iter+step+1));
          if (input.precoSide == "none") {
            Teuchos::TimeMonitor LocalTimer (*spmvTimer);
            A.apply (P, AP);
          }
          else if (input.precoSide == "right") {
            M.apply (P, MP);
            {
              Teuchos::TimeMonitor LocalTimer (*spmvTimer);
              A.apply (MP, AP);
            }
          }
          else {
            {
              Teuchos::TimeMonitor LocalTimer (*spmvTimer);
              A.apply (P, MP);
            }
            M.apply (MP, AP);
          }
          // Shift for Newton basis
          if ( int (output.ritzValues.size()) > step) {
            //AP.update (-output.ritzValues(step), P, one);
            const complex_type theta = output.ritzValues[step];
            UpdateNewton<SC, MV>::updateNewtonV(iter+step, Q, theta);
          }
          output.numIters++;
        }

        // Orthogonalization
        {
          Teuchos::TimeMonitor LocalTimer (*bortTimer);
          this->projectBelosOrthoManager (iter, step, Q, G);
        }
        int rank = 0;
        {
          Teuchos::TimeMonitor LocalTimer (*tsqrTimer);
          rank = recursiveCholQR (outPtr, iter, step, Q, G);
          if (useCholQR2_) {
            rank = recursiveCholQR (outPtr, iter, step, Q, G2);
            // merge R 
            dense_matrix_type Rfix (Teuchos::View, G2, step+1, step+1, iter, 0);
            dense_matrix_type Rold (Teuchos::View, G,  step+1, step+1, iter, 0);
            blas.TRMM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
                       Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
                       step+1, step+1,
                       one, Rfix.values(), Rfix.stride(),
                            Rold.values(), Rold.stride());
          }
          if (rank == 0) {
            // FIXME: Don't throw; report an error code.
            throw std::runtime_error("orthogonalization failed with rank = 0");
          }
        }
        updateHessenburg (iter, step, output.ritzValues, H, G);

        // Convergence check
        if (rank == step+1 && H(iter+step, iter+step-1) != zero) {
          // Copy H to T and apply Givens rotations to new columns of T and y
          for (int iiter = 0; iiter < step; iiter++) {
            // Check negative norm
            TEUCHOS_TEST_FOR_EXCEPTION
              (STS::real (H(iter+iiter+1, iter+iiter)) < STM::zero (),
               std::runtime_error, "At iteration " << output.numIters << ", H("
               << iter+iiter+1 << ", " << iter+iiter << ") = "
               << H(iter+iiter+1, iter+iiter) << " < 0.");

            for (int i = 0; i <= iter+iiter+1; i++) {
              T(i, iter+iiter) = H(i, iter+iiter);
            }
            this->reduceHessenburgToTriangular(iter+iiter, T, cs, sn, y);
            metric = this->getConvergenceMetric (STS::magnitude (y(iter+iiter+1)), b_norm, input);
            if (outPtr != nullptr) {
              *outPtr << " > implicit residual norm=(" << iter+iiter+1 << ")="
                      << STS::magnitude (y(iter+iiter+1))
                      << " metric=" << metric << endl;
            }
            if (STM::isnaninf (metric) || metric <= input.tol) {
              if (outPtr != nullptr) {
                *outPtr << " > break at step = " << iiter+1 << " (" << step << ")" << endl;
              }
              step = iiter+1;
              break;
            }
          }
          if (STM::isnaninf (metric)) {
              // metric is nan
              break;
          }
        }
        else {
          metric = STM::zero ();
        }

        // Optionally, compute Ritz values for generating Newton basis
        if (input.computeRitzValuesOnFly && int (output.ritzValues.size()) == 0
            && output.numIters >= input.stepSize) {
          for (int i = 0; i < input.stepSize; i++) {
            for (int iiter = 0; iiter < input.stepSize; iiter++) {
              G2(i, iiter) = H(i, iiter);
            }
          }
          computeRitzValues (input.stepSize, G2, output.ritzValues);
          sortRitzValues <LO, SC> (input.stepSize, output.ritzValues);
          if (outPtr != nullptr) {
            *outPtr << " > ComputeRitzValues: " << endl;
            for (int i = 0; i < input.stepSize; i++) {
              *outPtr << " > ritzValues[ " << i << " ] = " << output.ritzValues[i] << endl;
            }
          }
        }
      } // End of restart cycle

      // Update solution
      if (iter < restart) {
        // save the old solution, just in case explicit residual norm failed the convergence test
        Tpetra::deep_copy (Y, X);
        blas.COPY (1+iter, y.values(), 1, h.values(), 1);
      }
      r_norm_imp = STS::magnitude (y (iter)); // save implicit residual norm
      blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
                 Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
                 iter, 1, one,
                 T.values(), T.stride(), y.values(), y.stride());
      Teuchos::Range1D cols(0, iter-1);
      Teuchos::RCP<const MV> Qj = Q.subView(cols);
      dense_vector_type y_iter (Teuchos::View, y.values (), iter);
      if (input.precoSide == "right") {
        MVT::MvTimesMatAddMv (one, *Qj, y_iter, zero, R);
        M.apply (R, MP);
        X.update (one, MP, one);
      }
      else {
        MVT::MvTimesMatAddMv (one, *Qj, y_iter, one, X);
      }
      // Compute real residual (not-preconditioned)
      {
        Teuchos::TimeMonitor LocalTimer (*spmvTimer);
        A.apply (X, R);
      }
      R.update (one, B, -one);
      r_norm = R.norm2 (); // residual norm
      output.absResid = r_norm;
      output.relResid = r_norm / b0_norm;
      if (outPtr != nullptr) {
        *outPtr << "Implicit and explicit residual norms at restart: " << r_norm_imp << ", " << r_norm << endl;
      }

      metric = this->getConvergenceMetric (r_norm, b0_norm, input);
      if (metric <= input.tol) {
        // update solution
        output.converged = true;
      }
      else if (STM::isnaninf (metric)) {
        // failed with nan
        // Return residual norm as B
        Tpetra::deep_copy (B, R);
        return output;
      }
      else if (output.numIters < input.maxNumIters) {
        // Restart, only if max inner-iteration was reached.
        // Otherwise continue the inner-iteration.
        if (iter >= restart) {
          // Restart: Initialize starting vector for restart
          iter = 0;
          P0 = * (Q.getVectorNonConst (0));
          if (input.precoSide == "left") { // left-precond'd residual norm
            M.apply (R, P0);
            r_norm = P0.norm2 ();
          }
          else {
            // set the starting vector
            Tpetra::deep_copy (P0, R);
          }
          P0.scale (one / r_norm);
          y[0] = SC {r_norm};
          for (int i=1; i < restart+1; ++i) {
            y[i] = STS::zero ();
          }
          output.numRests++;
        }
        else {
          // reset to the old solution
          Tpetra::deep_copy (X, Y);
          blas.COPY (1+iter, h.values(), 1, y.values(), 1);
        }
      }
    }

    // Return residual norm as B
    Tpetra::deep_copy (B, R);

    if (outPtr != nullptr) {
      *outPtr << "At end of solve:" << endl;
      Indent indentInner (outPtr);
      *outPtr << output;
    }
    return output;
  }

protected:
  void
  updateHessenburg (const int n,
                    const int s,
                    std::vector<complex_type>& S,
                    dense_matrix_type& H,
                    dense_matrix_type& R) const
  {
    const SC one  = STS::one ();
    const SC zero = STS::zero ();

    // 1) multiply H with R(1:n+s+1, 1:n+s+1) from left
    //    where R(1:n, 1:n) = I and 
    //          H(n+1:n+s+1, n+1:n+s) = 0, except h(n+j+1,n+j)=1 for j=1,..,s
    // 1.1) copy: H(j:n-1, j:n-1) = R(j:n-1, j:n-1), i.e., H = R*B
    for (int j = 0; j < s; j++ ) {
      for (int i = 0; i <= n+j+1; i++) {
        H(i, n+j) = R(i, j+1);
        if (int (S.size ()) > j && i <= n+j) {
          //H(i, n+j) += S[j].real * R(i, j);
          H(i, n+j) += UpdateNewton<SC, MV>::updateNewtonH (i, j, R, S[j]);
        }
      }
      for(int i = n+j+2; i <= n+s; i++) {
        H(i, n+j) = zero;
      }
    }

    // submatrices
    dense_matrix_type r_diag (Teuchos::View, R, s,   s, n, 0);
    dense_matrix_type h_diag (Teuchos::View, H, s+1, s, n, n);
    Teuchos::BLAS<LO, SC> blas;

    if (n == 0) { // >> first matrix-power iteration <<
      // 2) multiply H with R(1:s, 1:s)^{-1} from right
      // diagonal block
      blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                 Teuchos::NON_UNIT_DIAG, s+1, s, one,
                 r_diag.values(), r_diag.stride(),
                 h_diag.values(), h_diag.stride());
    } else  { // >> rest of iterations <<
      // 1.2) update the starting vector
      for (int i = 0; i < n; i++ ) {
        H(i, n-1) += H(n, n-1) * R(i, 0);
      }
      H(n, n-1) *= R(n, 0);

      // 2) multiply H with R(1:n+s, 1:n+s)^{-1} from right,
      //    where R(1:n, 1:n) = I
      // 2.1) diagonal block
      for (int j = 0; j < s; j++ ) {
        H(n, n+j) -= H(n, n-1) * R(n-1, j);
      }
      // diagonal block
      blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                 Teuchos::NON_UNIT_DIAG, s, s, one,
                 r_diag.values(), r_diag.stride(),
                 h_diag.values(), h_diag.stride());
      H(n+s, n+s-1) /= R(n+s-1, s-1);

      // 2.2) upper off-diagonal block: H(0:j-1, j:j+n-2)
      dense_matrix_type r_off (Teuchos::View, R, n, s, 0, 0);
      dense_matrix_type h_off (Teuchos::View, H, n, s, 0, n);
      blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                n, s, n,
               -one, H.values(),      H.stride(),
                     r_off.values(), r_off.stride(),
                one, h_off.values(), h_off.stride());

      blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                Teuchos::NON_UNIT_DIAG, n, s, one,
                r_diag.values(), r_diag.stride(),
                h_off.values(),  h_off.stride() );
    }
  }

  //! Apply the orthogonalization using Belos' OrthoManager
  int
  normalizeCholQR (Teuchos::FancyOStream* outPtr,
                   const int n,
                   const int s,
                   MV& Q,
                   dense_matrix_type& R)
  {
    // vector to be orthogonalized
    Teuchos::Range1D index_prev(n, n+s);
    MV Qnew = * (Q.subView(index_prev));

    dense_matrix_type r_new (Teuchos::View, R, s+1, s+1, n, 0);

    int rank = 0;
    if (cholqr_ != Teuchos::null) {
      rank = cholqr_->factor (outPtr, Qnew, r_new);
    }
    else {
      rank = this->normalizeBelosOrthoManager (Qnew, r_new);
    }
    return rank;
  }

  //! Apply the orthogonalization using Belos' OrthoManager
  int
  recursiveCholQR (Teuchos::FancyOStream* outPtr,
                   const int n,
                   const int s,
                   MV& Q,
                   dense_matrix_type& R)
  {
    // vector to be orthogonalized
    Teuchos::Range1D index_prev(n, n+s);
    MV Qnew = * (Q.subView(index_prev));

    dense_matrix_type r_new (Teuchos::View, R, s+1, s+1, n, 0);

    int rank = 0;
    if (cholqr_ != Teuchos::null) {
      rank = cholqr_->reFactor (outPtr, Qnew, r_new);
    }
    else {
      rank = this->normalizeBelosOrthoManager (Qnew, r_new);
    }
    return rank;
  }

private:
  bool useCholQR2_;
  Teuchos::RCP<CholQR<SC, MV, OP> > cholqr_;
};


template<class SC, class MV, class OP,
         template<class, class, class> class KrylovSubclassType>
class SolverManager;

// This is the Belos::SolverManager subclass that gets registered with
// Belos::SolverFactory.
template<class SC, class MV, class OP>
using GmresSstepSolverManager = SolverManager<SC, MV, OP, GmresSstep>;

/// \brief Register GmresSstepSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_GmresSstep (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_GMRES_SSTEP_HPP
