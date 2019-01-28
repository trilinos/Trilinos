#ifndef BELOS_TPETRA_GMRES_SSTEP_HPP
#define BELOS_TPETRA_GMRES_SSTEP_HPP

#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_UpdateNewton.hpp"

// #include "wrapTpetraQR.hpp"
// #include "wrapTpetraCholQR.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::Operator<>::scalar_type,
         class MV = Tpetra::MultiVector<>,
         class OP = Tpetra::Operator<> >
class CholQR {
private:
  using LO = typename MV::local_ordinal_type;
  typedef Teuchos::BLAS<LO, SC> blas_type;
  typedef Teuchos::LAPACK<LO, SC> lapack_type;

public:
  /// \typedef FactorOutput
  /// \brief Return value of \c factor().
  ///
  /// Here, FactorOutput is just a minimal object whose value is
  /// irrelevant, so that this class' interface looks like that of
  /// \c CholQR.
  typedef int FactorOutput;
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType mag_type;
  typedef Teuchos::SerialDenseMatrix<LO, SC> dense_matrix_type;
  typedef Belos::MultiVecTraits<SC, MV> MVT;

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
  factor (MV& A, dense_matrix_type& R)
  {
    blas_type blas;
    lapack_type lapack;

    const SC zero = STS::zero ();
    const SC one  = STS::one ();

    LO ncols = A.getNumVectors ();
    LO nrows = A.getLocalLength ();

    // Compute R := A^T * A, using a single BLAS call.
    MVT::MvTransMv(one, A, A, R);

    // Compute the Cholesky factorization of R in place, so that
    // A^T * A = R^T * R, where R is ncols by ncols upper
    // triangular.
    int info = 0;
    lapack.POTRF ('U', ncols, R.values (), R.stride(), &info);
    if (info < 0) {
      ncols = -info;
      // FIXME (mfh 17 Sep 2018) Don't throw; report an error code.
      throw std::runtime_error("Cholesky factorization failed");
    }
    // TODO: replace with a routine to zero out lower-triangular
    //     : not needed, but done for testing
    for (int i=0; i<ncols; i++) {
      for (int j=0; j<i; j++) {
        R(i, j) = zero;
      }
    }

    // Compute A := A * R^{-1}.  We do this in place in A, using
    // BLAS' TRSM with the R factor (form POTRF) stored in the upper
    // triangle of R.

    // Compute A_cur / R (Matlab notation for A_cur * R^{-1}) in place.
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    A.template sync<Kokkos::HostSpace> ();
    A.template modify<Kokkos::HostSpace> ();
    auto A_lcl = A.template getLocalView<Kokkos::HostSpace> ();
#else
    A.sync_host ();
    A.modify_host ();
    auto A_lcl = A.getLocalViewHost ();
#endif
    SC* const A_lcl_raw = reinterpret_cast<SC*> (A_lcl.data ());
    const LO LDA = LO (A.getStride ());

    blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI,
               Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
               nrows, ncols, one, R.values(), R.stride(),
               A_lcl_raw, LDA);
    A.template sync<typename MV::device_type::memory_space> ();

    return (info > 0 ? info : ncols);
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
    stepSize_ (1),
    tsqr_ (Teuchos::null)
  {
    this->input_.computeRitzValues = true;
  }

  GmresSstep (const Teuchos::RCP<const OP>& A) :
    base_type::Gmres (A),
    stepSize_ (1),
    tsqr_ (Teuchos::null)
  {
    this->input_.computeRitzValues = true;
  }

  virtual ~GmresSstep ()
  {}

  virtual void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const
  {
    base_type::getParameters (params, defaultValues);

    const int stepSize = defaultValues ? 100 : stepSize_;

    params.set ("Step Size", stepSize );
  }

  virtual void
  setParameters (Teuchos::ParameterList& params) {
    base_type::setParameters (params);
    constexpr bool useCholQR_default = true;

    int stepSize = stepSize_;
    if (params.isParameter ("Step Size")) {
      stepSize = params.get<int> ("Step Size");
    }

    bool useCholQR = useCholQR_default;
    if (params.isParameter ("CholeskyQR")) {
      useCholQR = params.get<bool> ("CholeskyQR");
    }

    if (useCholQR && tsqr_.is_null ()) {
      tsqr_ = Teuchos::rcp (new CholQR<SC, MV, OP> ());
    }
    stepSize_ = stepSize;
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
    const int stepSize = stepSize_;
    int restart = input.resCycle;
    int step = stepSize;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();
    const bool computeRitzValues = input.computeRitzValues;

    // initialize output parameters
    SolverOutput<SC> output {};
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;

    if (outPtr != nullptr) {
      *outPtr << "GmresSstep" << endl;
    }
    Indent indent1 (outPtr);
    if (outPtr != nullptr) {
      *outPtr << "Solver input:" << endl;
      Indent indentInner (outPtr);
      *outPtr << input;
    }

    Teuchos::BLAS<LO, SC> blas;
    Teuchos::LAPACK<LO, SC> lapack;
    dense_matrix_type  H (restart+1, restart, true); // Hessenburg matrix
    dense_matrix_type  T (restart+1, restart, true); // H reduced to upper-triangular matrix
    dense_matrix_type  G (restart+1, step+1, true);  // Upper-triangular matrix from ortho process
    dense_vector_type  y (restart+1, true);
    dense_matrix_type  h (restart+1, 1, true); // used for reorthogonalization
    std::vector<mag_type> cs (restart);
    std::vector<SC> sn (restart);

    mag_type b_norm;  // initial residual norm
    mag_type b0_norm; // initial residual norm, not left preconditioned
    mag_type r_norm;
    mag_type metric;
    vec_type R (B.getMap ());
    vec_type MP (B.getMap ());
    MV  Q (B.getMap (), restart+1);
    vec_type P = * (Q.getVectorNonConst (0));

    // Compute initial residual (making sure R = B - Ax)
    A.apply (X, R);
    R.update (one, B, -one);
    b0_norm = R.norm2 (); // initial residual norm, not preconditioned
    if (input.precoSide == "left") {
      M.apply (R, P);
      r_norm = P.norm2 (); // initial residual norm, left-preconditioned
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
      Tpetra::deep_copy (B, P);
      return output;
    } else if (computeRitzValues) {
      // Invoke standard Gmres for the first restart cycle, to compute
      // Ritz values for use as Newton shifts
      if (outPtr != nullptr) {
        *outPtr << "Run standard GMRES for first restart cycle" << endl;
      }
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
        M.apply (R, P);
        r_norm = P.norm2 (); // residual norm
      }
      else {
        r_norm = output.absResid;
      }
      output.numRests++;
    }

    // Initialize starting vector
    if (input.precoSide != "left") {
      Tpetra::deep_copy (P, R);
    }
    P.scale (one / r_norm);
    y[0] = r_norm;

    // Main loop
    while (output.numIters < input.maxNumIters && ! output.converged) {
      if (outPtr != nullptr) {
        *outPtr << "Restart cycle " << output.numRests << ":" << endl;
      }
      Indent indent2 (outPtr);
      if (outPtr != nullptr) {
        *outPtr << output;
      }

      int iter = 0;
      if (input.maxNumIters < output.numIters+restart) {
        restart = input.maxNumIters-output.numIters;
      }

      // Restart cycle
      for (iter = 0; iter < restart && metric > input.tol; iter+=step) {
        if (outPtr != nullptr) {
          *outPtr << "Current iteration: iter=" << iter
                  << ", restart=" << restart
                  << ", step=" << step
                  << ", metric=" << metric << endl;
        }
        Indent indent3 (outPtr);

        // Compute matrix powers
        for (step=0; step < stepSize && iter+step < restart; step++) {
          if (outPtr != nullptr) {
            *outPtr << "step=" << step
                    << ", stepSize=" << stepSize
                    << ", iter+step=" << (iter+step)
                    << ", restart=" << restart << endl;
          }

          // AP = A*P
          vec_type P  = * (Q.getVectorNonConst (iter+step));
          vec_type AP = * (Q.getVectorNonConst (iter+step+1));
          if (input.precoSide == "none") {
            A.apply (P, AP);
          }
          else if (input.precoSide == "right") {
            M.apply (P, MP);
            A.apply (MP, AP);
          }
          else {
            A.apply (P, MP);
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
        this->projectBelosOrthoManager (iter, step, Q, G);
        const int rank = normalizeCholQR (iter, step, Q, G);
        if (outPtr != nullptr) {
          *outPtr << "Rank of s-step basis: " << rank << endl;
        }
        updateHessenburg (iter, step, output.ritzValues, H, G);

        // Check negative norm
        TEUCHOS_TEST_FOR_EXCEPTION
          (STS::real (H(iter+step, iter+step-1)) < STM::zero (),
           std::runtime_error, "At iteration " << output.numIters << ", H("
           << iter+step << ", " << iter+step-1 << ") = "
           << H(iter+step, iter+step-1) << " < 0.");

        // Convergence check
        if (rank == step+1 && H(iter+step, iter+step-1) != zero) {
          // Copy H to T and apply Givens rotations to new columns of T and y
          for (int iiter = 0; iiter < step; iiter++) {
            for (int i = 0; i <= iter+iiter+1; i++) {
              T(i, iter+iiter) = H(i, iter+iiter);
            }
            this->reduceHessenburgToTriangular(iter+iiter, T, cs, sn, y);
          }
          metric = this->getConvergenceMetric (STS::magnitude (y(iter+step)), b_norm, input);
        }
        else {
          metric = STM::zero ();
        }
      } // End of restart cycle

      // Update solution
      blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
                 Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
                 iter, 1, one,
                 T.values(), T.stride(), y.values(), y.stride());
      Teuchos::Range1D cols(0, iter-1);
      Teuchos::RCP<const MV> Qj = Q.subView(cols);
      if (input.precoSide == "right") {
        dense_vector_type y_iter (Teuchos::View, y.values (), iter);

        MVT::MvTimesMatAddMv (one, *Qj, y_iter, zero, R);
        M.apply (R, MP);
        X.update (one, MP, one);
      }
      else {
        dense_vector_type y_iter (Teuchos::View, y.values (), iter);

        MVT::MvTimesMatAddMv (one, *Qj, y_iter, one, X);
      }
      // Compute real residual (not-preconditioned)
      P = * (Q.getVectorNonConst (0));
      A.apply (X, P);
      P.update (one, B, -one);
      r_norm = P.norm2 (); // residual norm
      output.absResid = r_norm;
      output.relResid = r_norm / b0_norm;

      metric = this->getConvergenceMetric (r_norm, b0_norm, input);
      if (metric <= input.tol) {
        output.converged = true;
      }
      else if (output.numIters < input.maxNumIters) {
        // Initialize starting vector for restart
        if (input.precoSide == "left") { // left-precond'd residual norm
          Tpetra::deep_copy (R, P);
          M.apply (R, P);
          r_norm = P.norm2 ();
        }
        P.scale (one / r_norm);
        y[0] = SC {r_norm};
        for (int i=1; i < restart+1; ++i) {
          y[i] = STS::zero ();
        }
        output.numRests++;
      }
    }

    // Return residual norm as B
    Tpetra::deep_copy (B, P);

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

    // copy: H(j:n-1, j:n-1) = R(j:n-1, j:n-1), i.e., H = R*B
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
    dense_matrix_type r_diag (Teuchos::View, R, s+1, s+1, n, 0);
    dense_matrix_type h_diag (Teuchos::View, H, s+1, s,   n, n);
    Teuchos::BLAS<LO, SC> blas;

    // H = H*R^{-1}
    if (n == 0) { // >> first matrix-power iteration <<
      // diagonal block
      blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                 Teuchos::NON_UNIT_DIAG, s+1, s, one,
                 r_diag.values(), r_diag.stride(),
                 h_diag.values(), h_diag.stride());
    } else  { // >> rest of iterations <<
      for (int j = 1; j < s; j++ ) {
        H(n, n+j) -= H(n, n-1) * R(n-1, j);
      }
      // diagonal block
      blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                 Teuchos::NON_UNIT_DIAG, s, s, one,
                 r_diag.values(), r_diag.stride(),
                 h_diag.values(), h_diag.stride());
      H(n+s, n+s-1) /= R(n+s-1, s-1);

      // upper off-diagonal block: H(0:j-1, j:j+n-2)
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
  normalizeCholQR (const int n,
                   const int s,
                   MV& Q,
                   dense_matrix_type& R)
  {
    // vector to be orthogonalized
    Teuchos::Range1D index_prev(n, n+s);
    MV Qnew = * (Q.subView(index_prev));

    dense_matrix_type r_new (Teuchos::View, R, s+1, s+1, n, 0);

    int rank = 0;
    if (tsqr_ != Teuchos::null) {
      rank = tsqr_->factor (Qnew, r_new);
    }
    else {
      rank = this->normalizeBelosOrthoManager (Qnew, r_new);
    }
    return rank;
  }

private:
  int stepSize_;
  Teuchos::RCP<CholQR<SC, MV, OP> > tsqr_;
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
