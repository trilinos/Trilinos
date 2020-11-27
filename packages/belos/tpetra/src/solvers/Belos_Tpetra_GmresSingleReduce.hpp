#ifndef BELOS_TPETRA_GMRES_SINGLE_REDUCE_HPP
#define BELOS_TPETRA_GMRES_SINGLE_REDUCE_HPP

#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_UpdateNewton.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::Operator<>::scalar_type,
         class MV = Tpetra::MultiVector<SC>,
         class OP = Tpetra::Operator<SC>>
class GmresSingleReduce : public Gmres<SC, MV, OP> {
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

public:
  GmresSingleReduce () :
    base_type::Gmres (),
    stepSize_ (1)
  {
    this->input_.computeRitzValues = true;
  }

  GmresSingleReduce (const Teuchos::RCP<const OP>& A) :
    base_type::Gmres (A),
    stepSize_ (1)
  {
    this->input_.computeRitzValues = true;
  }

  virtual ~GmresSingleReduce () = default;

  virtual void
  setParameters (Teuchos::ParameterList& params) {
    Gmres<SC, MV, OP>::setParameters (params);

    int stepSize = params.get<int> ("Step Size", stepSize_);
    stepSize_ = stepSize;
  }

  void
  setStepSize(int stepSize) {
    stepSize_ = stepSize;
  }

  int
  getStepSize() {
    return stepSize_;
  }

protected:
  virtual void
  setOrthogonalizer (const std::string& ortho)
  {
    if (ortho == "MGS" || ortho == "CGS" || ortho == "CGS2") {
      this->input_.orthoType = ortho;
    } else {
      this->input_.orthoType = "CGS";
      //base_type::setOrthogonalizer (ortho);
    }
  }

private:
  //! Apply the orthogonalization using a single all-reduce
  int
  projectAndNormalizeSingleReduce (int n,
                                   const SolverInput<SC>& input, 
                                   MV& Q,
                                   dense_matrix_type& H,
                                   dense_matrix_type& WORK) const
  {
    Teuchos::RCP< Teuchos::Time > orthTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSingleReduce::LowSynch::Ortho");
    Teuchos::TimeMonitor OrthTimer (*orthTimer);

    int rank = 0;
    if (input.orthoType == "CGS") {
      // default, one-synch CGS, optionally with reortho
      rank = projectAndNormalizeSingleReduce_CGS (n, input, Q, H, WORK);
    } else {
      // one-synch MGS or CGS2, optionally with delayed renorm
      rank = projectAndNormalizeSingleReduce_GS (n, input, Q, H, WORK);
    }
    return rank;
  }


  //! MGS/CGS2 specialization
  int
  projectAndNormalizeSingleReduce_GS (int n,
                                      const SolverInput<SC>& input, 
                                      MV& Q,
                                      dense_matrix_type& H,
                                      dense_matrix_type& T) const
  {
    const SC zero = STS::zero ();
    const SC one  = STS::one  ();
    const SC two  = one + one;
    const mag_type eps  = STS::eps ();
    const mag_type tolOrtho = static_cast<mag_type> (10.0) * eps;

    Teuchos::RCP< Teuchos::Time > dotsTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSingleReduce::LowSynch::dot-prod");
    Teuchos::RCP< Teuchos::Time > projTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSingleReduce::LowSynch::project");

    int rank = 1;
    // ----------------------------------------------------------
    // dot-product for single-reduce orthogonalization
    // ----------------------------------------------------------
    Teuchos::Range1D index_next (n, n+1);
    MV Qnext = * (Q.subView (index_next));

    // vectors to be orthogonalized against
    Teuchos::Range1D index (0, n+1);
    Teuchos::RCP< const MV > Qi = MVT::CloneView (Q, index);

    // compute coefficient, T(:,n:n+1) = Q(:,0:n+1)'*Q(n:n+1)
    Teuchos::RCP< dense_matrix_type > tj
      = Teuchos::rcp (new dense_matrix_type (Teuchos::View, T, n+2, 2, 0, n));
    {
      Teuchos::TimeMonitor LocalTimer (*dotsTimer);
      MVT::MvTransMv (one, *Qi, Qnext, *tj);
    }

    // ----------------------------------------------------------
    // lagged/delayed re-orthogonalization
    // ----------------------------------------------------------
    mag_type prevNorm = STS::real (T(n, n));
    prevNorm = STM::squareroot (prevNorm);
    T(n, n) = SC {prevNorm};
    if (input.needToReortho) {
      Teuchos::Range1D index_old (n, n);
      MV Qold = * (Q.subView (index_old));
      MVT::MvScale (Qold, one / T (n, n)); // normalize

      // update coefficients after reortho
      T (n, n+1) /= T (n, n);
      for (int i = 0; i <= n; i++) {
        T (i, n) /= T (n, n);
      }
      if (n > 0) {
        H (n, n-1) *= T (n, n);
      }
    }

    // ----------------------------------------------------------
    // comopute new coefficients (one-synch MGS/CGS2)
    // ----------------------------------------------------------
    // extract new coefficients 
    for (int i = 0; i <= n+1; i++) {
      H (i, n) = T (i, n+1);
    }

    // update new coefficients 
    Teuchos::BLAS<LO ,SC> blas;
    dense_matrix_type Hnew (Teuchos::View, H, n+1, 1, 0, n);
    if (input.orthoType == "MGS") {
      // H := (I+T)^(-T) H, where T is upper-triangular
      blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
                 Teuchos::TRANS, Teuchos::UNIT_DIAG,
                 n+1, 1,
                 one, T.values(), T.stride(),
                      Hnew.values(), Hnew.stride());
    } else {
      // expand from triangular to full, conjugate
      for (int j = 0; j < n; j++) {
        T(n, j) = STS::conjugate(T(j, n));
      }

      // H := (2*I-L) H, where T is full symmetrix Q'*Q
      dense_matrix_type Told (Teuchos::View, T, n+1, 1, 0, n+1);
      blas.COPY (n+1, Hnew.values(), 1, Told.values(), 1);
      blas.GEMV(Teuchos::NO_TRANS,
                n+1, n+1,
                -one, T.values(), T.stride(),
                      Told.values(), 1,
                 two, Hnew.values(), 1);
    }

    // ----------------------------------------------------------
    // orthogonalize the new vectors against the previous columns
    // ----------------------------------------------------------
    Teuchos::Range1D index_new (n+1, n+1);
    MV Qnew = * (Q.subView (index_new));

    Teuchos::Range1D index_prev (0, n);
    Teuchos::RCP< const MV > Qprev = MVT::CloneView (Q, index_prev);
    {
      Teuchos::TimeMonitor LocalTimer (*projTimer);
      MVT::MvTimesMatAddMv (-one, *Qprev, Hnew, one, Qnew);
    }

    // ----------------------------------------------------------
    // normalize the new vector
    // ----------------------------------------------------------
    // fix the norm
    mag_type oldNorm = STS::real (H(n+1, n));
    for (int i = 0; i <= n; ++i) {
      H(n+1, n) -= (H(i, n)*H(i, n));
    }
    mag_type newNorm = STS::real (H(n+1, n));

    // check
    if (newNorm > oldNorm * tolOrtho) {
      // compute norm
      newNorm = STM::squareroot (newNorm);
      H(n+1, n) = SC {newNorm};

      // scale
      MVT::MvScale (Qnew, one / H (n+1, n)); // normalize
      rank = 1;
    }
    else {
      H(n+1, n) = zero;
      rank = 0;
    }


    return rank;
  }


  //! CGS1 specialization
  int
  projectAndNormalizeSingleReduce_CGS (int n,
                                       const SolverInput<SC>& input, 
                                       MV& Q,
                                       dense_matrix_type& H,
                                       dense_matrix_type& WORK) const
  {
    Teuchos::BLAS<LO, SC> blas;
    const SC one = STS::one ();
    const mag_type eps  = STS::eps ();
    const mag_type tolOrtho = static_cast<mag_type> (10.0) * eps;

    int rank = 0;
    Teuchos::Range1D index_all(0, n+1);
    Teuchos::Range1D index_prev(0, n);
    const MV Qall  = * (Q.subView(index_all));
    const MV Qprev = * (Q.subView(index_prev));

    dense_matrix_type h_all (Teuchos::View, H, n+2, 1, 0, n);
    dense_matrix_type h_prev (Teuchos::View, H, n+1, 1, 0, n);

    // Q(:,0:j+1)'*Q(:,j+1)
    vec_type AP = * (Q.getVectorNonConst (n+1));
    MVT::MvTransMv(one, Qall, AP, h_all);

    // orthogonalize (project)
    MVT::MvTimesMatAddMv (-one, Qprev, h_prev, one, AP);

    // save the norm before ortho
    mag_type oldNorm = STS::real (H(n+1, n));

    // reorthogonalize if requested
    if (input.needToReortho) {
      // Q(:,0:j+1)'*Q(:,j+1)
      dense_matrix_type w_all (Teuchos::View, WORK, n+2, 1, 0, n);
      dense_matrix_type w_prev (Teuchos::View, WORK, n+1, 1, 0, n);
      MVT::MvTransMv(one, Qall, AP, w_all);
      // orthogonalize (project)
      MVT::MvTimesMatAddMv (-one, Qprev, w_prev, one, AP);
      // recompute the norm
      for (int i = 0; i <= n; ++i) {
        w_all(n+1, 0) -= (w_prev(i, 0)*w_prev(i, 0));
      }

      // accumulate results
      blas.AXPY (n+1, one, w_prev.values (), 1, h_prev.values (), 1);
      H(n+1, n) = w_all(n+1, 0); 
    } else {
      for (int i = 0; i <= n; ++i) {
        H(n+1, n) -= (H(i, n)*H(i, n));
      }
    }

    // check for negative norm
    TEUCHOS_TEST_FOR_EXCEPTION
      (STS::real (H(n+1, n)) < STM::zero (), std::runtime_error, "At iteration "
       << n << ", H(" << n+1 << ", " << n << ") = " << H(n+1, n) << " < 0.");
    // Check for zero norm.  OK to take real part of H(n+1, n), since
    // this is a magnitude (cosine) anyway and therefore should always
    // be real.
    const mag_type H_np1_n = STS::real (H(n+1, n));
    if (H_np1_n > oldNorm * tolOrtho) {
      const mag_type H_np1_n_sqrt = STM::squareroot (H_np1_n);
      H(n+1, n) = SC {H_np1_n_sqrt};
      MVT::MvScale (AP, STS::one () / H(n+1, n)); // normalize
      rank = 1;
    }
    else {
      H(n+1, n) = STS::zero ();
      rank = 0;
    }
    return rank;
  }


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
    const SC zero = STS::zero ();
    const SC one  = STS::one ();
    const bool computeRitzValues = input.computeRitzValues;


    // timers
    Teuchos::RCP< Teuchos::Time > spmvTimer = Teuchos::TimeMonitor::getNewCounter ("GmresSingleReduce::Sparse Mat-Vec");

    SolverOutput<SC> output {};
    // initialize output parameters
    output.numRests = 0;
    output.numIters = 0;
    output.converged = false;

    mag_type b_norm;  // initial residual norm
    mag_type b0_norm; // initial residual norm, not left-preconditioned
    mag_type r_norm;
    mag_type r_norm_imp;

    bool zeroOut = false; 
    MV Q (B.getMap (), restart+1, zeroOut);
    vec_type R  (B.getMap (), zeroOut);
    vec_type Y  (B.getMap (), zeroOut);
    vec_type MP (B.getMap (), zeroOut);
    vec_type P = * (Q.getVectorNonConst (0));

    Teuchos::BLAS<LO ,SC> blas;
    dense_matrix_type H (restart+1, restart,   true);
    dense_matrix_type T (restart+1, restart+2, true);
    dense_vector_type y (restart+1, true);
    dense_vector_type z (restart+1, true);
    std::vector<mag_type> cs (restart);
    std::vector<SC> sn (restart);
    
    #ifdef HAVE_TPETRA_DEBUG
    dense_matrix_type H2 (restart+1, restart,   true);
    dense_matrix_type H3 (restart+1, restart,   true);
    #endif

    // initial residual (making sure R = B - Ax)
    {
      Teuchos::TimeMonitor LocalTimer (*spmvTimer);
      A.apply (X, R);
    }
    R.update (one, B, -one);
    b0_norm = R.norm2 (); // initial residual norm, no-preconditioned
    if (input.precoSide == "left") {
      M.apply (R, P);
      r_norm = P.norm2 (); // initial residual norm, left-preconditioned
    } else {
      r_norm = b0_norm;
    }
    b_norm = r_norm;

    mag_type metric = this->getConvergenceMetric (b0_norm, b0_norm, input);
    if (metric <= input.tol) {
      if (outPtr != NULL) {
        *outPtr << "Initial guess' residual norm " << b0_norm
                << " meets tolerance " << input.tol << endl;
      }
      output.absResid = b0_norm;
      output.relResid = STM::one ();
      output.converged = true;
      return output;
    } else if (outPtr != NULL) {
      *outPtr << "Initial guess' residual norm " << b0_norm << endl;
    }

    // compute Ritz values as Newton shifts
    if (computeRitzValues) {
      // Invoke ordinary GMRES for the first restart
      SolverInput<SC> input_gmres = input;
      input_gmres.maxNumIters = input.resCycle;
      input_gmres.computeRitzValues = computeRitzValues;
      if (input.orthoType == "MGS") {
        // MGS1 or MGS2
        input_gmres.orthoType = "IMGS";
        if (input.needToReortho) {
          input_gmres.maxOrthoSteps = 2;
        } else {
          input_gmres.maxOrthoSteps = 1;
        }
      } else if (input.orthoType == "CGS") {
        // CGS1 or CGS2
        input_gmres.orthoType = "ICGS";
        if (input.needToReortho) {
          input_gmres.maxOrthoSteps = 2;
        } else {
          input_gmres.maxOrthoSteps = 1;
        }
      } else if (input.orthoType == "CGS2") {
        // CGS2
        input_gmres.orthoType = "ICGS";
        input_gmres.maxOrthoSteps = 2;
      }
      if (outPtr != nullptr) {
        *outPtr << "Run standard GMRES for first restart cycle" << endl;
      }

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

    // initialize starting vector
    if (input.precoSide != "left") {
      Tpetra::deep_copy (P, R);
    }
    P.scale (one / r_norm);
    y[0] = SC {r_norm};
    const int s = getStepSize ();
    // main loop
    bool delayed_ortho = ((input.orthoType == "MGS" && input.needToReortho) ||
                          (input.orthoType == "CGS2"));
    int iter = 0;
    while (output.numIters < input.maxNumIters && ! output.converged) {
      if (input.maxNumIters < output.numIters+restart) {
        restart = input.maxNumIters-output.numIters;
      }
      // restart cycle
      for (; iter < restart && (metric > input.tol && !STS::isnaninf (metric)); ++iter) {
        // AP = A*P
        vec_type P  = * (Q.getVectorNonConst (iter));
        vec_type AP = * (Q.getVectorNonConst (iter+1));
        if (input.precoSide == "none") { // no preconditioner
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
        if (computeRitzValues) {
          //AP.update (-output.ritzValues[iter],  P, one);
          const complex_type theta = output.ritzValues[iter % s];
          UpdateNewton<SC, MV>::updateNewtonV (iter, Q, theta);
        }
        output.numIters++; 

        // Orthogonalization
        projectAndNormalizeSingleReduce (iter, input, Q, H, T);

        // Convergence check
        if (!delayed_ortho || iter > 0) {
          int check = (delayed_ortho ? iter-1 : iter);
          if (outPtr != nullptr) {
            *outPtr << "Current iteration: iter=" << iter
                    << ", restart=" << restart
                    << ", metric=" << metric << endl;
            Indent indent3 (outPtr);
          }

          // Shift back for Newton basis
          if (computeRitzValues) {
            const complex_type theta = output.ritzValues[check % s];
            UpdateNewton<SC, MV>::updateNewtonH (check, H, theta);
          }
          #ifdef HAVE_TPETRA_DEBUG
          this->checkNumerics (outPtr, iter, check, A, M, Q, X, B, y,
                               H, H2, H3, cs, sn, input);
          #endif

          if (H(check+1, check) != zero) {
            // Apply Givens rotations to new column of H and y
            this->reduceHessenburgToTriangular (check, H, cs, sn, y.values ());
            // Convergence check
            metric = this->getConvergenceMetric (STS::magnitude (y[check+1]),
                                                 b_norm, input);
          }
          else {
            H(check+1, check) = zero;
            metric = STM::zero ();
          }
        }
      } // end of restart cycle 

      if (delayed_ortho) {
        int check = iter-1;
        // Shift back for Newton basis
        if (computeRitzValues) {
          const complex_type theta = output.ritzValues[check % s];
          UpdateNewton<SC, MV>::updateNewtonH (check, H, theta);
        }
        #ifdef HAVE_TPETRA_DEBUG
        this->checkNumerics (outPtr, iter, check, A, M, Q, X, B, y,
                             H, H2, H3, cs, sn, input);
        #endif

        if (H(check+1, check) != zero) {
          // Apply Givens rotations to new column of H and y
          this->reduceHessenburgToTriangular (check, H, cs, sn, y.values ());
          // Convergence check
          metric = this->getConvergenceMetric (STS::magnitude (y[check+1]),
                                               b_norm, input);
        }
        else {
          H(check+1, check) = zero;
          metric = STM::zero ();
        }
      }

      if (iter < restart) {
        // save the old solution, just in case explicit residual norm failed the convergence test
        Tpetra::deep_copy (Y, X);
        blas.COPY (1+iter, y.values(), 1, z.values(), 1);
      }
      r_norm_imp = STS::magnitude (y (iter)); // save implicit residual norm
      if (iter > 0) {
        // Update solution
        blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                   Teuchos::NON_UNIT_DIAG, iter, 1, one,
                   H.values(), H.stride(), y.values (), y.stride ());
        Teuchos::Range1D cols(0, iter-1);
        Teuchos::RCP<const MV> Qj = Q.subView(cols);
        y.resize (iter);
        if (input.precoSide == "right") {
          MVT::MvTimesMatAddMv (one, *Qj, y, zero, R);
          M.apply (R, MP);
          X.update (one, MP, one);
        }
        else {
          MVT::MvTimesMatAddMv (one, *Qj, y, one, X);
        }
        y.resize (restart+1);
      }

      // Compute explicit residual vector in preparation for restart.
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
        output.converged = true;
        return output;
      }
      else if (output.numIters < input.maxNumIters) {
        // Restart, only if max inner-iteration was reached.
        // Otherwise continue the inner-iteration.
        if (iter >= restart || H(iter,iter-1) == zero) { // done with restart cycyle, or probably lost ortho
          // Initialize starting vector for restart
          iter = 0;
          P = * (Q.getVectorNonConst (0));
          if (input.precoSide == "left") {
            M.apply (R, P);
            // FIXME (mfh 14 Aug 2018) Didn't we already compute this above?
            r_norm = P.norm2 ();
          }
          else {
            // set the starting vector
            Tpetra::deep_copy (P, R);
          }
          P.scale (one / r_norm);
          y[0] = SC {r_norm};
          for (int i=1; i < restart+1; i++) {
            y[i] = zero;
          }
          output.numRests++; // restart
        }
        else {
          // reset to the old solution
          if (outPtr != nullptr) {
            *outPtr << " > not-restart with iter=" << iter << endl;
          }
          Tpetra::deep_copy (X, Y);
          blas.COPY (1+iter, z.values(), 1, y.values(), 1);
        }
      }
    }

    return output;
  }
  
  int stepSize_; // "step size" for Newton basis
};

template<class SC, class MV, class OP,
         template<class, class, class> class KrylovSubclassType>
class SolverManager;

// This is the Belos::SolverManager subclass that gets registered with
// Belos::SolverFactory.
template<class SC, class MV, class OP>
using GmresSingleReduceSolverManager = SolverManager<SC, MV, OP, GmresSingleReduce>;

/// \brief Register GmresSingleReduceSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_GmresSingleReduce (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_GMRES_SINGLE_REDUCE_HPP
