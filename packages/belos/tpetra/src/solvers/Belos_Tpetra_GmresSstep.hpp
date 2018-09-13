#ifndef BELOS_TPETRA_GMRES_SSTEP_HPP
#define BELOS_TPETRA_GMRES_SSTEP_HPP

#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_UpdateNewton.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::Operator<>::scalar_type,
	 class MV = Tpetra::MultiVector<SC>,
	 class OP = Tpetra::Operator<SC>>
class GmresSstep : public Gmres<SC, MV, OP>  {
private:
  using base_type = Gmres<SC, MV, OP>;
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using mag_type = real_type;
  using complex_type = std::complex<real_type>;
  using MVT = Belos::MultiVecTraits<SC, MV>;
  using LO = typename MV::local_ordinal_type;
  using STM = Teuchos::ScalarTraits<mag_type>;
  using vec_type = Tpetra::Vector<SC, LO,
				  typename MV::global_ordinal_type,
				  typename MV::node_type>;  
  using dense_matrix_type = Teuchos::SerialDenseMatrix<LO, SC>;
  using dense_vector_type = Teuchos::SerialDenseVector<LO, SC>;  

public:
  GmresSstep () :
    base_type::Gmres (),
    stepSize_ (1)
  {}

  GmresSstep (const Teuchos::RCP<const OP>& A) :
    base_type::Gmres (A),
    stepSize_ (1)
  {}

  virtual ~GmresSstep ()
  {}

  virtual void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const
  {
    Gmres<SC, MV, OP>::getParameters (params, defaultValues);

    const int stepSize = defaultValues ? 100 : stepSize_;

    params.set ("Step Size", stepSize );
  }

  virtual void
  setParameters (Teuchos::ParameterList& params) {
    Gmres<SC, MV, OP>::setParameters (params);

    int stepSize = stepSize_;
    if (params.isParameter ("Step Size")) {
      stepSize = params.get<int> ("Step Size");
    }
    stepSize_ = stepSize;
  }

private:
  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X,
               vec_type& B,
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input)
  {
    using std::endl;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();
    const bool computeRitzValues = true;    
    int restart = input.resCycle;
    int step = input.stepSize;
    SolverOutput<SC> output {};
    
    if (outPtr != nullptr) {
      *outPtr << "s-step GMRES" << endl;
    }
    Indent indent1 (outPtr);    
    if (outPtr != nullptr) {
      *outPtr << input << endl;
    }
    
    mag_type b_norm;  // initial residual norm
    mag_type b0_norm; // initial residual norm, not left-preconditioned
    mag_type r_norm;
    MV Q (B.getMap (), restart+1);
    vec_type P = * (Q.getVectorNonConst (0));    
    vec_type R (B.getMap ());
    vec_type MP (B.getMap ());

    if (outPtr != nullptr) {
      *outPtr << "Compute initial residual R = B - Ax" << endl;
    }
    if (input.precoSide == "right") {
      M.apply (X, MP);
      A.apply (MP, R);
    }
    else {
      A.apply (X, R);
    }
    R.update (one, B, -one);
    b0_norm = R.norm2 (); // residual norm, not preconditioned
    if (input.precoSide == "left") {
      M.apply (R, P);
      b_norm = P.norm2 (); // residual norm, left preconditioned
    }
    else {
      Tpetra::deep_copy (P, R);
      b_norm = b0_norm;
    }
    mag_type metric = this->getConvergenceMetric (b0_norm, b0_norm, input);
    if (outPtr != nullptr) {
      *outPtr << "Absolute residual:  " << b0_norm << endl
	      << "Relative residual:  " << STM::one () << endl
	      << "Convergence metric: " << metric << endl
	      << "Tolerance:          " << input.tol << endl;
    }

    if (metric <= input.tol) {
      if (outPtr != nullptr) {
        *outPtr << "Initial solution meets tolerance" << endl;
      }
      output.absResid = b_norm;
      output.relResid = STM::one ();
      output.numIters = 0;
      output.converged = true;
      Tpetra::deep_copy (B, P); // return residual norm as B
      return output;
    }    

    Teuchos::BLAS<LO ,SC> blas;
    Teuchos::LAPACK<LO ,SC> lapack;
    dense_matrix_type  H (restart+1, restart, true); // Upper Hessenburg matrix
    dense_matrix_type  T (restart+1, restart, true); // H reduced to upper tri
    dense_matrix_type  G (restart+1, step+1, true); // Upper-tri matrix from ortho
    dense_vector_type  y (restart+1, true);
    dense_matrix_type  h (restart+1, 1, true);
    std::vector<real_type> cs (restart);
    std::vector<SC> sn (restart);

    // initialize output parameters
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;

    if (computeRitzValues) {
      if (outPtr != nullptr) {
	*outPtr << "First restart cycle: Standard GMRES (for Ritz values)"
		<< endl;
      }
      SolverInput<SC> input_gmres = input;
      input_gmres.maxNumIters = input.resCycle;
      input_gmres.computeRitzValues = true;
      output = Gmres<SC, MV, OP>::solveOneVec (outPtr, X, R, A, M,
					       input_gmres);
      if (outPtr != nullptr) {
	*outPtr << "Standard GMRES results:" << endl;
	Indent indent2 (outPtr);
	*outPtr << output;
      }
      if (output.converged) {
	if (outPtr != nullptr) {
          *outPtr << "Converged to tolerance " << input.tol << " in "
		  << output.numIters << " standard GMRES iterations"
		  << endl;
        }
	Tpetra::deep_copy (B, R); // return residual vector as B
	return output; // ordinary GMRES converged
      }
      
      if (input.precoSide == "left") {
        M.apply (R, P);
        r_norm = P.norm2 (); // residual norm
      }
      else {
	Tpetra::deep_copy (P, R);
        r_norm = output.absResid;
      }
      output.numRests++;
    }

    // initialize starting vector
    P.scale (one / b_norm);
    y[0] = SC {b_norm};
    // main loop
    while (output.numIters < input.maxNumIters && ! output.converged) {
      int iter = 0;
      if (input.maxNumIters < output.numIters+restart) {
        restart = input.maxNumIters-output.numIters;
      }
      if (outPtr != nullptr) {
	*outPtr << "Main loop: numIters=" << output.numIters
		<< ", restart=" << restart << endl;
      }
      // restart cycle
      for (iter = 0; iter < restart && metric > input.tol; iter += step) {
	Indent indent2 (outPtr);
	if (outPtr != nullptr) {
	  *outPtr << "Restart cycle: iter=" << iter << ", step=" << step
		  << endl;
	}
	if (outPtr != nullptr) {
	  *outPtr << "Matrix powers: stepSize=" << input.stepSize << endl;
	}
        // compute matrix powers
        for (step = 0; step < input.stepSize && iter+step < restart; ++step) {
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
          if (computeRitzValues) {
            //AP.update (-output.ritzValues(step), P, one);
            const complex_type theta = output.ritzValues[step];
            UpdateNewton<SC, MV>::updateNewtonV (iter+step, Q, theta);
          }
          output.numIters++; 
        }

        this->projectBelosOrthoManager (iter, step, Q, G);
        int rank = this->normalize (iter, step, Q, G);
        this->updateHessenburg (iter, step, output.ritzValues, H, G);

        // Check for negative norm
        TEUCHOS_TEST_FOR_EXCEPTION
          (STS::real (H(iter+step, iter+step-1)) < STM::zero (),
	   std::runtime_error, "At iteration " << output.numIters
	   << ", H(" << iter+step << ", " << iter+step-1 << ") = "
	   << H(iter+step, iter+step-1) << " < 0.");

        if (rank == step+1 && H(iter+step, iter+step-1) != zero) {
          // Copy H to T and apply Givens rotations to new columns of T and y
          for (int iiter=0; iiter<step; iiter++) {
            for (int i=0; i<=iter+iiter+1; i++) {
	      T(i, iter+iiter) = H(i, iter+iiter);
	    }
            this->reduceHessenburgToTriangular (iter+iiter, T, cs, sn, y.values());
          }
          metric = this->getConvergenceMetric (STS::magnitude (y[iter+step]),
					       b_norm, input);
        }
	else {
          metric = STM::zero ();
        }
	if (outPtr != nullptr) {
	  *outPtr << "Convergence metric: " << metric << endl;
	}
      } // end of restart cycle

      if (outPtr != nullptr) {
	*outPtr << "Update solution after restart cycle" << endl;
      }
      blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                 Teuchos::NON_UNIT_DIAG, iter, 1, one,
                 T.values(), T.stride(), y.values(), y.stride());
      Teuchos::Range1D cols(0, iter-1);
      Teuchos::RCP<const MV> Qj = Q.subView(cols);
      if (input.precoSide == "right") {
	// FIXME (mfh 16 Sep 2018) not quite sure if this is needed.
	dense_vector_type y_restart (Teuchos::View, y.values (), iter);
        MVT::MvTimesMatAddMv (one, *Qj, y_restart, zero, R);
        M.apply (R, MP);
        X.update (one, MP, one);
      }
      else {
	// FIXME (mfh 16 Sep 2018) not quite sure if this is needed.
	dense_vector_type y_restart (Teuchos::View, y.values (), iter);
        MVT::MvTimesMatAddMv (one, *Qj, y_restart, one, X);
      }

      if (outPtr != nullptr) {
	*outPtr << "Compute explicit unpreconditioned residual" << endl;
      }
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
        if (input.precoSide == "left") {
	  Tpetra::deep_copy (R, P);
          M.apply (R, P);
          r_norm = P.norm2 ();
        }
        P.scale (one / r_norm);
        y[0] = SC {r_norm};
        for (int i=1; i < restart+1; i++) {
	  y[i] = zero;
	}
        output.numRests++;
      }
      if (outPtr != nullptr) {
	*outPtr << "End of restart loop:" << endl;
	Indent lastIndent (outPtr);
	*outPtr << "Converged:          " << (output.converged ? "true" : "false")
		<< endl
		<< "Absolute residual:  " << output.absResid << endl
		<< "Relative residual:  " << output.relResid << endl
		<< "Convergence metric: " << metric << endl
		<< "Tolerance:          " << input.tol << endl;
      }
    }

    if (outPtr != nullptr) {
      *outPtr << "End of solve:" << endl;
      Indent lastIndent (outPtr);
      *outPtr << "Converged:          " << (output.converged ? "true" : "false")
	      << endl
	      << "Absolute residual:  " << output.absResid << endl
	      << "Relative residual:  " << output.relResid << endl
	      << "Convergence metric: " << metric << endl
	      << "Tolerance:          " << input.tol << endl;
    }
    // return residual norm as B
    Tpetra::deep_copy (B, P);
    return output;
  }

private:
  void
  updateHessenburg (const int n,
		    const int s,
		    std::vector<std::complex<real_type> >& S,
		    dense_matrix_type& H,
		    dense_matrix_type& R) const
  {
    const SC one  = STS::one ();
    const SC zero = STS::zero ();

    // copy: H(j:n-1, j:n-1) = R(j:n-1, j:n-1), i.e., H = R*B
    for (int j = 0; j < s; ++j) {
      for (int i = 0; i <= n+j+1; ++i) {
        H(i, n+j) = R(i, j+1);
        if (int (S.size()) > j) {
          //H(i, n+j) += S[j].real * R(i, j);
          H(i, n+j) += UpdateNewton<SC, MV>::updateNewtonH (i, j, R, S[j]);
        }
      }
      for (int i=n+j+2; i<=n+s; ++i) {
	H(i, n+j) = zero;
      }
    }

    dense_matrix_type r_diag (Teuchos::View, R, s+1, s+1, n, 0);
    dense_matrix_type h_diag (Teuchos::View, H, s+1, s,   n, n);
    Teuchos::BLAS<LO ,SC> blas;

    // H = H*R^{-1}
    if (n == 0) { // >> first matrix-power iteration <<
      // diagonal block
      blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                 Teuchos::NON_UNIT_DIAG, s+1, s, one,
                 r_diag.values(), r_diag.stride(),
                 h_diag.values(), h_diag.stride());
    }
    else { // >> rest of iterations <<
      for (int j = 1; j < s; ++j) {
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

  int
  normalize (int n, int s, MV &Q, dense_matrix_type &R)
  {
    // vector to be orthogonalized
    Teuchos::Range1D index_prev (n, n+s);
    auto Qnew = Q.subViewNonConst (index_prev);
    dense_matrix_type r_new (Teuchos::View, R, s+1, s+1, n, 0);
    return this->normalizeBelosOrthoManager (*Qnew, r_new);
  }

private:
  int stepSize_;
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
