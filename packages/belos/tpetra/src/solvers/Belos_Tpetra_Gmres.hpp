#ifndef BELOS_TPETRA_GMRES_HPP
#define BELOS_TPETRA_GMRES_HPP

#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "Belos_Tpetra_Krylov.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace BelosTpetra {
namespace Impl {

namespace { // (anonymous)

// LAPACK's _GSEQR routine returns eigenvalues as (array of real
// parts, array of imaginary parts) for real Scalar types, but returns
// eigenvalues as (array of pairs (real, imag)) for complex Scalar
// types.  computeRitzValues() provides a single interface for both
// cases.
template<class LO, class SC, bool isComplex = Teuchos::ScalarTraits<SC>::isComplex>
struct ComputeRitzValues {
  static void
  run (const int iter,
       Teuchos::SerialDenseMatrix<LO, SC>& G,
       std::vector<SC>& ritzValues);
};

// Real Scalar (SC) case.
template<class LO, class SC>
struct ComputeRitzValues<LO, SC, false> {
  static void
  run (const int iter,
       Teuchos::SerialDenseMatrix<LO, SC>& G,
       std::vector<SC>& ritzValues)
  {
    using STS = Teuchos::ScalarTraits<SC>;
    using mag_type = typename STS::magnitudeType;
    
    if (ritzValues.size () < static_cast<std::size_t> (iter)) {
      ritzValues.resize (iter);
    }

    // FIXME (mfh 14 Aug 2018) The code as I found it ignored the
    // imaginary Ritz values.  I'm just going to leave that in place
    // for now.
    
    //std::vector<mag_type> WR (iter);
    std::vector<mag_type>& WR = ritzValues;
    std::vector<mag_type> WI (iter);
    Teuchos::LAPACK<LO, SC> lapack;
    SC TEMP = STS::zero ();    
    LO info = 0;
    LO lwork = -1;
    lapack.HSEQR ('E', 'N', iter, 1, iter, G.values (), G.stride (), 
		  WR.data (), WI.data (), nullptr, 1, &TEMP, lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION
      (info != 0, std::runtime_error, "LAPACK {D,S}SEQR LWORK query failed "
       "with INFO = " << info << " != 0.");
    lwork = static_cast<LO> (TEMP);
    TEUCHOS_TEST_FOR_EXCEPTION
      (lwork < LO {0}, std::runtime_error, "LAPACK {C,Z}SEQR LWORK query "
       "returned LWORK = " << lwork << " < 0.");
    
    std::vector<SC> WORK (lwork);
    lapack.HSEQR ('E', 'N', iter, 1, iter, G.values (), G.stride (), 
		  WR.data (), WI.data (), nullptr, 1, WORK.data (),
		  lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION
      (info != 0, std::runtime_error, "LAPACK {D,S}SEQR failed "
       "with INFO = " << info << " != 0.");
    // FIXME (mfh 14 Aug 2018) The code as I found it ignored the
    // imaginary Ritz values.  I'm just going to leave that in place
    // for now.
  }
};

// Complex Scalar (SC) case.  
template<class LO, class SC>
struct ComputeRitzValues<LO, SC, true> {
  static void
  run (const int iter,
       Teuchos::SerialDenseMatrix<LO, SC>& G,
       std::vector<SC>& ritzValues)
  {
    using STS = Teuchos::ScalarTraits<SC>;

    if (ritzValues.size () < static_cast<std::size_t> (iter)) {
      ritzValues.resize (iter);
    }

    Teuchos::LAPACK<LO, SC> lapack;
    SC TEMP = STS::zero ();
    LO lwork = -1;
    LO info = 0;
    lapack.HSEQR ('E', 'N', iter, 1, iter, G.values (), G.stride (), 
		  ritzValues.data (), nullptr, 1, &TEMP, lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION
      (info != 0, std::runtime_error, "LAPACK {C,Z}SEQR LWORK query failed "
       "with INFO = " << info << " != 0.");
    lwork = static_cast<LO> (STS::real (TEMP));
    TEUCHOS_TEST_FOR_EXCEPTION
      (lwork < LO {0}, std::runtime_error, "LAPACK {C,Z}SEQR LWORK query "
       "returned LWORK = " << lwork << " < 0.");
    
    std::vector<SC> WORK (lwork);
    lapack.HSEQR ('E', 'N', iter, 1, iter, G.values (), G.stride (), 
		  ritzValues.data (), nullptr, 1, WORK.data (), lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION
      (info != 0, std::runtime_error, "LAPACK {C,Z}SEQR failed "
       "with INFO = " << info << " != 0.");
  }
};

template<class LO, class SC>
void
computeRitzValues (const int iter,
		   Teuchos::SerialDenseMatrix<LO, SC>& G,
		   std::vector<SC>& ritzValues)
{
  ComputeRitzValues<LO, SC>::run (iter, G, ritzValues);
}
  
} // namespace (anonymous)

template<class SC = Tpetra::Operator<>::scalar_type,
	 class MV = Tpetra::MultiVector<SC>,
	 class OP = Tpetra::Operator<SC>>
class Gmres : public Krylov<SC, MV, OP> {
public:
  using base_type = Krylov<SC, MV, OP>;

protected:
  using vec_type = typename base_type::vec_type;

private:
  using MVT = Belos::MultiVecTraits<SC, MV>;
  using LO = typename MV::local_ordinal_type;
  using STS = Teuchos::ScalarTraits<SC>;
  using mag_type = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<mag_type>;
  using device_type = typename MV::device_type;

  using ortho_type = Belos::OrthoManager<SC, MV>;
  using dense_matrix_type = Teuchos::SerialDenseMatrix<LO, SC>;
  using dense_vector_type = Teuchos::SerialDenseVector<LO, SC>;

public:
  Gmres () :
    Krylov<SC, MV, OP>::Krylov ()
  {}

  Gmres (const Teuchos::RCP<const OP>& A) :
    Krylov<SC, MV, OP>::Krylov (A)
  {}

  virtual ~Gmres()
  {}

  virtual void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const override
  {
    Krylov<SC, MV, OP>::getParameters (params, defaultValues);

    const int resCycle = defaultValues ? 30 : this->input_.resCycle;
    params.set ("Num Blocks", resCycle);
    params.set ("Maximum Restarts", this->input_.maxNumIters);
  }

  virtual void 
  setParameters (Teuchos::ParameterList& params) override
  {    
    Krylov<SC, MV, OP>::setParameters (params);

    bool computeRitzValues = this->input_.computeRitzValues;
    if (params.isParameter ("Compute Ritz Values")) {
      this->input_.computeRitzValues = params.get<bool> ("Compute Ritz Values");
    }

    bool needToReortho = this->input_.needToReortho;
    if (params.isParameter ("Reorthogonalize Blocks")) {
      this->input_.needToReortho = params.get<bool> ("Reorthogonalize Blocks");
    }

    int resCycle = this->input_.resCycle;
    if (params.isParameter ("Num Blocks")) {
      const int resCycle = params.get<int> ("Num Blocks");
      TEUCHOS_TEST_FOR_EXCEPTION
        (resCycle < 0, std::invalid_argument,
         "\"Num Blocks\" (restart length) = " << resCycle << " < 0.");
      this->input_.resCycle = resCycle;
    }

    if (params.isParameter ("Orthogonalization")) {
      this->input_.orthoType = params.get<std::string> ("Orthogonalization");
    }
  }

private:
  //! Create Belos::OrthoManager instance.
  void
  setOrthogonalizer (const std::string& ortho)
  {
    if (ortho_.get () == nullptr || this->input_.orthoType != ortho) {
      // Since setOrthogonalizer only gets called on demand, we know
      // the preconditioner (if any) at this point.  Thus, we can use
      // Belos::OrthoManagerFactory here.
      Belos::OrthoManagerFactory<SC, MV, OP> factory;
      Teuchos::RCP<const OP> M = this->getPreconditioner ();
      Teuchos::RCP<Belos::OutputManager<SC>> outMan; // can be null
      Teuchos::RCP<Teuchos::ParameterList> params; // can be null
      ortho_ = factory.makeMatOrthoManager (ortho, M, outMan, "Belos", params);
      TEUCHOS_TEST_FOR_EXCEPTION
	(ortho_.get () == nullptr, std::runtime_error, "Gmres: Failed to "
	 "create (Mat)OrthoManager of type \"" << ortho << "\".");
      this->input_.orthoType = ortho;
    }
  }

  int
  projectAndNormalize (const int n,
		       const SolverInput<SC>& /* input */,
		       MV& Q,
		       dense_matrix_type& H,
		       dense_matrix_type& /* WORK */)
  {
    return this->projectAndNormalizeBelosOrthoManager (n, Q, H);
  }

  // ! Apply the orthogonalization using Belos' OrthoManager
  int
  projectAndNormalizeBelosOrthoManager (int n, MV &Q, dense_matrix_type &H)
  {
    using Teuchos::RCP;    
    using Teuchos::rcp;
    
    // vector to be orthogonalized
    vec_type AP = * (Q.getVectorNonConst (n+1));
    // vectors to be orthogonalized against
    Teuchos::Range1D index_prev(0, n);
    RCP<const MV> Q_prev = MVT::CloneView (Q, index_prev);
    Teuchos::Array<RCP<const MV>> Q_array( 1, Q_prev );

    // orthogonalize scalars
    auto h_new = rcp (new dense_matrix_type (Teuchos::View, H, n+1, 1, 0, n));
    Teuchos::Array<RCP<dense_matrix_type>> h_array (1, h_new);
    auto r_new = rcp (new dense_matrix_type (Teuchos::View, H, 1, 1, n+1, n));

    if (ortho_.get () == nullptr) {
      setOrthogonalizer (this->input_.orthoType);
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (ortho_.get () == nullptr, std::logic_error, "Gmres: Failed to create "
       "(Mat)OrthoManager.  This should never happen.  Please report this bug "
       "to the Belos developers.");
    // FIXME (mfh 14 Aug 2018) For left preconditioning, we may need
    // to use the MatOrthoManager interface.
    return ortho_->projectAndNormalize (AP, h_array, r_new, Q_array);
  }

  //! Reduce a column of Henssenburg matrix to triangular form
  void
  reduceHessenburgToTriangular(const int j,
			       dense_matrix_type& H,
			       std::vector<mag_type>& cs,
			       std::vector<SC>& sn,
			       dense_vector_type& y) const
  {
    Teuchos::BLAS<LO, SC> blas;
    // Apply previous Givens rotations to new column of H
    for (int i = 0; i < j; ++i) {
      blas.ROT (1, &H(i, j), 1, &H(i+1, j), 1, &cs[i], &sn[i]);
    }
    // Calculate new Givens rotation
    blas.ROTG (&H(j, j), &H(j+1, j), &cs[j], &sn[j]);
    H(j+1, j) = STS::zero ();
    // Update RHS w/ new transformation
    blas.ROT (1, &y(j), 1, &y(j+1), 1, &cs[j], &sn[j]);
  }

  //! Sort Ritz values, using the Leja ordering.
  void
  sortRitzValues (const LO m, std::vector<SC>& RR)
  {
    LO index = -1;
    SC value = STS::zero();
    for (int i = 0; i < m; i++) {
      LO next_index = index;
      SC next_value = value;
      for (int k = i;  k < m; k++ ) {
        if( std::abs( next_value - value ) <= std::abs( RR[k] - value ) ) {
          next_index = k;
          next_value = RR[k];
        }
      }
      RR[next_index] = RR[i];
      RR[i] = next_value;
      index = next_index;
      value = next_value;
    }
  }

protected:  
  SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in X/out X
               vec_type& B, // in B/out R (not left-preconditioned)
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input) override
  {
    using std::endl;
    int restart = input.resCycle;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();

    SolverOutput<SC> output {};

    mag_type b_norm;  // initial residual norm
    mag_type b0_norm; // initial residual norm (before left-preconditioner)
    mag_type r_norm;
    vec_type R (B.getMap ());
    vec_type MP (B.getMap ());
    MV Q (B.getMap (), restart+1);
    vec_type P = * (Q.getVectorNonConst (0));

    // initial residual (making sure R = B - Ax)
    if (input.precoSide == "right") {
      // right-preconditioner
      M.apply (X, MP);
      A.apply (MP, R);
    }
    else {
      A.apply (X, R);
    }
    R.update (one, B, -one);
    b0_norm = R.norm2 (); // residual norm, not-preconditioned
    if (input.precoSide == "left") {
      M.apply (R, P);
      b_norm = P.norm2 (); // residual norm, left-preconditioned
    }
    else {
      Tpetra::deep_copy (P, R);
      b_norm = b0_norm;
    }
    mag_type metric = this->getConvergenceMetric (b0_norm, b0_norm, input);

    if (metric <= input.tol) {
      if (outPtr != NULL) {
        *outPtr << "Initial guess' residual norm " << b_norm
                << " meets tolerance " << input.tol << endl;
      }
      output.absResid = b_norm;
      output.relResid = STM::one ();
      output.numIters = 0;
      output.converged = true;
      // return residual norm as B
      Tpetra::deep_copy (B, P);
      return output;
    }
    
    Teuchos::BLAS<LO ,SC> blas;
    Teuchos::LAPACK<LO ,SC> lapack;
    dense_matrix_type  H (restart+1, restart, true);
    dense_matrix_type  G (restart+1, restart, true); // only for Ritz values
    dense_vector_type  y (restart+1, true);
    dense_matrix_type  h (restart+1, 1, true); // for reorthogonalization
    std::vector<mag_type> cs (restart);
    std::vector<SC> sn (restart);

    // initialize output parameters
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;

    // initialize starting vector
    P.scale (one / b_norm);
    y[0] = SC {b_norm};
    
    // main loop
    while (output.numIters < input.maxNumIters && ! output.converged) {
      int iter = 0;
      if (input.maxNumIters < output.numIters+restart) {
        restart = input.maxNumIters-output.numIters;
      }
      // restart cycle
      for (iter = 0; iter < restart && metric > input.tol; ++iter) {
        // AP = A*P
        vec_type P  = * (Q.getVectorNonConst (iter));
        vec_type AP = * (Q.getVectorNonConst (iter+1));
        if (input.precoSide == "none") {
          A.apply (P, AP);
        }
	else if (input.precoSide == "right") {
          M.apply (P, MP);
          A.apply (MP, AP);
        }
	else { // left
          A.apply (P, MP);
          M.apply (MP, AP);
        }
        output.numIters++; 

        const int rank = this->projectAndNormalize (iter, input, Q, H, h);
        // Save H if Ritz values are requested
        if (input.computeRitzValues && output.numRests == 0) {
          blas.COPY (iter+2, &H(0, iter), 1, &G(0, iter), 1);
        }
        // Check for negative norm
        TEUCHOS_TEST_FOR_EXCEPTION
          (STS::real (H(iter+1, iter)) < STM::zero (), std::runtime_error,
	   "At iteration " << iter << ", H(" << iter+1 << ", " << iter << ") = "
	   << H(iter+1, iter) << " < 0.");
        // Convergence check
        if (rank == 1 && H(iter+1, iter) != zero) {
          // Apply Givens rotations to new column of H and y
          this->reduceHessenburgToTriangular(iter, H, cs, sn, y);
          metric = this->getConvergenceMetric (std::abs (y(iter+1)), b_norm, input);
        }
	else {
          metric = STM::zero ();
        }
      } // end of restart-cycle 
      if (iter > 0) {
        // Compute Ritz values, if requested
        if (input.computeRitzValues && output.numRests == 0) {
	  computeRitzValues (iter, G, output.ritzValues);
          this->sortRitzValues (iter, output.ritzValues);
        }
        // Update solution
        blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                   Teuchos::NON_UNIT_DIAG, iter, 1, one,
                   H.values(), H.stride(), y.values(), y.stride());
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
      // Compute explicit unpreconditioned residual
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
          r_norm = P.norm2 (); // norm 
        }
        P.scale (one / r_norm);
        y[0] = SC {r_norm};
        for (int i=1; i < restart+1; i++) {
	  y[i] = STS::zero ();
	}
        output.numRests++;
      }
    }

    // return residual norm as B
    Tpetra::deep_copy (B, P);
    return output;
  }

private:
  Teuchos::RCP<ortho_type> ortho_;
};

template<class SC, class MV, class OP,
	 template<class, class, class> class KrylovSubclassType>
class SolverManager;

// This is the Belos::SolverManager subclass that gets registered with
// Belos::SolverFactory.
template<class SC, class MV, class OP>
using GmresSolverManager = SolverManager<SC, MV, OP, Gmres>;

/// \brief Register GmresSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_Gmres (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_GMRES_HPP
