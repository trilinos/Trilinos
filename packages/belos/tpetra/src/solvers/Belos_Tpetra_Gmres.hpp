// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename STS::magnitudeType;
  using complex_type = std::complex<real_type>;

  static void
  run (const int iter,
       Teuchos::SerialDenseMatrix<LO, SC>& G,
       std::vector<complex_type>& ritzValues);
};

// Real Scalar (SC) case.
template<class LO, class SC>
struct ComputeRitzValues<LO, SC, false> {
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename STS::magnitudeType;
  using complex_type = std::complex<real_type>;

  static void
  run (const int iter,
       Teuchos::SerialDenseMatrix<LO, SC>& G,
       std::vector<complex_type>& ritzValues)
  {
    if (ritzValues.size () < std::size_t (iter)) {
      ritzValues.resize (iter);
    }

    std::vector<real_type> WR (iter);
    std::vector<real_type> WI (iter);
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
    for (int i = 0; i < iter; ++i) {
      ritzValues[i] = complex_type { WR[i], WI[i] };
    }
  }
};

// Complex Scalar (SC) case.
template<class LO, class SC>
struct ComputeRitzValues<LO, SC, true> {
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename STS::magnitudeType;
  using complex_type = std::complex<real_type>;

  static void
  run (const int iter,
       Teuchos::SerialDenseMatrix<LO, SC>& G,
       std::vector<complex_type>& ritzValues)
  {
    if (ritzValues.size () < std::size_t (iter)) {
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

/// \brief Compute eigenvalues of the upper Hessenberg matrix from GMRES.
///
/// \param iter [in] The current GMRES iteration number
/// \param G [in] Upper Hessenberg matrix from GMRES.
/// \param ritzValues [out] The eigenvalues of G.
template<class LO, class SC>
void
computeRitzValues (const int iter,
                   Teuchos::SerialDenseMatrix<LO, SC>& G,
                   std::vector<std::complex<typename Teuchos::ScalarTraits<SC>::magnitudeType>>& ritzValues)
{
  ComputeRitzValues<LO, SC>::run (iter, G, ritzValues);
}

/// \brief Sort Ritz values computed by computeRitzValues, using the
///   Leja ordering.
template<class LO, class SC>
void
sortRitzValues (const LO m,
                std::vector<std::complex<typename Teuchos::ScalarTraits<SC>::magnitudeType>>& RR)
{
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename STS::magnitudeType;
  using complex_type = std::complex<real_type>;
  std::vector<complex_type> ritzValues (m);

  // the first Ritz is the largest
  LO i = 0;
  LO next_index = 0;
  real_type next_value = std::abs (RR[0]);
  for (i = 1; i < m; ++i) {
    if (next_value < std::abs (RR[i])) {
      next_index = i;
      next_value = std::abs (RR[i]);
    }
  }
  ritzValues[0] = RR[next_index];
  RR[next_index] = RR[0];

  i = 0;
  if (! STS::isComplex) {
    if (RR[i].imag() != 0.0) {

      if (next_index == 0) {
        ritzValues[1] = ritzValues[0];
        RR[next_index+1] = RR[1];
      } else if (next_index == m-1) {
        ritzValues[1] = ritzValues[0];
        RR[next_index-1] = RR[1];
      } else {
        real_type val1 = std::abs(std::conj(RR[next_index-1]) - RR[next_index]);
        real_type val2 = std::abs(std::conj(RR[next_index+1]) - RR[next_index]);
        if (val1 < val2) {
          ritzValues[1] = ritzValues[0];
          RR[next_index-1] = RR[1];
        } else {
          ritzValues[1] = ritzValues[0];
          RR[next_index+1] = RR[1];
        }
      }

      ritzValues[0].imag( 0.0 );
      i++;
    }
  }
  i++;

  // sort the rest of Ritz values
  for (; i < m; i++) {
    next_index = i;
    next_value = std::abs( RR[i] - ritzValues[0] );
    for (int j = 1;  j < i; j++ ) {
      next_value *= std::abs( RR[i] - ritzValues[j] );
    }

    for (int k = i+1;  k < m; k++ ) {
      real_type value = std::abs( RR[k] - ritzValues[0] );
      for (int j = 1;  j < i; j++ ) {
        value *= std::abs( RR[k] - ritzValues[j] );
      }
      if (next_value < value) {
        next_value = value;
        next_index = k;
      }
    }
    ritzValues[i] = RR[next_index];
    RR[next_index] = RR[i];

    if (! STS::isComplex) {
      if (RR[i].imag() != 0.0) {
        if (next_index == 0) {
          ritzValues[i+1] = ritzValues[i];
          RR[next_index+1] = RR[i+1];
        } else if (next_index == m-1) {
          if (i + 1 < m) {
            ritzValues[i+1] = ritzValues[i];
            RR[next_index-1] = RR[i+1];
          }
        } else {
          real_type val1 = std::abs(std::conj(RR[next_index-1]) - ritzValues[i]);
          real_type val2 = std::abs(std::conj(RR[next_index+1]) - ritzValues[i]);
          if (val1 < val2) {
            ritzValues[i+1] = ritzValues[i];
            RR[next_index-1] = RR[i+1];
          } else {
            ritzValues[i+1] = ritzValues[i];
            RR[next_index+1] = RR[i+1];
          }
        }
        ritzValues[i].imag( 0.0 );
        i++;
      }
    }
  }

  // copy back in the Leja order
  for (i = 0; i < m; i++) {
    RR[i] = ritzValues[i];
  }
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
  using real_type = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<real_type>;
  using device_type = typename MV::device_type;

  using ortho_type = Belos::OrthoManager<SC, MV>;
  using dense_matrix_type = Teuchos::SerialDenseMatrix<LO, SC>;
  using dense_vector_type = Teuchos::SerialDenseVector<LO, SC>;
  using real_vector_type = Teuchos::SerialDenseVector<LO, real_type>;

public:
  Gmres () = default;

  Gmres (const Teuchos::RCP<const OP>& A) :
    Krylov<SC, MV, OP>::Krylov (A)
  {}

  virtual ~Gmres () = default;

  virtual void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const override
  {
    Krylov<SC, MV, OP>::getParameters (params, defaultValues);

    const int resCycle = defaultValues ? 30 : this->input_.resCycle;
    params.set ("Num Blocks", resCycle);
    // FIXME (mfh 16 Sep 2018) This isn't quite right -- "Maximum
    // Restarts" in Belos refers to the number of restart cycles.
    // These Krylov solvers appear to use maxNumIters as the maximum
    // total number of iterations over all restart cycles.
    params.set ("Maximum Restarts", this->input_.maxNumIters);
  }

  virtual void
  setParameters (Teuchos::ParameterList& params) override
  {
    Krylov<SC, MV, OP>::setParameters (params);

    bool computeRitzValues 
      = params.get<bool> ("Compute Ritz Values", this->input_.computeRitzValues);

    bool needToReortho
      = params.get<bool> ("Reorthogonalize Blocks", this->input_.needToReortho);

    int resCycle = params.get<int> ("Num Blocks", this->input_.resCycle);
    TEUCHOS_TEST_FOR_EXCEPTION
        (resCycle < 0, std::invalid_argument,
         "\"Num Blocks\" (restart length) = " << resCycle << " < 0.");

    std::string orthoType
      = params.get<std::string> ("Orthogonalization", this->input_.orthoType);

    int maxOrthoSteps
      = params.get<int> ("Max Orthogonalization Passes", this->input_.maxOrthoSteps);

    this->input_.computeRitzValues = computeRitzValues;
    this->input_.needToReortho = needToReortho;
    this->input_.resCycle = resCycle;
    this->input_.orthoType = orthoType;
    this->input_.maxOrthoSteps = maxOrthoSteps;
  }

protected:
  //! Create Belos::OrthoManager instance.
  void
  setOrthogonalizer (const std::string& ortho)
  {
    if (ortho_.get () == nullptr || this->input_.orthoType != ortho) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->input_.orthoType == "", std::runtime_error,
         "Gmres: Failed to specify \"Orthogonalization\" parameter.");
      // Since setOrthogonalizer only gets called on demand, we know
      // the preconditioner (if any) at this point.  Thus, we can use
      // Belos::OrthoManagerFactory here.
      Belos::OrthoManagerFactory<SC, MV, OP> factory;
      Teuchos::RCP<Belos::OutputManager<SC>> outMan; // can be null
      Teuchos::RCP<Teuchos::ParameterList> params;   // can be null
      if (this->input_.maxOrthoSteps > 0) {
        params = Teuchos::rcp (new Teuchos::ParameterList());
        params->set ("maxNumOrthogPasses", this->input_.maxOrthoSteps);
      }
      ortho_ = factory.makeMatOrthoManager (ortho, Teuchos::null, outMan, "Belos", params);
      TEUCHOS_TEST_FOR_EXCEPTION
        (ortho_.get () == nullptr, std::runtime_error, "Gmres: Failed to "
         "create (Mat)OrthoManager of type \"" << ortho << "\".");
      this->input_.orthoType = ortho;
    }
  }

private:
  int
  projectAndNormalize (const int n,
                       const SolverInput<SC>& input,
                       MV& Q,
                       dense_matrix_type& H,
                       dense_matrix_type& /* WORK */)
  {
    return this->projectAndNormalizeBelosOrthoManager (n, Q, H, input);
  }

  // ! Apply the orthogonalization using Belos' OrthoManager
  int
  projectAndNormalizeBelosOrthoManager (int n, MV &Q, dense_matrix_type &H,
                                        const SolverInput<SC>& input)
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
      setOrthogonalizer (input.orthoType);
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (ortho_.get () == nullptr, std::logic_error, "Gmres: Failed to create "
       "(Mat)OrthoManager.  This should never happen.  Please report this bug "
       "to the Belos developers.");
    // FIXME (mfh 14 Aug 2018) For left preconditioning, we may need
    // to use the MatOrthoManager interface.
    return ortho_->projectAndNormalize (AP, h_array, r_new, Q_array);
  }

protected:
  //! Apply the orthogonalization using Belos' OrthoManager
  void
  projectBelosOrthoManager (const int n,
                            const int s,
                            MV& Q,
                            dense_matrix_type& R)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    if (n > 0) {
      // vector to be orthogonalized
      Teuchos::Range1D index(n, n+s);
      RCP<MV> Qnew = Q.subViewNonConst (index);

      // vectors to be orthogonalized against
      Teuchos::Range1D index_prev (0, n-1);
      RCP<const MV> Q_prev = MVT::CloneView (Q, index_prev);
      Teuchos::Array<RCP<const MV> > Qprev (1, Q_prev);

      // orthogonalize scalars
      RCP<dense_matrix_type> r_new =
        rcp (new dense_matrix_type (Teuchos::View, R, n, s+1, 0, 0));
      Teuchos::Array<RCP<dense_matrix_type>> r_array (1, r_new);

      if (this->ortho_ == Teuchos::null) {
        this->setOrthogonalizer (this->input_.orthoType);
      }
      this->ortho_->project (*Qnew, r_array, Qprev);
    }
  }

  int
  normalizeBelosOrthoManager (MV& Q, dense_matrix_type& R)
  {
    if (ortho_.get () == nullptr) {
      setOrthogonalizer (this->input_.orthoType);
    }
    Teuchos::RCP<dense_matrix_type> R_ptr = Teuchos::rcpFromRef (R);
    return ortho_->normalize (Q, R_ptr);
  }

  //! Reduce a column of Henssenburg matrix to triangular form
  void
  reduceHessenburgToTriangular(const int j,
                               dense_matrix_type& H,
                               std::vector<real_type>& cs,
                               std::vector<SC>& sn,
                               SC y[]) const
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
    blas.ROT (1, &y[j], 1, &y[j+1], 1, &cs[j], &sn[j]);
  }

  void
  reduceHessenburgToTriangular(const int j,
                               dense_matrix_type& H,
                               std::vector<real_type>& cs,
                               std::vector<SC>& sn,
                               dense_vector_type& y) const
  {
    this->reduceHessenburgToTriangular (j, H, cs, sn, y.values ());
  }

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

    // timers
    Teuchos::RCP< Teuchos::Time > spmvTimer  = Teuchos::TimeMonitor::getNewCounter ("Gmres::matrix-apply ");
    Teuchos::RCP< Teuchos::Time > precTimer  = Teuchos::TimeMonitor::getNewCounter ("Gmres::precondition ");
    Teuchos::RCP< Teuchos::Time > orthTimer  = Teuchos::TimeMonitor::getNewCounter ("Gmres::orthogonalize");

    Teuchos::RCP< Teuchos::Time > totalTimer = Teuchos::TimeMonitor::getNewCounter ("Gmres::total        ");
    Teuchos::TimeMonitor GmresTimer (*totalTimer);

    // initialize output parameters
    SolverOutput<SC> output {};
    output.converged = false;
    output.numRests = 0;
    output.numIters = 0;


    if (outPtr != nullptr) {
      *outPtr << "Gmres" << endl;
    }
    Indent indent1 (outPtr);
    if (outPtr != nullptr) {
      *outPtr << "Solver input:" << endl;
      Indent indentInner (outPtr);
      *outPtr << input;
    }

    real_type b_norm;  // initial residual norm
    real_type b0_norm; // initial residual norm, not left preconditioned
    real_type r_norm;
    real_type r_norm_imp;

    bool zeroOut = false; // Kokkos::View:init can take a long time on GPU?
    vec_type R (B.getMap (), zeroOut);
    vec_type Y (B.getMap (), zeroOut);
    vec_type MP (B.getMap (), zeroOut);
    MV Q (B.getMap (), restart+1, zeroOut);
    vec_type P0 = * (Q.getVectorNonConst (0));

    // initial residual (making sure R = B - Ax)
    {
      Teuchos::TimeMonitor LocalTimer (*spmvTimer);
      A.apply (X, R);
    }
    R.update (one, B, -one);
    b0_norm = R.norm2 (); // residual norm, not-preconditioned
    if (input.precoSide == "left") {
      {
        Teuchos::TimeMonitor LocalTimer (*precTimer);
        M.apply (R, P0);
      }
      b_norm = P0.norm2 (); // residual norm, left-preconditioned
    }
    else {
      Tpetra::deep_copy (P0, R);
      b_norm = b0_norm;
    }
    real_type metric = this->getConvergenceMetric (b0_norm, b0_norm, input);

    if (metric <= input.tol) {
      if (outPtr != nullptr) {
        *outPtr << "Initial guess' residual norm " << b_norm
                << " meets tolerance " << input.tol << endl;
      }
      output.absResid = b_norm;
      output.relResid = STM::one ();
      output.numIters = 0;
      output.converged = true;
      // return residual norm as B
      Tpetra::deep_copy (B, P0);
      return output;
    } else if (outPtr != NULL) {
      *outPtr << "Initial guess' residual norm " << b0_norm << endl;
    }

    Teuchos::BLAS<LO ,SC> blas;
    Teuchos::LAPACK<LO ,SC> lapack;
    dense_matrix_type  H (restart+1, restart, true);
    dense_matrix_type  G (restart+1, restart, true); // only for Ritz values
    dense_vector_type  y (restart+1, true);
    dense_matrix_type  h (restart+1, 1, true); // for reorthogonalization
    std::vector<real_type> cs (restart);
    std::vector<SC> sn (restart);

    //#define HAVE_TPETRA_DEBUG
    #ifdef HAVE_TPETRA_DEBUG
    dense_matrix_type H2 (restart+1, restart,   true);
    dense_matrix_type H3 (restart+1, restart,   true);
    #endif

    // initialize starting vector
    P0.scale (one / b_norm);
    y[0] = SC {b_norm};

    // main loop
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

      // restart cycle
      for (; iter < restart && metric > input.tol; ++iter) {
        if (outPtr != nullptr) {
          *outPtr << "Current iteration: iter=" << iter
                  << ", restart=" << restart
                  << ", metric=" << metric << endl;
          Indent indent3 (outPtr);
        }

        // AP = A*P
        vec_type P  = * (Q.getVectorNonConst (iter));
        vec_type AP = * (Q.getVectorNonConst (iter+1));
        if (input.precoSide == "none") {
          Teuchos::TimeMonitor LocalTimer (*spmvTimer);
          A.apply (P, AP);
        }
        else if (input.precoSide == "right") {
          {
            Teuchos::TimeMonitor LocalTimer (*precTimer);
            M.apply (P, MP);
          }
          {
            Teuchos::TimeMonitor LocalTimer (*spmvTimer);
            A.apply (MP, AP);
          }
        }
        else { // left
          {
            Teuchos::TimeMonitor LocalTimer (*spmvTimer);
            A.apply (P, MP);
          }
          {
            Teuchos::TimeMonitor LocalTimer (*precTimer);
            M.apply (MP, AP);
          }
        }
        output.numIters++;

        int rank = 0;
        {
          Teuchos::TimeMonitor LocalTimer (*orthTimer);
          rank = this->projectAndNormalize (iter, input, Q, H, h);
        }
        // Save H if Ritz values are requested
        if (input.computeRitzValues && output.numRests == 0) {
          blas.COPY (iter+2, &H(0, iter), 1, &G(0, iter), 1);
        }
        // Check for negative norm
        TEUCHOS_TEST_FOR_EXCEPTION
          (STS::real (H(iter+1, iter)) < STM::zero (), std::runtime_error,
           "At iteration " << iter << ", H(" << iter+1 << ", " << iter << ") = "
           << H(iter+1, iter) << " < 0.");

        #ifdef HAVE_TPETRA_DEBUG
        // Numeric check
        this->checkNumerics (outPtr, iter, iter, A, M, Q, X, B, y,
                             H, H2, H3, cs, sn, input);
        #endif

        // Convergence check
        if (rank == 1 && H(iter+1, iter) != zero) {
          // Apply Givens rotations to new column of H and y
          this->reduceHessenburgToTriangular(iter, H, cs, sn, y);
          metric = this->getConvergenceMetric (STS::magnitude (y(iter+1)), b_norm, input);
        }
        else {
          metric = STM::zero ();
        }
      } // end of restart cycle

      if (iter < restart) {
        // save the old solution, just in case explicit residual norm failed the convergence test
        Tpetra::deep_copy (Y, X);
        blas.COPY (1+iter, y.values(), 1, h.values(), 1);
      }
      r_norm_imp = STS::magnitude (y (iter)); // save implicit residual norm
      if (iter > 0) {
        // Compute Ritz values, if requested
        if (input.computeRitzValues && output.numRests == 0) {
          computeRitzValues (iter, G, output.ritzValues);
          sortRitzValues <LO, SC> (iter, output.ritzValues);
        }
        // Update solution
        blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
                   Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
                   iter, 1, one,
                   H.values(), H.stride(), y.values(), y.stride());
        Teuchos::Range1D cols(0, iter-1);
        Teuchos::RCP<const MV> Qj = Q.subView(cols);
        dense_vector_type y_iter (Teuchos::View, y.values (), iter);
        if (input.precoSide == "right") {
          //MVT::MvTimesMatAddMv (one, *Qj, y, zero, R);
          MVT::MvTimesMatAddMv (one, *Qj, y_iter, zero, R);
          {
            Teuchos::TimeMonitor LocalTimer (*precTimer);
            M.apply (R, MP);
          }
          X.update (one, MP, one);
        }
        else {
          //MVT::MvTimesMatAddMv (one, *Qj, y, one, X);
          MVT::MvTimesMatAddMv (one, *Qj, y_iter, one, X);
        }
        //y.resize (restart+1);
      }
      // Compute explicit unpreconditioned residual
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
      }
      else if (output.numIters < input.maxNumIters) {
        // Restart, only if max inner-iteration was reached.
        // Otherwise continue the inner-iteration.
        if (iter >= restart) {
          // Restart: Initialize starting vector for restart
          iter = 0;
          P0 = * (Q.getVectorNonConst (0));
          if (input.precoSide == "left") {
            {
              Teuchos::TimeMonitor LocalTimer (*precTimer);
              M.apply (R, P0);
            }
            r_norm = P0.norm2 (); // norm
          }
          else {
            // set the starting vector
            Tpetra::deep_copy (P0, R);
          }
          P0.scale (one / r_norm);
          y[0] = SC {r_norm};
          for (int i=1; i < restart+1; i++) {
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

    // return residual norm as B
    Tpetra::deep_copy (B, R);

    if (outPtr != nullptr) {
      *outPtr << "At end of solve:" << endl;
      Indent indentInner (outPtr);
      *outPtr << output;
    }
    return output;
  }

  // ! compute matrix norm
  real_type
  computeNorm(dense_matrix_type &T)
  {
    const LO ione = 1;

    LO m = T.numRows ();
    LO n = T.numCols ();
    LO minmn = (m < n ? m : n);

    LO INFO, LWORK;
    SC  U, VT, TEMP;
    real_type RWORK;
    real_vector_type S (minmn, true);
    LWORK = -1;
    Teuchos::LAPACK<LO ,SC> lapack;
    lapack.GESVD('N', 'N', m, n, T.values (), T.stride (),
                 S.values (), &U, ione, &VT, ione,
                 &TEMP, LWORK, &RWORK, &INFO);
    LWORK = Teuchos::as<LO> (STS::real (TEMP));
    dense_vector_type WORK (LWORK, true);
    lapack.GESVD('N', 'N', m, n, T.values (), T.stride (),
                 S.values (), &U, ione, &VT, ione,
                 WORK.values (), LWORK, &RWORK, &INFO);

    return S(0);
  }

  // ! Check numerics
  void 
  checkNumerics (Teuchos::FancyOStream* outPtr,
                 const int iter,
                 const int check,
                 const OP& A,
                 const OP& M,
                 const MV& Q,
                 const vec_type& X,
                 const vec_type& B,
                 const dense_vector_type& y,
                 const dense_matrix_type& H,
                       dense_matrix_type& H2,
                       dense_matrix_type& H3,
                       std::vector<real_type>& cs,
                       std::vector<SC>& sn,
                 const SolverInput<SC>& input)
  {
    Teuchos::BLAS<LO ,SC> blas;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();

    // quick return
    if (outPtr == nullptr) {
      return;
    }

    // save H (before convert it to triangular form)
    for (int i = 0; i <= check+1; i++) {
      H2 (i, check) = H (i, check);
      H3 (i, check) = H (i, check);
    }
    // reduce H3 to triangular
    dense_vector_type y2 (check+2, true); // to save y
    blas.COPY (check+1, y.values(), 1, y2.values(), 1);
    this->reduceHessenburgToTriangular (check, H3, cs, sn, y2.values ());

    #if 0
    printf( " > checkNumeric(iter = %d, check = %d)\n",iter,check );
    /*auto Q_lcl = Q.getLocalViewHost();
    printf( " Q = [\n" );
    for (int i = 0; i < (int)Q_lcl.extent(0); i++) {
      for (int j = 0; j <= check+1; j++) printf( "%.16e ",Q_lcl(i,j) );
      printf("\n" );
    }
    printf("];\n" );*/

    printf( " H2 = [\n" );
    for (int i = 0; i <= check+1; i++) {
      for (int j = 0; j <= check; j++) {
        printf( "%.16e ", H2 (i, j) );
      }
      printf( "\n" );
    }
    printf( "];\n" );
    /*printf( " H3 = [\n" );
    for (int i = 0; i <= check+1; i++) {
      for (int j = 0; j <= check; j++) {
        printf( "%e ", H3 (i, j) );
      }
      printf( "\n" );
    }
    printf( "];\n" );*/
    //for (int i = 0; i <= check+1; i++) printf( "%d %e -> %e\n",i,y(i),y2(i) );
    fflush(stdout);
    #endif

    // save X
    vec_type X2 (X.getMap (), false);    // to save X
    Tpetra::deep_copy (X2, X);

    // implicit residual norm
    real_type r_norm_imp = STS::magnitude (y2 (check+1));
    y2.resize (check+1);

    // --------------------------------------
    // compute explicit residual norm
    // > Update solution, X += Q*y
    vec_type R  (X.getMap (), false);  // to save X
    vec_type MP (X.getMap (), false);  // to save X
    blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
               Teuchos::NON_UNIT_DIAG, check+1, 1, one,
               H3.values(), H3.stride(), y2.values (), y2.stride ());
    Teuchos::Range1D cols(0, check);
    Teuchos::RCP<const MV> Qj = Q.subView(cols);
    if (input.precoSide == "right") {
      MVT::MvTimesMatAddMv (one, *Qj, y2, zero, R);
      M.apply (R, MP);
      X2.update (one, MP, one);
    }
    else {
      MVT::MvTimesMatAddMv (one, *Qj, y2, one, X2);
    }
    // > Compute explicit residual vector
    A.apply (X2, R);
    R.update (one, B, -one);
    real_type r_norm = R.norm2 (); // residual norm
    real_type b_norm = B.norm2 (); // originial noorm

    // --------------------------------------
    // compute orthogonality error, norm (Q'*Q-I)
    real_type ortho_error (0.0);
    {
      Teuchos::Range1D index_prev(0, check+1);
      Teuchos::RCP<const MV> Q_prev = MVT::CloneView (Q, index_prev);
      dense_matrix_type T (check+2, check+2, true);
      MVT::MvTransMv(one, *Q_prev, *Q_prev, T);
      for (int i = 0; i < check+2; i++) {
        T (i, i) -= one;
      }
      /*printf("Y=[\n");
      for (int i = 0; i < check+2; i++) {
        for (int j = 0; j < check+2; j++) {
          printf("%e ",T (i, j));
        }
        printf("\n");
      }
      printf("];\n");*/
      ortho_error = computeNorm(T);
    }

    // --------------------------------------
    // compute Arnoldi representation error
    real_type repre_error (0.0);
    real_type proje_error (0.0);
    {
      // > compute AQ=A*Q
      MV AQ (Q.getMap (), check+1);
      {
        Teuchos::Range1D index_prev(0, check);
        Teuchos::RCP<const MV> Q_prev = MVT::CloneView (Q, index_prev);
        if (input.precoSide == "left") {
          MV AM (Q.getMap (), check+1);
          A.apply (*Q_prev, AM);
          M.apply (AM, AQ);
        } else if (input.precoSide == "right") {
          MV AM (Q.getMap (), check+1);
          M.apply (*Q_prev, AM);
          A.apply (AM, AQ);
        } else {
          A.apply (*Q_prev, AQ);
        }
      }

      // > compute HH = H - Q'*A*Q and AQ = Q*H - AQ
      dense_matrix_type HH(check+2, check+1, true);
      {
        Teuchos::Range1D index_prev(0, check+1);
        Teuchos::RCP<const MV> Q_prev = MVT::CloneView (Q, index_prev);
        auto H_new = rcp (new dense_matrix_type (Teuchos::View, H2, check+2, check+1, 0, 0));

        // > compute H - Q'*A*Q
        MVT::MvTransMv(one, *Q_prev, AQ, HH);
        for (int j = 0; j < check+1; j++) {
          blas.AXPY (check+2, -one, &(H2(0, j)), 1, &(HH(0, j)), 1);
        }

        // > compute AQ = Q*H - AQ
        MVT::MvTimesMatAddMv (one, *Q_prev, *H_new, -one, AQ);
      }

      // > compute norm, norm(HH) and sqrt(norm(AQ'*AQ))
      {
        proje_error = computeNorm(HH);

        dense_matrix_type T (check+1, check+1, true);
        MVT::MvTransMv(one, AQ, AQ, T);
        repre_error = STM::squareroot (computeNorm(T));
      }
    }

    *outPtr << " > iter = " << iter
            << ", check = " << check
            << ", Right-hand side norm: "
            << b_norm
            << ", Implicit and explicit residual norms: "
            << r_norm_imp << ", " << r_norm
            << " -> "
            << r_norm_imp/b_norm << ", " << r_norm/b_norm
            << ", Ortho error: "
            << ortho_error
            << ", Arnoldi representation error: "
            << repre_error
            << ", Projection error: "
            << proje_error
            << std::endl;
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
