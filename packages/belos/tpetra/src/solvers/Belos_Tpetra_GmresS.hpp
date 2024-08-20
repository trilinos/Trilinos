// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_GMRES_S_HPP
#define BELOS_TPETRA_GMRES_S_HPP

#include "Belos_Tpetra_Gmres.hpp"
#include "Belos_Tpetra_UpdateNewton.hpp"

namespace BelosTpetra {
namespace Impl {

template<class SC = Tpetra::Operator<>::scalar_type,
         class MV = Tpetra::MultiVector<SC>,
         class OP = Tpetra::Operator<SC>>
class GmresS : public Gmres<SC, MV, OP> {
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
  GmresS () = default;

  GmresS (const Teuchos::RCP<const OP>& A) :
    base_type::Gmres (A)
  {}

  virtual ~GmresS () = default;

  virtual void
  setParameters (Teuchos::ParameterList& params) {
    base_type::setParameters (params);

    int stepSize = this->input_.stepSize;
    if (params.isParameter ("Step Size")) {
      stepSize = params.get<int> ("Step Size");
    }
    this->input_.stepSize = stepSize;
  }

private:
  // If no preconditioning:    AP := A*P.
  // If right preconditioning: MP := M(P), AP = A*MP  (= A*M(P)).
  // if left preconditioning:  MP := A*P,  AP = M(MP) (= M(A*P)).
  void
  applyPreconditionedOperator (MV& AP,
                               MV& MP,
                               const OP& A,
                               const OP& M,
                               const MV& P,
                               const SolverInput<SC>& input) const
  {
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
    const int s = input.stepSize;
    const SC zero = STS::zero ();
    const SC one  = STS::one ();
    SolverOutput<SC> output {};

    if (outPtr != nullptr) {
      *outPtr << "GMRES(s)" << endl;
    }
    Indent indent1 (outPtr);
    if (outPtr != nullptr) {
      *outPtr << "Solver input:" << endl;
      Indent indentInner (outPtr);
      *outPtr << input;
    }

    // Initial unpreconditioned residual norm.
    // This is only useful to help us judge convergence.
    //
    // FIXME (mfh 18 Sep 2018) It's a waste to recompute this, since
    // standard GMRES already does.
    mag_type b0_norm = STM::zero ();
    { // FIXME (could get this from standard GMRES)
      vec_type R (B, Teuchos::Copy);
      A.apply (X, R, Teuchos::NO_TRANS, STS::one (), -STS::one ());
      b0_norm = R.norm2 ();
    }

    // Invoke standard Gmres for the first restart cycle, to compute
    // Ritz values for use as Newton shifts
    vec_type R (B, Teuchos::Copy);
    {
      if (outPtr != nullptr) {
        *outPtr << "Run standard GMRES for first restart cycle" << endl;
      }
      SolverInput<SC> input_gmres = input;
      input_gmres.maxNumIters = s;
      input_gmres.computeRitzValues = true;
      output = Gmres<SC, MV, OP>::solveOneVec (outPtr, X, R, A, M,
                                               input_gmres);
      if (output.converged) {
        return output; // standard GMRES converged
      }
      output.numRests++;
    }

    // Standard GMRES did some work for us.
    // output.absResid is _our_ initial residual norm.
    mag_type b_norm = output.absResid;

    vec_type MP (B.getMap ());
    MV Q (B.getMap (), s+1);
    vec_type P0 = * (Q.getVectorNonConst (0));

    mag_type r_norm = STM::zero ();    
    if (input.precoSide == "left") {
      M.apply (R, P0);
      r_norm = P0.norm2 (); // residual norm
    }
    else {
      Tpetra::deep_copy (P0, R);
      r_norm = output.absResid;
    }
    mag_type metric = this->getConvergenceMetric (r_norm, b0_norm, input);

    Teuchos::BLAS<LO, SC> blas;
    Teuchos::LAPACK<LO, SC> lapack;
    dense_matrix_type H (s+1, s, true); // upper Hessenberg matrix
    dense_matrix_type G (s+1, s+1, true); // Upper-tri matrix from ortho
    dense_matrix_type BB (s, s, true); // Change-of-basis matrix
    dense_vector_type y (s+1, true);
    std::vector<mag_type> cs (s);
    std::vector<SC> sn (s);

    // Construct s x s upper part of change-of-basis matrix.
    if (outPtr != nullptr) {
      *outPtr << "Construct B" << endl;
    }
    for (int j = 0; j < s; ++j) {
      const auto newtonShift = output.ritzValues[j];
      BB(j, j) = SC {newtonShift.real ()};
      // Correction for real-arithmetic Newton basis
      if (j != 0) {
	BB(j-1, j) = -SC {newtonShift.imag () * newtonShift.imag ()};
      }
      if (j + 1 < s) {
	BB(j+1, j) = STS::one ();
      }
    }

    // initialize starting vector
    P0.scale (one / b_norm);
    y[0] = SC {b_norm};

    if (outPtr != nullptr) {
      *outPtr << "Main loop" << endl;
    }
    for (int k = 0; k < s; ++k) {
      vec_type P = * (Q.getVectorNonConst (k));
      vec_type AP = * (Q.getVectorNonConst (k+1));

      applyPreconditionedOperator (AP, MP, A, M, P, input);
      const auto newtonShift = output.ritzValues[k];
      // AP := AP - real(theta_k) * P
      AP.update (-newtonShift.real (), P, STS::one ());
      // Correction for real-arithmetic Newton basis
      if (k > 0 && newtonShift.imag () != STM::zero ()) {
	vec_type P_prev = * (Q.getVectorNonConst (k-1));
	AP.update (newtonShift.imag () * newtonShift.imag (),
		   P_prev, STS::one ());
      }
    }
    const int rank = this->normalizeBelosOrthoManager (Q, G);
    if (outPtr != nullptr) {
      *outPtr << "Rank of s-step basis: " << rank << endl;
    }

    // For notation, see Equation (3.6) in Hoemmen's dissertation.
    // Replace R in the dissertation with G here.

    const SC rho = G(s,s);
    const SC rho_tilde = G(s-1,s-1);
    if (outPtr != nullptr) {
      *outPtr << "G(s-1,s-1)=" << rho_tilde << ", G(s,s)=" << rho << endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (rho_tilde == STS::zero () || rho == STS::zero (),
       std::runtime_error, "Basis is rank deficient");

    // G * B * G^{-1}
    dense_matrix_type H_sxs (Teuchos::View, H.values (), H.stride (), s, s);
    dense_matrix_type G_sxs (Teuchos::View, G.values (), G.stride (), s, s);
    blas.GEMM (Teuchos::NO_TRANS, Teuchos::NO_TRANS, s, s, s, STS::one (),
	       G.values (), G.stride (),
	       BB.values (), BB.stride (),
	       STS::one (), H_sxs.values (), H_sxs.stride ());
    blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
	       Teuchos::NON_UNIT_DIAG, s, s, STS::one (),
	       G_sxs.values (), G_sxs.stride (),
	       H_sxs.values (), H_sxs.stride ());

    dense_vector_type z (Teuchos::View, &G(0, s), s);
    const SC b = STS::one (); // for Newton basis
    const SC z_scalingFactor = b / rho_tilde;
      
    for (int i = 0; i < s; ++i) {
      H_sxs(i, s-1) += z_scalingFactor * z(i);
    }
    H(s, s-1) = z_scalingFactor * rho;

    if (outPtr != nullptr) {
      *outPtr << "Finished computing H; H(s,s-1)=" << H(s,s-1) << endl;
    }

    // Convergence check
    if (H(s, s-1) != STS::zero ()) {
      // Reduce H to upper tri, and apply rotations to y
      for (int j = 0; j < s; ++j) {
	this->reduceHessenburgToTriangular (j, H, cs, sn, y);
      }
      if (outPtr != nullptr) {
	*outPtr << "y(s)=" << y(s) << ", b0_norm=" << b0_norm << endl;
      }
      metric = this->getConvergenceMetric (STS::magnitude (y(s)),
					   b0_norm, input);
    }
    else {
      metric = STM::zero ();
    }

    dense_vector_type y_s (Teuchos::View, y.values (), s);
    // y_s := H(1:s,1:s)^{-1} y_s
    blas.TRSM (Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI,
               Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
	       s, 1, one, H_sxs.values (), H_sxs.stride (),
	       y_s.values (), y_s.stride ());
    Teuchos::Range1D cols (0, s-1);
    Teuchos::RCP<const MV> Qj = Q.subView (cols);

    if (input.precoSide == "right") {
      MVT::MvTimesMatAddMv (one, *Qj, y_s, zero, R);
      M.apply (R, MP);
      X.update (one, MP, one);
    }
    else {
      MVT::MvTimesMatAddMv (one, *Qj, y_s, one, X);
    }

    // Compute explicit unpreconditioned residual
    P0 = * (Q.getVectorNonConst (0));
    A.apply (X, P0);
    P0.update (one, B, -one);
    r_norm = P0.norm2 (); // residual norm

    output.numIters += s;
    output.numRests++;
    output.absResid = r_norm;
    output.relResid = r_norm / b0_norm;
    metric = this->getConvergenceMetric (r_norm, b0_norm, input);
    output.converged = (metric <= input.tol);

    Tpetra::deep_copy (B, P0); // return residual norm as B

    if (outPtr != nullptr) {
      *outPtr << "At end of GMRES(s) cycle:" << endl;
      Indent indentInner (outPtr);
      *outPtr << output;
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
using GmresSSolverManager = SolverManager<SC, MV, OP, GmresS>;

/// \brief Register GmresSSolverManager for all enabled Tpetra
///   template parameter combinations.
void register_GmresS (const bool verbose);

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_GMRES_S_HPP
