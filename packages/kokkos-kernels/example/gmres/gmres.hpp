/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Jennifer Loe (jloe@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <math.h>
#include "Kokkos_ArithTraits.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas3_trsm.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_Preconditioner.hpp>
#include "KokkosKernels_Error.hpp"

////////////////////////////////////////////////////////////////////////////////
// libstdc++ half_t overloads
////////////////////////////////////////////////////////////////////////////////
namespace Kokkos {
namespace Issue1172WorkAround {
#if !KOKKOS_HALF_T_IS_FLOAT
Kokkos::Experimental::half_t fabs(Kokkos::Experimental::half_t arg) {
  using AT = Kokkos::Details::ArithTraits<Kokkos::Experimental::half_t>;
  return AT::abs(arg);
}

Kokkos::Experimental::half_t fabs(
    Kokkos::complex<Kokkos::Experimental::half_t> arg) noexcept {
  return Kokkos::Experimental::half_t(Kokkos::abs(
      Kokkos::complex<double>((double)arg.real(), (double)arg.imag())));
}
#endif  // KOKKOS_HALF_T_IS_FLOAT

#if !KOKKOS_BHALF_T_IS_FLOAT
Kokkos::Experimental::bhalf_t fabs(Kokkos::Experimental::bhalf_t arg) {
  using AT = Kokkos::Details::ArithTraits<Kokkos::Experimental::bhalf_t>;
  return AT::abs(arg);
}

Kokkos::Experimental::bhalf_t fabs(
    Kokkos::complex<Kokkos::Experimental::bhalf_t> arg) noexcept {
  return Kokkos::Experimental::bhalf_t(Kokkos::abs(
      Kokkos::complex<double>((double)arg.real(), (double)arg.imag())));
}
#endif  // KOKKOS_BHALF_T_IS_FLOAT

// This fabs wrapper was added to resolve:
// https://github.com/kokkos/kokkos-kernels/issues/1172
template <class T>
KOKKOS_INLINE_FUNCTION T fabs(const Kokkos::complex<T> &x) {
  return Kokkos::abs(x);
}
}  // namespace Issue1172WorkAround
}  // namespace Kokkos

// This struct is returned to the user to give solver
// statistics and convergence status.
struct GmresStats {
  int numIters;
  double endRelRes;
  enum FLAG { Conv, NoConv, LOA };
  FLAG convFlagVal;
  std::string convFlag() {
    switch (convFlagVal) {
      case Conv: return "Converged";
      case NoConv: return "Not Converged";
      case LOA: return "Solver has had loss of accuracy.";
      default: return "Flag not defined.";
    }
  }
};

// This struct allows the user to pass in several
// options to the solver.
template <class ScalarType>
struct GmresOpts {
  typename Kokkos::Details::ArithTraits<ScalarType>::mag_type tol;
  int m;
  int maxRestart;
  std::string ortho;
  std::string precSide;

  GmresOpts<ScalarType>() : tol(1e-8), m(50), maxRestart(50), ortho("CGS2") {}
};

template <class ScalarType, class Layout, class EXSP, class OrdinalType = int>
GmresStats gmres(
    const KokkosSparse::CrsMatrix<ScalarType, OrdinalType, EXSP> &A,
    const Kokkos::View<ScalarType *, Layout, EXSP> &B,
    Kokkos::View<ScalarType *, Layout, EXSP> &X,
    const GmresOpts<ScalarType> &opts,
    const KokkosSparse::Experimental::Preconditioner<
        ScalarType, Layout, EXSP, OrdinalType> *const M = NULL) {
  using namespace Kokkos::Issue1172WorkAround;  // For 'fabs' wrappers above
  Kokkos::Profiling::pushRegion("GMRES::TotalTime:");
  typedef Kokkos::Details::ArithTraits<ScalarType> AT;
  typedef typename AT::val_type
      ST;  // So this code will run with ScalarType = std::complex<T>.
  typedef typename AT::mag_type MT;
  ST one  = AT::one();
  ST zero = AT::zero();

  typedef Kokkos::View<ST *, Layout, EXSP> ViewVectorType;
  typedef Kokkos::View<ST *, Kokkos::LayoutRight, Kokkos::HostSpace>
      ViewHostVectorType;
  typedef Kokkos::View<ST **, Layout, EXSP> ViewMatrixType;

  unsigned int n = A.numRows();

  // Store solver options:
  const int m = opts.m;

  // Check compatibility of dimensions at run time.
  if (n != unsigned(A.numCols())) {
    std::ostringstream os;
    os << "gmres: A must be a square matrix: "
       << "numRows: " << n << "  numCols: " << A.numCols();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (X.extent(0) != B.extent(0) || X.extent(0) != n) {
    std::ostringstream os;
    os << "gmres: Dimensions of A, X, and B do not match: "
       << "A: " << n << " x " << n << ", X: " << X.extent(0)
       << "x 1, B: " << B.extent(0) << " x 1";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  // Check parameter validity:
  if (m <= 0) {
    throw std::invalid_argument(
        "gmres: Please choose restart size m greater than zero.");
  }
  if (opts.maxRestart < 0) {
    throw std::invalid_argument(
        "gmres: Please choose maxRestart greater than zero.");
  }

  bool converged = false;
  int cycle      = 0;  // How many times have we restarted?
  int numIters   = 0;  // Number of iterations within the cycle before
                       // convergence.
  MT nrmB, trueRes, relRes, shortRelRes;
  GmresStats myStats;

  std::cout << "Convergence tolerance is: " << opts.tol << std::endl;

  ViewVectorType Xiter(
      "Xiter", n);  // Intermediate solution at iterations before restart.
  ViewVectorType Res(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Res"),
                     n);  // Residual vector
  ViewVectorType Wj(Kokkos::view_alloc(Kokkos::WithoutInitializing, "W_j"),
                    n);  // Tmp work vector 1
  ViewVectorType Wj2(Kokkos::view_alloc(Kokkos::WithoutInitializing, "W_j2"),
                     n);  // Tmp work vector 2
  ViewHostVectorType GVec_h(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "GVec"), m + 1);
  ViewMatrixType GLsSoln(
      "GLsSoln", m,
      1);  // LS solution vec for Givens Rotation. Must be 2-D for trsm.
  typename ViewMatrixType::HostMirror GLsSoln_h = Kokkos::create_mirror_view(
      GLsSoln);  // This one is needed for triangular solve.
  ViewHostVectorType CosVal_h("CosVal", m);
  ViewHostVectorType SinVal_h("SinVal", m);
  ViewMatrixType V(Kokkos::view_alloc(Kokkos::WithoutInitializing, "V"), n,
                   m + 1);
  ViewMatrixType VSub;  // Subview of 1st m cols for updating soln.
  ViewVectorType orthoTmp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "orthoTmp"), m);

  ViewMatrixType H("H", m + 1,
                   m);  // H matrix on device. Also used in Arn Rec debug.
  typename ViewMatrixType::HostMirror H_h =
      Kokkos::create_mirror_view(H);  // Make H into a host view of H.

  // Compute initial residuals:
  nrmB = KokkosBlas::nrm2(B);
  Kokkos::deep_copy(Res, B);

  // This is initial true residual, so don't need prec here.
  KokkosSparse::spmv("N", one, A, X, zero, Wj);  // wj = Ax
  KokkosBlas::axpy(-one, Wj, Res);               // res = res-Wj = b-Ax.
  trueRes = KokkosBlas::nrm2(Res);
  if (nrmB != 0) {
    relRes = trueRes / nrmB;
  } else if (trueRes == 0) {
    relRes = trueRes;
  } else {  // B is zero, but X has wrong initial guess.
    Kokkos::deep_copy(X, 0.0);
    relRes = 0;
  }
  shortRelRes = relRes;
  std::cout << "Initial relative residual is: " << relRes << std::endl;
  if (relRes < opts.tol) {
    converged = true;
  }

  while (!converged && cycle <= opts.maxRestart && shortRelRes >= 1e-14) {
    GVec_h(0) = trueRes;

    // Run Arnoldi iteration:
    auto Vj = Kokkos::subview(V, Kokkos::ALL, 0);
    Kokkos::deep_copy(Vj, Res);
    KokkosBlas::scal(Vj, one / trueRes, Vj);  // V0 = V0/norm(V0)

    for (int j = 0; j < m; j++) {
      if (M != NULL) {                                   // Apply Right prec
        M->apply(Vj, Wj2);                               // wj2 = M*Vj
        KokkosSparse::spmv("N", one, A, Wj2, zero, Wj);  // wj = A*MVj = A*Wj2
      } else {
        KokkosSparse::spmv("N", one, A, Vj, zero, Wj);  // wj = A*Vj
      }
      Kokkos::Profiling::pushRegion("GMRES::Orthog:");
      if (opts.ortho == "MGS") {
        for (int i = 0; i <= j; i++) {
          auto Vi   = Kokkos::subview(V, Kokkos::ALL, i);
          H_h(i, j) = KokkosBlas::dot(Vi, Wj);   // Vi^* Wj
          KokkosBlas::axpy(-H_h(i, j), Vi, Wj);  // wj = wj-Hij*Vi
        }
        auto Hj_h = Kokkos::subview(H_h, Kokkos::make_pair(0, j + 1), j);
      } else if (opts.ortho == "CGS2") {
        auto V0j = Kokkos::subview(V, Kokkos::ALL, Kokkos::make_pair(0, j + 1));
        auto Hj  = Kokkos::subview(H, Kokkos::make_pair(0, j + 1), j);
        auto Hj_h = Kokkos::subview(H_h, Kokkos::make_pair(0, j + 1), j);
        KokkosBlas::gemv("C", one, V0j, Wj, zero, Hj);  // Hj = Vj^T * wj
        KokkosBlas::gemv("N", -one, V0j, Hj, one, Wj);  // wj = wj - Vj * Hj

        // Re-orthog CGS:
        auto orthoTmpSub =
            Kokkos::subview(orthoTmp, Kokkos::make_pair(0, j + 1));
        KokkosBlas::gemv("C", one, V0j, Wj, zero,
                         orthoTmpSub);  // tmp (Hj) = Vj^T * wj
        KokkosBlas::gemv("N", -one, V0j, orthoTmpSub, one,
                         Wj);                    // wj = wj - Vj * tmp
        KokkosBlas::axpy(one, orthoTmpSub, Hj);  // Hj = Hj + tmp
        Kokkos::deep_copy(Hj_h, Hj);
      } else {
        throw std::invalid_argument(
            "Invalid argument for 'ortho'.  Please use 'CGS2' or 'MGS'.");
      }

      MT tmpNrm     = KokkosBlas::nrm2(Wj);
      H_h(j + 1, j) = tmpNrm;
      if (tmpNrm > 1e-14) {
        Vj = Kokkos::subview(V, Kokkos::ALL, j + 1);
        KokkosBlas::scal(Vj, one / H_h(j + 1, j), Wj);  // Vj = Wj/H(j+1,j)
      }
      Kokkos::Profiling::popRegion();

      // Givens for real and complex (See Alg 3 in "On computing Givens
      // rotations reliably and efficiently" by Demmel, et. al. 2001) Apply
      // Givens rotation and compute shortcut residual:
      for (int i = 0; i < j; i++) {
        ST tempVal = CosVal_h(i) * H_h(i, j) + SinVal_h(i) * H_h(i + 1, j);
        H_h(i + 1, j) =
            -AT::conj(SinVal_h(i)) * H_h(i, j) + CosVal_h(i) * H_h(i + 1, j);
        H_h(i, j) = tempVal;
      }
      ST f          = H_h(j, j);
      ST g          = H_h(j + 1, j);
      MT f2         = AT::real(f) * AT::real(f) + AT::imag(f) * AT::imag(f);
      MT g2         = AT::real(g) * AT::real(g) + AT::imag(g) * AT::imag(g);
      ST fg2        = f2 + g2;
      ST D1         = one / AT::sqrt(f2 * fg2);
      CosVal_h(j)   = f2 * D1;
      fg2           = fg2 * D1;
      H_h(j, j)     = f * fg2;
      SinVal_h(j)   = f * D1 * AT::conj(g);
      H_h(j + 1, j) = zero;

      GVec_h(j + 1) = GVec_h(j) * (-AT::conj(SinVal_h(j)));
      GVec_h(j)     = GVec_h(j) * CosVal_h(j);
      shortRelRes   = fabs(GVec_h(j + 1)) / nrmB;

      std::cout << "Shortcut relative residual for iteration "
                << j + (cycle * m) << " is: " << shortRelRes << std::endl;
      if (tmpNrm <= 1e-14 && shortRelRes >= opts.tol) {
        throw std::runtime_error(
            "GMRES has experienced lucky breakdown, but the residual has not converged.\n\
                                  Solver terminated without convergence.");
      }
      if (AT::isNan(ST(shortRelRes))) {
        throw std::runtime_error(
            "gmres: Relative residual is nan. Terminating solver.");
      }

      // If short residual converged, or time to restart, check true residual
      if (shortRelRes < opts.tol || j == m - 1) {
        // Compute least squares soln with Givens rotation:
        auto GLsSolnSub_h = Kokkos::subview(
            GLsSoln_h, Kokkos::ALL,
            0);  // Original view has rank 2, need a rank 1 here.
        auto GVecSub_h = Kokkos::subview(GVec_h, Kokkos::make_pair(0, m));
        Kokkos::deep_copy(GLsSolnSub_h,
                          GVecSub_h);  // Copy LS rhs vec for triangle solve.
        auto GLsSolnSub2_h = Kokkos::subview(
            GLsSoln_h, Kokkos::make_pair(0, j + 1), Kokkos::ALL);
        auto H_Sub_h = Kokkos::subview(H_h, Kokkos::make_pair(0, j + 1),
                                       Kokkos::make_pair(0, j + 1));
        KokkosBlas::trsm("L", "U", "N", "N", one, H_Sub_h,
                         GLsSolnSub2_h);  // GLsSoln = H\GLsSoln
        Kokkos::deep_copy(GLsSoln, GLsSoln_h);

        // Update solution and compute residual with Givens:
        VSub = Kokkos::subview(V, Kokkos::ALL, Kokkos::make_pair(0, j + 1));
        Kokkos::deep_copy(Xiter,
                          X);  // Can't overwrite X with intermediate solution.
        auto GLsSolnSub3 =
            Kokkos::subview(GLsSoln, Kokkos::make_pair(0, j + 1), 0);
        if (M != NULL) {  // Apply right prec to correct soln.
          KokkosBlas::gemv("N", one, VSub, GLsSolnSub3, zero,
                           Wj);                // wj = V(1:j+1)*lsSoln
          M->apply(Wj, Xiter, "N", one, one);  // Xiter = M*wj + X
        } else {
          KokkosBlas::gemv("N", one, VSub, GLsSolnSub3, one,
                           Xiter);  // x_iter = x + V(1:j+1)*lsSoln
        }
        KokkosSparse::spmv("N", one, A, Xiter, zero, Wj);  // wj = Ax
        Kokkos::deep_copy(Res, B);                         // Reset r=b.
        KokkosBlas::axpy(-one, Wj, Res);                   // r = b-Ax.
        trueRes = KokkosBlas::nrm2(Res);
        relRes  = trueRes / nrmB;
        std::cout << "True relative residual for iteration " << j + (cycle * m)
                  << " is : " << relRes << std::endl;
        numIters = j + 1;

        if (relRes < opts.tol) {
          converged = true;
          Kokkos::deep_copy(
              X, Xiter);  // Final solution is the iteration solution.
          break;          // End Arnoldi iteration.
        } else if (shortRelRes < 1e-30) {
          std::cout << "Short residual has converged to machine zero, but true "
                       "residual is not converged.\n"
                    << "You may have given GMRES a singular matrix. Ending the "
                       "GMRES iteration."
                    << std::endl;
          break;  // End Arnoldi iteration; we can't make any more progress.
        }
      }

    }  // end Arnoldi iter.

    cycle++;

    // This is the end, or it's time to restart. Update solution to most
    // recent vector.
    Kokkos::deep_copy(X, Xiter);
  }

  std::cout << "Ending relative residual is: " << relRes << std::endl;
  myStats.endRelRes = static_cast<double>(relRes);
  if (converged) {
    std::cout << "Solver converged! " << std::endl;
    myStats.convFlagVal = GmresStats::FLAG::Conv;
  } else if (shortRelRes < opts.tol) {
    std::cout << "Shortcut residual converged, but solver experienced a loss "
                 "of accuracy."
              << std::endl;
    myStats.convFlagVal = GmresStats::FLAG::LOA;
  } else {
    std::cout << "Solver did not converge. :( " << std::endl;
    myStats.convFlagVal = GmresStats::FLAG::NoConv;
  }
  if (cycle > 0) {
    myStats.numIters = (cycle - 1) * m + numIters;
  } else {
    myStats.numIters = 0;
  }
  std::cout << "The solver completed " << myStats.numIters << " iterations."
            << std::endl;

  Kokkos::Profiling::popRegion();
  return myStats;
}
