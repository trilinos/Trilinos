//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DEFAULTSPARSESOLVE_KERNELOPS_HPP
#define KOKKOS_DEFAULTSPARSESOLVE_KERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif


namespace Kokkos {

  // 
  // Matrix formatting and mat-vec options
  // Applies to all four operations below
  // 
  // unitDiag indicates whether we neglect the diagonal row entry and scale by it
  // or utilize all row entries and implicitly scale by a unit diagonal (i.e., don't need to scale)
  // upper (versus lower) will determine the ordering of the solve and the location of the diagonal
  // 
  // upper -> diagonal is first entry on row
  // lower -> diagonal is last entry on row
  // 

  template <class Scalar, class OffsetType, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseSolveOp {
    // mat data
    const OffsetType *offs;
    const Ordinal    *inds;
    const Scalar     *vals;
    Ordinal numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar  *x;
    const RangeScalar *y;
    Ordinal numRHS, xstride, ystride;

    inline KERNEL_PREFIX void execute() {
      // solve for X in A * X = Y
      // 
      // upper triangular requires backwards substition, solving in reverse order
      // must unroll the last iteration, because decrement results in wrap-around
      // 
      if (upper && unitDiag) {
        // upper + unit
        for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+numRows-1] = (DomainScalar)y[j*ystride+numRows-1];
        for (Ordinal r=2; r < numRows+1; ++r) {
          const Ordinal row = numRows - r; // for row=numRows-2 to 0 step -1
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] = (DomainScalar)y[j*ystride+row];
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] -= (DomainScalar)val * x[j*xstride+ind];
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit
        DomainScalar diag = vals[offs[numRows-1]];
        for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+numRows-1] = (DomainScalar)( y[j*ystride+numRows-1] / diag );
        for (Ordinal r=2; r < numRows+1; ++r) {
          const Ordinal row = numRows - r; // for row=numRows-2 to 0 step -1
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          // capture and skip the diag
          ++i;
          diag = *v++;
          //
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] = (DomainScalar) y[j*ystride+row];
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] -= (DomainScalar) val * x[j*xstride+ind];
          }
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] /= diag;
        }
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (Ordinal j=0; j<numRHS; ++j) x[j*xstride] = (DomainScalar) y[j*ystride];
        for (Ordinal row=1; row < numRows; ++row) {
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] = (DomainScalar) y[j*ystride+row];
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] -= (DomainScalar) val * x[j*xstride+ind];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit
        DomainScalar diag = vals[0];
        for (Ordinal j=0; j<numRHS; ++j) x[j*xstride] = (DomainScalar)( y[j*ystride] / (RangeScalar)diag );
        for (Ordinal row=1; row < numRows; ++row) {
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          // skip the diag
          --ie;
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] = (DomainScalar) y[j*ystride+row];
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] -= (DomainScalar) val * x[j*xstride+ind];
          }
          diag = *v;
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] /= diag;
        }
      }
    }
  };

  template <class Scalar, class OffsetType, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeSolveOp {
    // mat data
    const OffsetType *offs;
    const Ordinal    *inds;
    const Scalar     *vals;
    Ordinal numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar  *x;
    const RangeScalar *y;
    Ordinal numRHS, xstride, ystride;

    inline KERNEL_PREFIX void execute() {
      // solve for X in A^H * X = Y
      // 
      // put y into x and solve system in-situ
      // this is copy-safe, in the scenario that x and y point to the same location.
      //
      for (Ordinal rhs=0; rhs < numRHS; ++rhs) {
        for (Ordinal row=0; row < numRows; ++row) {
          x[rhs*xstride+row] = y[rhs*xstride+row];
        }
      }
      // 
      if (upper && unitDiag) {
        // upper + unit
        for (Ordinal row=0; row < numRows-1; ++row) {
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+ind] -= (DomainScalar)val * x[j*xstride+row];
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit; diag is first element in row
        DomainScalar diag;
        for (Ordinal row=0; row < numRows-1; ++row) {
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          // capture and skip the diag
          ++i;
          diag = *v++;
          //
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] /= diag;
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+ind] -= (DomainScalar)val * x[j*xstride+row];
          }
        }
        diag = vals[offs[numRows-1]];
        for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+numRows-1] /= diag;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (Ordinal row=numRows-1; row > 0; --row) {
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+ind] -= (DomainScalar)val * x[j*xstride+row];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diag is last element in row
        DomainScalar diag;
        for (Ordinal row=numRows-1; row > 0; --row) {
          const Ordinal *i = inds+offs[row],
                       *ie = inds+offs[row+1];
          const Scalar  *v = vals+offs[row];
          // capture and skip the diag
          diag = v[ie-i-1];
          --ie;
          //
          for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+row] /= diag;
          while (i != ie) {
            const Ordinal ind = *i++;
            const Scalar  val = *v++;
            for (Ordinal j=0; j<numRHS; ++j) x[j*xstride+ind] -= (DomainScalar)val * x[j*xstride+row];
          }
        }
        // last row
        diag = (DomainScalar)vals[0];
        for (Ordinal j=0; j<numRHS; ++j) x[j*xstride] /= diag;
      }
    }
  };

} // namespace Kokkos

#endif
