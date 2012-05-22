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

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseSolveOp1 {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar  *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      // 
      // upper triangular requires backwards substition, solving in reverse order
      // must unroll the last iteration, because decrement results in wrap-around
      // 
      if (upper && unitDiag) {
        // upper + unit
        xj[numRows-1] = (DomainScalar)yj[numRows-1];
        for (size_t r=2; r < numRows+1; ++r) {
          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
          const size_t begin = ptrs[row], end = ptrs[row+1];
          xj[row] = (DomainScalar)yj[row];
          for (size_t c=begin; c != end; ++c) {
            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit
        xj[numRows-1] = (DomainScalar)( yj[numRows-1] / (RangeScalar)vals[ptrs[numRows-1]] );
        for (size_t r=2; r < numRows+1; ++r) {
          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
          const size_t diag = ptrs[row], end = ptrs[row+1];
          const DomainScalar dval = (DomainScalar)vals[diag];
          xj[row] = (DomainScalar)yj[row];
          for (size_t c=diag+1; c != end; ++c) {
            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
      }
      else if (!upper && unitDiag) {
        // lower + unit
        xj[0] = (DomainScalar)yj[0];
        for (size_t row=1; row < numRows; ++row) {
          const size_t begin = ptrs[row], end = ptrs[row+1];
          xj[row] = (DomainScalar)yj[row];
          for (size_t c=begin; c != end; ++c) {
            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit
        xj[0] = (DomainScalar)( yj[0] / (RangeScalar)vals[0] );
        for (size_t row=1; row < numRows; ++row) {
          const size_t begin = ptrs[row], diag = ptrs[row+1]-1;
          const DomainScalar dval = vals[diag];
          xj[row] = (DomainScalar)yj[row];
          for (size_t c=begin; c != diag; ++c) {
            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
      }
    }
  };

//  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
//  struct DefaultSparseSolveOp1 {
//    // mat data
//    const size_t  *begs;
//    const size_t  *ends;
//    const Ordinal *inds;
//    const Scalar  *vals;
//    size_t numRows;
//    // matvec params
//    bool unitDiag, upper;
//    // mv data
//    DomainScalar  *x;
//    const RangeScalar *y;
//    size_t xstride, ystride;
//
//    inline KERNEL_PREFIX void execute(size_t i) {
//      // solve rhs i for lhs i
//      const size_t rhs = i;
//      DomainScalar      *xj = x + rhs * xstride;
//      const RangeScalar *yj = y + rhs * ystride;
//      // 
//      // upper triangular requires backwards substition, solving in reverse order
//      // must unroll the last iteration, because decrement results in wrap-around
//      // 
//      if (upper && unitDiag) {
//        // upper + unit
//        xj[numRows-1] = (DomainScalar)yj[numRows-1];
//        for (size_t r=2; r < numRows+1; ++r) {
//          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
//          const size_t begin = begs[row], end = ends[row];
//          xj[row] = (DomainScalar)yj[row];
//          for (size_t c=begin; c != end; ++c) {
//            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
//          }
//        }
//      }
//      else if (upper && !unitDiag) {
//        // upper + non-unit
//        xj[numRows-1] = (DomainScalar)( yj[numRows-1] / (RangeScalar)vals[begs[numRows-1]] );
//        for (size_t r=2; r < numRows+1; ++r) {
//          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
//          const size_t diag = begs[row], end = ends[row];
//          const DomainScalar dval = (DomainScalar)vals[diag];
//          xj[row] = (DomainScalar)yj[row];
//          for (size_t c=diag+1; c != end; ++c) {
//            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
//          }
//          xj[row] /= dval;
//        }
//      }
//      else if (!upper && unitDiag) {
//        // lower + unit
//        xj[0] = (DomainScalar)yj[0];
//        for (size_t row=1; row < numRows; ++row) {
//          const size_t begin = begs[row], end = ends[row];
//          xj[row] = (DomainScalar)yj[row];
//          for (size_t c=begin; c != end; ++c) {
//            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
//          }
//        }
//      }
//      else if (!upper && !unitDiag) {
//        // lower + non-unit
//        xj[0] = (DomainScalar)( yj[0] / (RangeScalar)vals[0] );
//        for (size_t row=1; row < numRows; ++row) {
//          const size_t begin = begs[row], diag = ends[row]-1;
//          const DomainScalar dval = vals[diag];
//          xj[row] = (DomainScalar)yj[row];
//          for (size_t c=begin; c != diag; ++c) {
//            xj[row] -= (DomainScalar)vals[c] * xj[inds[c]];
//          }
//          xj[row] /= dval;
//        }
//      }
//    }
//  };
//
//
//  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
//  struct DefaultSparseSolveOp2 {
//    // mat data
//    const Ordinal * const * inds_beg;
//    const Scalar  * const * vals_beg;
//    const size_t  *         numEntries;
//    size_t numRows;
//    // matvec params
//    bool unitDiag, upper;
//    // mv data
//    DomainScalar      *x;
//    const RangeScalar *y;
//    size_t xstride, ystride;
//
//    inline KERNEL_PREFIX void execute(size_t i) {
//      // solve rhs i for lhs i
//      const size_t rhs = i;
//      DomainScalar      *xj = x + rhs * xstride;
//      const RangeScalar *yj = y + rhs * ystride;
//      const Scalar  *rowvals;
//      const Ordinal *rowinds;
//      DomainScalar dval;
//      size_t nE;
//      // 
//      // upper triangular requires backwards substition, solving in reverse order
//      // must unroll the last iteration, because decrement results in wrap-around
//      // 
//      if (upper && unitDiag) {
//        // upper + unit
//        xj[numRows-1] = (DomainScalar)yj[numRows-1];
//        for (size_t row=numRows-2; row != 0; --row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          xj[row] = yj[row];
//          for (size_t j=0; j != nE; ++j) {
//            xj[row] -= (DomainScalar)rowvals[j] * xj[rowinds[j]];
//          }
//        }
//        nE = numEntries[0];
//        rowvals = vals_beg[0];
//        rowinds = inds_beg[0];
//        xj[0] = (DomainScalar)yj[0];
//        for (size_t j=0; j != nE; ++j) {
//          xj[0] -= (DomainScalar)rowvals[j] * xj[rowinds[j]];
//        }
//      }
//      else if (upper && !unitDiag) {
//        // upper + non-unit: diagonal is first entry
//        dval = (DomainScalar)vals_beg[numRows-1][0];
//        xj[numRows-1] = (DomainScalar)yj[numRows-1] / dval;
//        for (size_t row=numRows-2; row != 0; --row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          xj[row] = (DomainScalar)yj[row];
//          Scalar dval_inner = rowvals[0];
//          for (size_t j=1; j < nE; ++j) {
//            xj[row] -= (DomainScalar)rowvals[j] * xj[rowinds[j]];
//          }
//          xj[row] /= dval_inner;
//        }
//        nE = numEntries[0];
//        rowvals = vals_beg[0];
//        rowinds = inds_beg[0];
//        xj[0] = (DomainScalar)yj[0];
//        DomainScalar dval_inner = (DomainScalar)rowvals[0];
//        for (size_t j=1; j < nE; ++j) {
//          xj[0] -= (DomainScalar)rowvals[j] * xj[rowinds[j]];
//        }
//        xj[0] /= dval_inner;
//      }
//      else if (!upper && unitDiag) {
//        // lower + unit
//        xj[0] = (DomainScalar)yj[0];
//        for (size_t row=1; row < numRows; ++row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          xj[row] = (DomainScalar)yj[row];
//          for (size_t j=0; j < nE; ++j) {
//            xj[row] -= (DomainScalar)rowvals[j] * xj[rowinds[j]];
//          }
//        }
//      }
//      else if (!upper && !unitDiag) {
//        // lower + non-unit; diagonal is last entry
//        nE = numEntries[0];
//        rowvals = vals_beg[0];
//        dval = (DomainScalar)rowvals[0];
//        xj[0] = yj[0];
//        for (size_t row=1; row < numRows; ++row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          dval = (DomainScalar)rowvals[nE-1];
//          xj[row] = (DomainScalar)yj[row];
//          if (nE > 1) {
//            for (size_t j=0; j < nE-1; ++j) {
//              xj[row] -= (DomainScalar)rowvals[j] * xj[rowinds[j]];
//            }
//          }
//          xj[row] /= dval;
//        }
//      }
//    }
//  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeSolveOp1 {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar  *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      // 
      // put y into x and solve system in-situ
      // this is copy-safe, in the scenario that x and y point to the same location.
      //
      for (size_t row=0; row < numRows; ++row) {
        xj[row] = yj[row];
      }
      // 
      if (upper && unitDiag) {
        // upper + unit
        size_t beg, endplusone;
        for (size_t row=0; row < numRows-1; ++row) {
          beg = ptrs[row]; 
          endplusone = ptrs[row+1];
          for (size_t j=beg; j < endplusone; ++j) {
            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit; diag is first element in row
        size_t diag, endplusone;
        DomainScalar dval;
        for (size_t row=0; row < numRows-1; ++row) {
          diag = ptrs[row]; 
          endplusone = ptrs[row+1];
          dval = (DomainScalar)vals[diag];
          xj[row] /= dval;
          for (size_t j=diag+1; j < endplusone; ++j) {
            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
          }
        }
        diag = ptrs[numRows-1];
        dval = (DomainScalar)vals[diag];
        xj[numRows-1] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (size_t row=numRows-1; row > 0; --row) {
          size_t beg = ptrs[row], endplusone = ptrs[row+1];
          for (size_t j=beg; j < endplusone; ++j) {
            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diag is last element in row
        DomainScalar dval;
        for (size_t row=numRows-1; row > 0; --row) {
          size_t beg = ptrs[row], diag = ptrs[row+1]-1;
          dval = (DomainScalar)vals[diag];
          xj[row] /= dval;
          for (size_t j=beg; j < diag; ++j) {
            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
          }
        }
        // last row
        dval = (DomainScalar)vals[0];
        xj[0] /= dval;
      }
    }
  };

//  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
//  struct DefaultSparseTransposeSolveOp1 {
//    // mat data
//    const size_t  *begs;
//    const size_t  *ends;
//    const Ordinal *inds;
//    const Scalar  *vals;
//    size_t numRows;
//    // matvec params
//    bool unitDiag, upper;
//    // mv data
//    DomainScalar  *x;
//    const RangeScalar *y;
//    size_t xstride, ystride;
//
//    inline KERNEL_PREFIX void execute(size_t i) {
//      // solve rhs i for lhs i
//      const size_t rhs = i;
//      DomainScalar      *xj = x + rhs * xstride;
//      const RangeScalar *yj = y + rhs * ystride;
//      // 
//      // put y into x and solve system in-situ
//      // this is copy-safe, in the scenario that x and y point to the same location.
//      //
//      for (size_t row=0; row < numRows; ++row) {
//        xj[row] = yj[row];
//      }
//      // 
//      if (upper && unitDiag) {
//        // upper + unit
//        size_t beg, endplusone;
//        for (size_t row=0; row < numRows-1; ++row) {
//          beg = begs[row]; 
//          endplusone = ends[row];
//          for (size_t j=beg; j < endplusone; ++j) {
//            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
//          }
//        }
//      }
//      else if (upper && !unitDiag) {
//        // upper + non-unit; diag is first element in row
//        size_t diag, endplusone;
//        DomainScalar dval;
//        for (size_t row=0; row < numRows-1; ++row) {
//          diag = begs[row]; 
//          endplusone = ends[row];
//          dval = (DomainScalar)vals[diag];
//          xj[row] /= dval;
//          for (size_t j=diag+1; j < endplusone; ++j) {
//            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
//          }
//        }
//        diag = begs[numRows-1];
//        dval = (DomainScalar)vals[diag];
//        xj[numRows-1] /= dval;
//      }
//      else if (!upper && unitDiag) {
//        // lower + unit
//        for (size_t row=numRows-1; row > 0; --row) {
//          size_t beg = begs[row], endplusone = ends[row];
//          for (size_t j=beg; j < endplusone; ++j) {
//            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
//          }
//        }
//      }
//      else if (!upper && !unitDiag) {
//        // lower + non-unit; diag is last element in row
//        DomainScalar dval;
//        for (size_t row=numRows-1; row > 0; --row) {
//          size_t beg = begs[row], diag = ends[row]-1;
//          dval = (DomainScalar)vals[diag];
//          xj[row] /= dval;
//          for (size_t j=beg; j < diag; ++j) {
//            xj[inds[j]] -= (DomainScalar)vals[j] * xj[row];
//          }
//        }
//        // last row
//        dval = (DomainScalar)vals[0];
//        xj[0] /= dval;
//      }
//    }
//  };
//
//
//  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
//  struct DefaultSparseTransposeSolveOp2 {
//    // mat data
//    const Ordinal * const * inds_beg;
//    const Scalar  * const * vals_beg;
//    const size_t  *         numEntries;
//    size_t numRows;
//    // matvec params
//    bool unitDiag, upper;
//    // mv data
//    DomainScalar      *x;
//    const RangeScalar *y;
//    size_t xstride, ystride;
//
//    inline KERNEL_PREFIX void execute(size_t i) {
//      // solve rhs i for lhs i
//      const size_t rhs = i;
//      DomainScalar      *xj = x + rhs * xstride;
//      const RangeScalar *yj = y + rhs * ystride;
//      const Scalar  *rowvals;
//      const Ordinal *rowinds;
//      DomainScalar dval;
//      size_t nE;
//      // 
//      // put y into x and solve system in-situ
//      // this is copy-safe, in the scenario that x and y point to the same location.
//      //
//      for (size_t row=0; row < numRows; ++row) {
//        xj[row] = (DomainScalar)yj[row];
//      }
//      // 
//      if (upper && unitDiag) {
//        // upper + unit
//        for (size_t row=0; row < numRows-1; ++row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          for (size_t j=0; j < nE; ++j) {
//            xj[rowinds[j]] -= (DomainScalar)rowvals[j] * xj[row];
//          }
//        }
//      }
//      else if (upper && !unitDiag) {
//        // upper + non-unit; diag is first element in row
//        for (size_t row=0; row < numRows-1; ++row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          dval = (DomainScalar)rowvals[0];
//          xj[row] /= dval;
//          for (size_t j=1; j < nE; ++j) {
//            xj[rowinds[j]] -= (DomainScalar)rowvals[j] * xj[row];
//          }
//        }
//        rowvals = vals_beg[numRows-1];
//        dval = (DomainScalar)rowvals[0];
//        xj[numRows-1] /= dval;
//      }
//      else if (!upper && unitDiag) {
//        // lower + unit
//        for (size_t row=numRows-1; row > 0; --row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          for (size_t j=0; j < nE; ++j) {
//            xj[rowinds[j]] -= (DomainScalar)rowvals[j] * xj[row];
//          }
//        }
//      }
//      else if (!upper && !unitDiag) {
//        // lower + non-unit; diag is last element in row
//        for (size_t row=numRows-1; row > 0; --row) {
//          nE = numEntries[row];
//          rowvals = vals_beg[row];
//          rowinds = inds_beg[row];
//          dval = (DomainScalar)rowvals[nE-1];
//          xj[row] /= dval;
//          for (size_t j=0; j < nE-1; ++j) {
//            xj[rowinds[j]] -= (DomainScalar)rowvals[j] * xj[row];
//          }
//        }
//        rowvals = vals_beg[0];
//        dval = (DomainScalar)rowvals[0];
//        xj[0] /= dval;
//      }
//    }
//  };

} // namespace Kokkos

#endif
