//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DEFAULTSPARSESOLVE_KERNELOPS_HPP
#define KOKKOS_DEFAULTSPARSESOLVE_KERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#include <Teuchos_ScalarTraits.hpp>

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
    const size_t  *offsets;
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
        xj[numRows-1] = yj[numRows-1];
        for (size_t r=2; r < numRows+1; ++r) {
          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
          const size_t begin = offsets[row], end = offsets[row+1];
          xj[row] = yj[row];
          for (size_t c=begin; c != end; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit
        xj[numRows-1] = yj[numRows-1] / vals[offsets[numRows-1]];
        for (size_t r=2; r < numRows+1; ++r) {
          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
          const size_t diag = offsets[row], end = offsets[row+1];
          const Scalar dval = vals[diag];
          xj[row] = yj[row];
          for (size_t c=diag+1; c != end; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
      }
      else if (!upper && unitDiag) {
        // lower + unit
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          const size_t begin = offsets[row], end = offsets[row+1];
          xj[row] = yj[row];
          for (size_t c=begin; c != end; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit
        xj[0] = yj[0] / vals[0];
        for (size_t row=1; row < numRows; ++row) {
          const size_t begin = offsets[row], diag = offsets[row+1]-1;
          const Scalar dval = vals[diag];
          xj[row] = yj[row];
          for (size_t c=begin; c != diag; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseSolveOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar      *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      const Scalar  *rowvals;
      const Ordinal *rowinds;
      Scalar dval;
      size_t nE;
      // 
      // upper triangular requires backwards substition, solving in reverse order
      // must unroll the last iteration, because decrement results in wrap-around
      // 
      if (upper && unitDiag) {
        // upper + unit
        xj[numRows-1] = yj[numRows-1];
        for (size_t row=numRows-2; row != 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          xj[row] = yj[row];
          for (size_t j=0; j != nE; ++j) {
            xj[row] -= rowvals[j] * xj[rowinds[j]];
          }
        }
        nE = numEntries[0];
        rowvals = vals_beg[0];
        rowinds = inds_beg[0];
        xj[0] = yj[0];
        for (size_t j=0; j != nE; ++j) {
          xj[0] -= rowvals[j] * xj[rowinds[j]];
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit: diagonal is first entry
        dval = vals_beg[numRows-1][0];
        xj[numRows-1] = yj[numRows-1] / dval;
        for (size_t row=numRows-2; row != 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          xj[row] = yj[row];
          Scalar dval = rowvals[0];
          for (size_t j=1; j < nE; ++j) {
            xj[row] -= rowvals[j] * xj[rowinds[j]];
          }
          xj[row] /= dval;
        }
        nE = numEntries[0];
        rowvals = vals_beg[0];
        rowinds = inds_beg[0];
        xj[0] = yj[0];
        Scalar dval = rowvals[0];
        for (size_t j=1; j < nE; ++j) {
          xj[0] -= rowvals[j] * xj[rowinds[j]];
        }
        xj[0] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          xj[row] = yj[row];
          for (size_t j=0; j < nE; ++j) {
            xj[row] -= rowvals[j] * xj[rowinds[j]];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diagonal is last entry
        nE = numEntries[0];
        rowvals = vals_beg[0];
        dval = rowvals[0];
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          dval = rowvals[nE-1];
          xj[row] = yj[row];
          if (nE > 1) {
            for (size_t j=0; j < nE-1; ++j) {
              xj[row] -= rowvals[j] * xj[rowinds[j]];
            }
          }
          xj[row] /= dval;
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeSolveOp1 {
    // mat data
    const size_t  *offsets;
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
          beg = offsets[row]; 
          endplusone = offsets[row+1];
          for (size_t j=beg; j < endplusone; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit; diag is first element in row
        size_t diag, endplusone;
        Scalar dval;
        for (size_t row=0; row < numRows-1; ++row) {
          diag = offsets[row]; 
          endplusone = offsets[row+1];
          dval = vals[diag];
          xj[row] /= dval;
          for (size_t j=diag+1; j < endplusone; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
        diag = offsets[numRows-1];
        dval = vals[diag];
        xj[numRows-1] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (size_t row=numRows-1; row > 0; --row) {
          size_t beg = offsets[row], endplusone = offsets[row+1];
          for (size_t j=beg; j < endplusone; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diag is last element in row
        Scalar dval;
        for (size_t row=numRows-1; row > 0; --row) {
          size_t beg = offsets[row], diag = offsets[row+1]-1;
          dval = vals[diag];
          xj[row] /= dval;
          for (size_t j=beg; j < diag; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
        // last row
        dval = vals[0];
        xj[0] /= dval;
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeSolveOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar      *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      const Scalar  *rowvals;
      const Ordinal *rowinds;
      Scalar dval;
      size_t nE;
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
        for (size_t row=0; row < numRows-1; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          for (size_t j=0; j < nE; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit; diag is first element in row
        for (size_t row=0; row < numRows-1; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          dval = rowvals[0];
          xj[row] /= dval;
          for (size_t j=1; j < nE; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
        rowvals = vals_beg[numRows-1];
        dval = rowvals[0];
        xj[numRows-1] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (size_t row=numRows-1; row > 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          for (size_t j=0; j < nE; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diag is last element in row
        for (size_t row=numRows-1; row > 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          dval = rowvals[nE-1];
          xj[row] /= dval;
          for (size_t j=0; j < nE-1; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
        rowvals = vals_beg[0];
        dval = rowvals[0];
        xj[0] /= dval;
      }
    }
  };

} // namespace Kokkos

#endif
