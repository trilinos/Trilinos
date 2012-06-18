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

#ifndef KOKKOS_DEFAULTSPARSEMULTIPLY_KERNELOPS_HPP
#define KOKKOS_DEFAULTSPARSEMULTIPLY_KERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif

namespace Kokkos {

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseMultiplyOp {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t numRHS, xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t row) {
      const Scalar  *v = vals + ptrs[row];
      const Ordinal *i = inds + ptrs[row],
                   *ie = inds + ptrs[row+1];
      if (NO_BETA_AND_OVERWRITE) {
        for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] = Teuchos::ScalarTraits<RangeScalar>::zero();
      }
      else {
        for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] *= beta;
      }
      // save the extra multiplication if possible
      if (alpha == Teuchos::ScalarTraits<RangeScalar>::one()) {
        while (i != ie)
        {
          const  Scalar val = *v++;
          const Ordinal ind = *i++;
          for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] += (RangeScalar)val * (RangeScalar)x[j*xstride+ind];
        }
      }
      else { // alpha != one
        while (i != ie)
        {
          const  Scalar val = *v++;
          const Ordinal ind = *i++;
          for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] += alpha * (RangeScalar)val * (RangeScalar)x[j*xstride+ind];
        }
      }
    }
  };


  // mfh 15 June 2012: I added the ScalarIsComplex Boolean template parameter to
  // avoid build errors due to attempts to cast from Scalar=int to
  // RangeScalar=std::complex<double>.  This additional template parameter is a
  // standard technique to specialize code for the real or complex arithmetic
  // case, when attempts to write a single code for both cases would result in
  // syntax errors.  The template parameter refers to RangeScalar, because we
  // have to decide whether the RangeScalar constructor takes one or two
  // arguments.
  //
  // If you wish to optimize further, you might like to add another Boolean
  // template parameter for whether Scalar is real or complex.  First, this
  // would let you avoid calling Teuchos::ScalarTraits<Scalar>::conjugate().
  // This should inline to nothing, but not all compilers are good at inlining.
  // Second, this would simplify the code for the case when both RangeScalar and
  // Scalar are both complex.  However, this second Boolean template parameter
  // is not necessary for correctness.
  template <class Scalar, // the type of entries in the sparse matrix
            class Ordinal, // the type of indices in the sparse matrix
            class DomainScalar, // the type of entries in the input (multi)vector
            class RangeScalar, // the type of entries in the output (multi)vector
            int NO_BETA_AND_OVERWRITE,
            bool RangeScalarIsComplex=Teuchos::ScalarTraits<RangeScalar>::isComplex> // whether RangeScalar is complex valued.
  struct DefaultSparseTransposeMultiplyOp;

  // Partial specialization for when RangeScalar is real (not complex) valued.
  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp<Scalar, Ordinal, DomainScalar, RangeScalar, NO_BETA_AND_OVERWRITE, false> {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t numRHS, xstride, ystride;

    inline void execute() {
      using Teuchos::ScalarTraits;

      if (NO_BETA_AND_OVERWRITE) {
        for (size_t j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (size_t row=0; row<numCols; ++row) {
            yp[row] = ScalarTraits<RangeScalar>::zero();
          }
        }
      }
      else {
        for (size_t j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (size_t row=0; row<numCols; ++row) {
            yp[row] *= beta;
          }
        }
      }
      // save the extra multiplication if possible
      if (alpha == ScalarTraits<RangeScalar>::one()) {
        for (size_t colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + ptrs[colAt];
          const Ordinal *i  = inds + ptrs[colAt];
          const Ordinal *ie = inds + ptrs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind]
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (size_t j = 0; j < numRHS; ++j) {
              // mfh 15 June 2012: Casting Scalar to RangeScalar may produce a
              // build warning, e.g., if Scalar is int and RangeScalar is
              // double.  The way to get it to work is not to rely on casting
              // Scalar to RangeScalar.  Instead, rely on C++'s type promotion
              // rules.  For now, we just cast the input vector's value to
              // RangeScalar, to handle the common iterative refinement case of
              // an output vector with higher precision.
              y[j*ystride+ind] += val * RangeScalar (x[j*xstride+colAt]);
            }
          }
        }
      }
      else { // alpha != one
        for (size_t colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + ptrs[colAt];
          const Ordinal *i  = inds + ptrs[colAt];
          const Ordinal *ie = inds + ptrs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (size_t j=0; j<numRHS; ++j) {
              // mfh 15 June 2012: See notes above about build warnings when
              // casting val from Scalar to RangeScalar.
              y[j*ystride+ind] += alpha * val * RangeScalar (x[j*xstride+colAt]);
            }
          }
        }
      }
    }
  };

  // Partial specialization for when RangeScalar is complex valued.  The most
  // common case is that RangeScalar is a specialization of std::complex, but
  // this is not necessarily the case.  For example, CUDA has its own complex
  // arithmetic type, which is necessary because std::complex's methods are not
  // marked as device methods.  This code assumes that RangeScalar, being
  // complex valued, has a constructor which takes two arguments, each of which
  // can be converted to Teuchos::ScalarTraits<RangeScalar>::magnitudeType.
  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp<Scalar, Ordinal, DomainScalar, RangeScalar, NO_BETA_AND_OVERWRITE, true> {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t numRHS, xstride, ystride;

    inline void execute() {
      using Teuchos::ScalarTraits;
      typedef typename ScalarTraits<RangeScalar>::magnitudeType RSMT;

      if (NO_BETA_AND_OVERWRITE) {
        for (size_t j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (size_t row=0; row<numCols; ++row) {
            yp[row] = ScalarTraits<RangeScalar>::zero();
          }
        }
      }
      else {
        for (size_t j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (size_t row=0; row<numCols; ++row) {
            yp[row] *= beta;
          }
        }
      }
      // save the extra multiplication if possible
      if (alpha == ScalarTraits<RangeScalar>::one()) {
        for (size_t colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + ptrs[colAt];
          const Ordinal *i  = inds + ptrs[colAt];
          const Ordinal *ie = inds + ptrs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind]
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (size_t j = 0; j < numRHS; ++j) {
              // mfh 15 June 2012: Casting Scalar to RangeScalar via a
              // static_cast won't work if Scalar is int and RangeScalar is
              // std::complex<double>.  (This is valid code.)  This is because
              // std::complex<double> doesn't have a constructor that takes one
              // int argument.  Furthermore, the C++ standard library doesn't
              // define operator*(int, std::complex<double>), so we can't rely
              // on C++'s type promotion rules.  However, any reasonable complex
              // arithmetic type should have a two-argument constructor that
              // takes arguments convertible to RSMT
              // (ScalarTraits<RangeScalar>::magnitudeType), and we should also
              // be able to cast from Scalar to RSMT.
              //
              // The mess with taking the real and imaginary components of val
              // is because Scalar could be either real or complex, but the
              // two-argument constructor of RangeScalar expects two real
              // arguments.  You can get rid of this by adding another template
              // parameter for whether Scalar is real or complex.
              y[j*ystride+ind] += RangeScalar (RSMT (ScalarTraits<Scalar>::real (val)),
                                               RSMT (ScalarTraits<Scalar>::imag (val))) *
                RangeScalar (x[j*xstride+colAt]);
            }
          }
        }
      }
      else { // alpha != one
        for (size_t colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + ptrs[colAt];
          const Ordinal *i  = inds + ptrs[colAt];
          const Ordinal *ie = inds + ptrs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind]
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (size_t j=0; j<numRHS; ++j) {
              // mfh 15 June 2012: See notes above about it sometimes being
              // invalid to cast val from Scalar to RangeScalar.
              y[j*ystride+ind] += alpha *
                RangeScalar (RSMT (ScalarTraits<Scalar>::real (val)),
                             RSMT (ScalarTraits<Scalar>::imag (val))) *
                RangeScalar (x[j*xstride+colAt]);
            }
          }
        }
      }
    }
  };

} // namespace Kokkos

#endif
