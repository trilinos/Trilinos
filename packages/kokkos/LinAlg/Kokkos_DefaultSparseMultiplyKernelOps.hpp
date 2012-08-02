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

#include <Teuchos_ScalarTraits.hpp>

namespace Kokkos {

  template <class Scalar, class OffsetType, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseMultiplyOp {
    // mat data
    const OffsetType *offs;
    const Ordinal    *inds;
    const Scalar     *vals;
    // matvec params
    RangeScalar        alpha, beta;
    Ordinal numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    Ordinal numRHS, xstride, ystride;

    // CrT 27Jul12: performance improvement 2x-8x.
    // Changed the "while" loop to a "for" loop,
    // Reversed order of loops (outside is the loop over RHS vectors, inside the loop over the row),
    // Introduced temporary variables,
    // Unrolled the loop over RHS vectors by 4, and
    // Added specific cases for 1,2,3,4 RHS vectors.
    //
    // Unrolling and specialization are invisible to code that calls
    // DefaultSparseMultiplyOP, since these operations are handled via
    // an if-else within the execute function.  Without unrolling one
    // would lose about 1/3 in performance (on Intel Sandy Bridge).
    //
    // When using the Intel compiler, the code below may be vectorized
    // through the "always vectorize" directive #pragma simd (see
    // comments in the code).  The resulting performance gain on Intel
    // Sandy Bridge is 10-20%.  However, this currently only works if
    // the temporary variable tmp (of type RangeScalar) is made a
    // double/float explicitly.  This could be because the Intel
    // compiler (most recent version, not yet released) does not know
    // to ignore the SIMD reduction directive when the template
    // parameter type does not support the SIMD vector instructions.
    // This results in a compile error.  One workaround would be to
    // specialize DefaultSparseMultiplyOp for various RangeScalar
    // types, and only insert the pragma in the specialized versions.
    //
    // mfh (27 Jul 2012): Note that with the Intel compiler, #pragma
    // simd vectorizes the loop according to the "fast" floating-point
    // model (as in -fp-model=fast), even if that is not the
    // prevailing model of the rest of the compilation unit.  Sparse
    // matrix-vector multiply is a less sensitive kernel than others,
    // so this should not be a risky optimization in terms of
    // accuracy.  The main concern here is that the loop may flush
    // denormal floating-point values to zero.  For a discussion, see
    // e.g.,
    // http://software.intel.com/en-us/articles/consistency-of-floating-point-results-using-the-intel-compiler/
    //
    // mfh (27 Jul 2012): Reviewed and briefly edited changes and comments.
    inline KERNEL_PREFIX void execute(size_t row) {
      const OffsetType start = offs[row];
      const OffsetType end = offs[row+1];

      // CrT: Unroll by 4 over numRHS; specialize for numRHS <= 4.
      if (numRHS > 4) {
        Ordinal j = 0;
        if (alpha == Teuchos::ScalarTraits<RangeScalar>::one()) {
          // Strip-mined portion of the loop over RHS vectors.
          for (; j < numRHS - 3; j+=4) {
            RangeScalar* const ys = y+j*ystride + row;
            RangeScalar tmp[4] = {
              Teuchos::ScalarTraits<RangeScalar>::zero(),
              Teuchos::ScalarTraits<RangeScalar>::zero(),
              Teuchos::ScalarTraits<RangeScalar>::zero(),
              Teuchos::ScalarTraits<RangeScalar>::zero()
            };

            // CrT: Adding the pragma below improves performance by
            // 15% on Intel SandyBridge, but with the (new beta)
            // version of the Intel compiler that I tested, this
            // doesn't work when RangeScalar is a template parameter.
            // You have to change "RangeScalar tmp;" to "double tmp;".
            //
            //#pragma simd reduction(+: tmp)
            for (OffsetType i = start; i < end; i++) {
              const RangeScalar Aij = (RangeScalar)vals[i];
              const DomainScalar* const xs = x+j*xstride + inds[i];
              tmp[0] += Aij * (RangeScalar)xs[0 * xstride];
              tmp[1] += Aij * (RangeScalar)xs[1 * xstride];
              tmp[2] += Aij * (RangeScalar)xs[2 * xstride];
              tmp[3] += Aij * (RangeScalar)xs[3 * xstride];
            }
            // The compiler should remove the branch, since the test
            // value is a compile-time constant.
            if (NO_BETA_AND_OVERWRITE) {
              ys[0 * ystride] = tmp[0];
              ys[1 * ystride] = tmp[1];
              ys[2 * ystride] = tmp[2];
              ys[3 * ystride] = tmp[3];
            } else {
              ys[0 * ystride] = tmp[0] + beta * ys[0 * ystride];
              ys[1 * ystride] = tmp[1] + beta * ys[1 * ystride];
              ys[2 * ystride] = tmp[2] + beta * ys[2 * ystride];
              ys[3 * ystride] = tmp[3] + beta * ys[3 * ystride];
            }
          }
          // Clean-up portion of the loop over RHS vectors.
          for (; j < numRHS; ++j) {
            RangeScalar* const ys = y+j*ystride + row;
            const DomainScalar* const xs = x+j*xstride;
            RangeScalar tmp = Teuchos::ScalarTraits<RangeScalar>::zero();

            // CrT: Adding the pragma below improves performance by
            // 15% on Intel SandyBridge, but with the (new beta)
            // version of the Intel compiler that I tested, this
            // doesn't work when RangeScalar is a template parameter.
            // You have to change "RangeScalar tmp;" to "double tmp;".
            //
            //#pragma simd reduction(+: tmp)
            for (OffsetType i = start; i < end; i++) {
              tmp += (RangeScalar)vals[i] * (RangeScalar)xs[inds[i]];
            }
            // The compiler should remove the branch, since the test
            // value is a compile-time constant.
            if (NO_BETA_AND_OVERWRITE)
              ys[0] = tmp;
            else
              ys[0] = tmp + beta * ys[0];
          }
        }
        else { // alpha != one
          for (; j < numRHS - 3; j += 4) {
            RangeScalar* const ys = y+j*ystride + row;
            RangeScalar tmp[4] = {
              Teuchos::ScalarTraits<RangeScalar>::zero(),
              Teuchos::ScalarTraits<RangeScalar>::zero(),
              Teuchos::ScalarTraits<RangeScalar>::zero(),
              Teuchos::ScalarTraits<RangeScalar>::zero()
            };

            // CrT: +15% on SandyBridge but only if you change "RangeScalar tmp;" to "double tmp;"
            //#pragma simd reduction(+: tmp)
            for (OffsetType i = start; i < end; i++) {
              const RangeScalar Aij = (RangeScalar)vals[i];
              const DomainScalar* const xs = x+j*xstride + inds[i];
              tmp[0] += Aij * (RangeScalar)xs[0 * xstride];
              tmp[1] += Aij * (RangeScalar)xs[1 * xstride];
              tmp[2] += Aij * (RangeScalar)xs[2 * xstride];
              tmp[3] += Aij * (RangeScalar)xs[3 * xstride];
            }

            tmp[0] *= alpha;
            tmp[1] *= alpha;
            tmp[2] *= alpha;
            tmp[3] *= alpha;

            if (NO_BETA_AND_OVERWRITE) {
              ys[0 * ystride] = tmp[0];
              ys[1 * ystride] = tmp[1];
              ys[2 * ystride] = tmp[2];
              ys[3 * ystride] = tmp[3];
            } else {
              ys[0 * ystride] = tmp[0] + beta * ys[0 * ystride];
              ys[1 * ystride] = tmp[1] + beta * ys[1 * ystride];
              ys[2 * ystride] = tmp[2] + beta * ys[2 * ystride];
              ys[3 * ystride] = tmp[3] + beta * ys[3 * ystride];
            }
          }
          for (; j < numRHS; ++j) {
            RangeScalar* const ys = y+j*ystride + row;
            const DomainScalar* const xs = x+j*xstride;
            RangeScalar tmp = Teuchos::ScalarTraits<RangeScalar>::zero();

            // CrT: +15% on SandyBridge but only if you change "RangeScalar tmp;" to "double tmp;"
            //#pragma simd reduction(+: tmp)
            for (OffsetType i = start; i < end; i++) {
              tmp += (RangeScalar)vals[i] * (RangeScalar)xs[inds[i]];
            }
            tmp *= alpha;
            if (NO_BETA_AND_OVERWRITE)
              ys[0] = tmp;
            else
              ys[0] = tmp + beta * ys[0];
          }
        }
      }
      else if (numRHS == 1) {
        RangeScalar tmp = Teuchos::ScalarTraits<RangeScalar>::zero();

        //+15% on SandyBridge but only if you change "RangeScalar tmp;" to "double tmp;"
        //#pragma simd reduction(+: tmp)
        for (OffsetType i = start; i < end; i++) {
          tmp += (RangeScalar)vals[i] * (RangeScalar)x[inds[i]];
        }
        if (alpha != Teuchos::ScalarTraits<RangeScalar>::one())
          tmp *= alpha;
        if (NO_BETA_AND_OVERWRITE)
          y[row] = tmp;
        else
          y[row] = tmp + beta * y[row];
      }
      else if (numRHS == 2) {
        RangeScalar* const ys = y + row;
        RangeScalar tmp[2] = {
          Teuchos::ScalarTraits<RangeScalar>::zero(),
          Teuchos::ScalarTraits<RangeScalar>::zero()
        };

        //+15% on SandyBridge but only if you change "RangeScalar tmp;" to "double tmp;"
        //#pragma simd reduction(+: tmp)
        for (OffsetType i = start; i < end; i++) {
          const RangeScalar Aij = (RangeScalar)vals[i];
          const DomainScalar* const xs = x + inds[i];
          tmp[0] += Aij * (RangeScalar)xs[0 * xstride];
          tmp[1] += Aij * (RangeScalar)xs[1 * xstride];
        }
        if (alpha != Teuchos::ScalarTraits<RangeScalar>::one()) {
          tmp[0] *= alpha;
          tmp[1] *= alpha;
        }
        if (NO_BETA_AND_OVERWRITE) {
          ys[0 * ystride] = tmp[0];
          ys[1 * ystride] = tmp[1];
        } else {
          ys[0 * ystride] = tmp[0] + beta * ys[0 * ystride];
          ys[1 * ystride] = tmp[1] + beta * ys[1 * ystride];
        }
      }
      else if (numRHS == 3) {
        RangeScalar* const ys = y + row;
        RangeScalar tmp[3] = {
          Teuchos::ScalarTraits<RangeScalar>::zero(),
          Teuchos::ScalarTraits<RangeScalar>::zero(),
          Teuchos::ScalarTraits<RangeScalar>::zero()
        };

        //+15% on SandyBridge but only if you change "RangeScalar tmp;" to "double tmp;"
        //#pragma simd reduction(+: tmp)
        for (OffsetType i = start; i < end; i++) {
          const RangeScalar Aij = (RangeScalar)vals[i];
          const DomainScalar* const xs = x + inds[i];
          tmp[0] += Aij * (RangeScalar)xs[0 * xstride];
          tmp[1] += Aij * (RangeScalar)xs[1 * xstride];
          tmp[2] += Aij * (RangeScalar)xs[2 * xstride];
        }
        if (alpha != Teuchos::ScalarTraits<RangeScalar>::one()) {
          tmp[0] *= alpha;
          tmp[1] *= alpha;
          tmp[2] *= alpha;
        }
        if (NO_BETA_AND_OVERWRITE) {
          ys[0 * ystride] = tmp[0];
          ys[1 * ystride] = tmp[1];
          ys[2 * ystride] = tmp[2];
        } else {
          ys[0 * ystride] = tmp[0] + beta * ys[0 * ystride];
          ys[1 * ystride] = tmp[1] + beta * ys[1 * ystride];
          ys[2 * ystride] = tmp[2] + beta * ys[2 * ystride];
        }
      }
      else { // if (numRHS == 4)
        RangeScalar* const ys = y + row;
        RangeScalar tmp[4] = {
          Teuchos::ScalarTraits<RangeScalar>::zero(),
          Teuchos::ScalarTraits<RangeScalar>::zero(),
          Teuchos::ScalarTraits<RangeScalar>::zero(),
          Teuchos::ScalarTraits<RangeScalar>::zero()
        };

        //+15% on SandyBridge but only if you change "RangeScalar tmp;" to "double tmp;"
        //#pragma simd reduction(+: tmp)
        for (OffsetType i = start; i < end; i++) {
          const RangeScalar Aij = (RangeScalar)vals[i];
          const DomainScalar* const xs = x + inds[i];
          tmp[0] += Aij * (RangeScalar)xs[0 * xstride];
          tmp[1] += Aij * (RangeScalar)xs[1 * xstride];
          tmp[2] += Aij * (RangeScalar)xs[2 * xstride];
          tmp[3] += Aij * (RangeScalar)xs[3 * xstride];
        }
        if (alpha != Teuchos::ScalarTraits<RangeScalar>::one()) {
          tmp[0] *= alpha;
          tmp[1] *= alpha;
          tmp[2] *= alpha;
          tmp[3] *= alpha;
        }
        if (NO_BETA_AND_OVERWRITE) {
          ys[0 * ystride] = tmp[0];
          ys[1 * ystride] = tmp[1];
          ys[2 * ystride] = tmp[2];
          ys[3 * ystride] = tmp[3];
        } else {
          ys[0 * ystride] = tmp[0] + beta * ys[0 * ystride];
          ys[1 * ystride] = tmp[1] + beta * ys[1 * ystride];
          ys[2 * ystride] = tmp[2] + beta * ys[2 * ystride];
          ys[3 * ystride] = tmp[3] + beta * ys[3 * ystride];
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
            class OffsetType, // the type of the row offsets
            class Ordinal, // the type of indices in the sparse matrix
            class DomainScalar, // the type of entries in the input (multi)vector
            class RangeScalar, // the type of entries in the output (multi)vector
            int NO_BETA_AND_OVERWRITE,
            bool RangeScalarIsComplex=Teuchos::ScalarTraits<RangeScalar>::isComplex> // whether RangeScalar is complex valued.
  struct DefaultSparseTransposeMultiplyOp;


  // Partial specialization for when RangeScalar is real (not complex) valued.
  template <class Scalar, class OffsetType, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp<Scalar, OffsetType, Ordinal, DomainScalar, RangeScalar, NO_BETA_AND_OVERWRITE, false> {
    // mat data
    const OffsetType *offs;
    const Ordinal    *inds;
    const Scalar     *vals;
    // matvec params
    RangeScalar        alpha, beta;
    Ordinal numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    Ordinal numRHS, xstride, ystride;

    inline void execute() {
      using Teuchos::ScalarTraits;
      
      if (NO_BETA_AND_OVERWRITE) {
        for (Ordinal j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (Ordinal row=0; row<numCols; ++row) {
            yp[row] = ScalarTraits<RangeScalar>::zero();
          }
        }
      }
      else {
        for (Ordinal j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (Ordinal row=0; row<numCols; ++row) {
            yp[row] *= beta;
          }
        }
      }
      // save the extra multiplication if possible
      if (alpha == ScalarTraits<RangeScalar>::one()) {
        for (Ordinal colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + offs[colAt];
          const Ordinal *i  = inds + offs[colAt];
          const Ordinal *ie = inds + offs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind]
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (Ordinal j = 0; j < numRHS; ++j) {
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
        for (Ordinal colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + offs[colAt];
          const Ordinal *i  = inds + offs[colAt];
          const Ordinal *ie = inds + offs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (Ordinal j=0; j<numRHS; ++j) {
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
  template <class Scalar, class OffsetType, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp<Scalar, OffsetType, Ordinal, DomainScalar, RangeScalar, NO_BETA_AND_OVERWRITE, true> {
    // mat data
    const OffsetType *offs;
    const Ordinal    *inds;
    const Scalar     *vals;
    // matvec params
    RangeScalar        alpha, beta;
    Ordinal numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    Ordinal numRHS, xstride, ystride;

    inline void execute() {
      using Teuchos::ScalarTraits;
      typedef typename ScalarTraits<RangeScalar>::magnitudeType RSMT;
      if (NO_BETA_AND_OVERWRITE) {
        for (Ordinal j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (Ordinal row=0; row<numCols; ++row) {
            yp[row] = ScalarTraits<RangeScalar>::zero();
          }
        }
      }
      else {
        for (Ordinal j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (Ordinal row=0; row<numCols; ++row) {
            yp[row] *= beta;
          }
        }
      }
      // save the extra multiplication if possible
      if (alpha == ScalarTraits<RangeScalar>::one()) {
        for (Ordinal colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + offs[colAt];
          const Ordinal *i  = inds + offs[colAt];
          const Ordinal *ie = inds + offs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind]
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (Ordinal j = 0; j < numRHS; ++j) {
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
        for (Ordinal colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v  = vals + offs[colAt];
          const Ordinal *i  = inds + offs[colAt];
          const Ordinal *ie = inds + offs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind]
          while (i != ie) {
            const  Scalar val = ScalarTraits<Scalar>::conjugate (*v++);
            const Ordinal ind = *i++;
            for (Ordinal j=0; j<numRHS; ++j) {
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
