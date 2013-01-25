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

#ifndef KOKKOS_MULTIVECTOR_KERNELOPS_HPP
#define KOKKOS_MULTIVECTOR_KERNELOPS_HPP

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif


#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <typename Scalar>
  struct InitOp {
    Scalar *x;
    Scalar alpha;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = alpha;
    }
  };

  template <typename Scalar>
  struct AssignOp {
    Scalar *x;
    const Scalar *y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = y[i];
    }
  };

  template <typename Scalar>
  struct SingleScaleOp {
    Scalar alpha;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = x[i];
      x[i] = alpha*tmp;
    }
  };

  template <typename Scalar>
  struct MVScaleOp {
    Scalar alpha;
    const Scalar *y;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = alpha*y[i];
    }
  };

  template <typename Scalar>
  struct MVElemMultOp {
    Scalar scalarX;
    Scalar scalarYZ;
    const Scalar *y;
    const Scalar *z;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = scalarX*x[i] + scalarYZ*y[i]*z[i];
    }
  };

  template <typename Scalar>
  struct AbsOp {
    const Scalar *y;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = Teuchos::ScalarTraits<Scalar>::magnitude(y[i]);
    }
  };

  template <typename Scalar>
  struct RecipOp {
    const Scalar *scale;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = Teuchos::ScalarTraits<Scalar>::one() / scale[i];
    }
  };

  /// \class ReciprocalThresholdOp
  /// \brief Compute x[i] = 1.0 / x[i] when the magnitude of x[i] is
  ///   big enough, else replace x[i] with a given value.
  template <typename Scalar>
  struct ReciprocalThresholdOp {
    typedef Scalar scalar_type;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    scalar_type* const x_;
    const scalar_type minDiagVal_;
    const magnitude_type minDiagMag_;

    /// \brief Constructor.
    ///
    /// \param x [in/out] The vector of entries to invert.
    /// \param minDiagVal [in] The value to use to replace x[i], if
    ///   x[i] is smaller in magnitude than the magnitude of
    ///   minDiagVal.
    ReciprocalThresholdOp (scalar_type* const x,
			   const scalar_type& minDiagVal) :
      x_ (x),
      minDiagVal_ (minDiagVal),
      minDiagMag_ (STS::magnitude (minDiagVal))
    {}

    inline KERNEL_PREFIX void execute (int i) const {
      if (STS::magnitude (x_[i]) < minDiagMag_) {
	x_[i] = minDiagVal_;
      }
      else {
	x_[i] = STS::one() / x_[i];
      }
    }
  };

  template <typename Scalar>
  struct GESUMOp {
    const Scalar *x;
    Scalar *y;
    Scalar alpha, beta;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = y[i];
      y[i] = alpha * x[i] + beta * tmp;
    }
  };

  template <typename Scalar>
  struct GESUMOp3 {
    const Scalar *x, *y;
    Scalar *z;
    Scalar alpha, beta, gamma;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = z[i];
      z[i] = alpha * x[i] + beta * y[i] + gamma * tmp;
    }
  };

  template <typename Scalar>
  struct SumAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude KERNEL_PREFIX identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    inline static Magnitude KERNEL_PREFIX reduce(Magnitude x, Magnitude y) {return x+y;}
    inline        Magnitude KERNEL_PREFIX generate(int i) {
      return SCT::magnitude((x[i]));
    }
  };

  template <typename Scalar>
  struct WeightNormOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x, *w;
    typedef  Magnitude ReductionType;
    inline static Magnitude KERNEL_PREFIX identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    inline static Magnitude KERNEL_PREFIX reduce(Magnitude x, Magnitude y) {return x+y;}
    inline        Magnitude KERNEL_PREFIX generate(int i) {
      Scalar tmp = x[i] / w[i];
      return SCT::real( SCT::conjugate(tmp)*tmp );
    }
  };

  template <typename Scalar>
  struct SumOp {
    const Scalar *x;
    typedef  Scalar ReductionType;
    inline static Scalar KERNEL_PREFIX identity() {return Teuchos::ScalarTraits<Scalar>::zero();}
    inline static Scalar KERNEL_PREFIX reduce(Scalar x, Scalar y) {return x+y;}
    inline        Scalar KERNEL_PREFIX generate(int i) { return x[i]; }
  };

  template <typename Scalar>
  struct MaxAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude KERNEL_PREFIX identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    inline static Magnitude KERNEL_PREFIX reduce(Magnitude x, Magnitude y) {return TEUCHOS_MAX(x,y);}
    inline        Magnitude KERNEL_PREFIX generate(int i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <typename Scalar>
  struct DotOp1 {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude KERNEL_PREFIX identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    inline static Magnitude KERNEL_PREFIX reduce(Magnitude x, Magnitude y) {return x+y;}
    inline        Magnitude KERNEL_PREFIX generate(int i) {
      Scalar xi = x[i]; 
      return SCT::real( SCT::conjugate(xi)*xi );
    }
  };

  template <typename Scalar>
  struct DotOp2 {
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    const Scalar *x, *y;
    typedef  Scalar ReductionType;
    inline static Scalar KERNEL_PREFIX identity() {return SCT::zero();}
    inline static Scalar KERNEL_PREFIX reduce(Scalar x, Scalar y) {return x+y;}
    inline        Scalar KERNEL_PREFIX generate(int i) {
      Scalar xi = x[i]; Scalar yi = y[i];
      return SCT::conjugate(xi)*yi;
    }
  };

}

#endif
