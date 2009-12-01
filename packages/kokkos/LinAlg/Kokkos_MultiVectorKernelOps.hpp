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
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(int i) {
      return SCT::magnitude((x[i]));
    }
  };

  template <typename Scalar>
  struct WeightNormOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x, *w;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(int i) {
      Scalar tmp = x[i] / w[i];
      return SCT::real( SCT::conjugate(tmp)*tmp );
    }
  };

  template <typename Scalar>
  struct SumOp {
    const Scalar *x;
    typedef  Scalar ReductionType;
    inline static Scalar identity() {return Teuchos::ScalarTraits<Scalar>::zero();}
    Scalar reduce(Scalar x, Scalar y) {return x+y;}
    Scalar generate(int i) { return x[i]; }
  };

  template <typename Scalar>
  struct MaxAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return TEUCHOS_MAX(x,y);}
    Magnitude generate(int i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <typename Scalar>
  struct DotOp1 {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(int i) {
      Scalar xi = x[i]; 
      return SCT::real( SCT::conjugate(xi)*xi );
    }
  };

  template <typename Scalar>
  struct DotOp2 {
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    const Scalar *x, *y;
    typedef  Scalar ReductionType;
    inline static Scalar identity() {return SCT::zero();}
    Scalar reduce(Scalar x, Scalar y) {return x+y;}
    Scalar generate(int i) {
      Scalar xi = x[i]; Scalar yi = y[i];
      return SCT::real( SCT::conjugate(xi)*yi );
    }
  };

}

#endif
