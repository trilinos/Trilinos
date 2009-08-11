//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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

#ifndef KOKKOS_DEFAULTARITHMETIC_H
#define KOKKOS_DEFAULTARITHMETIC_H

#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <stdexcept>

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar>
  struct InitOp {
    Scalar *x;
    Scalar alpha;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = alpha;
    }
  };

  template <class Scalar>
  struct AssignOp {
    Scalar *x;
    const Scalar *y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = y[i];
    }
  };

  template <class Scalar>
  struct SingleScaleOp {
    Scalar alpha;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = x[i];
      x[i] = alpha*tmp;
    }
  };

  template <class Scalar>
  struct MVScaleOp {
    Scalar alpha;
    const Scalar *y;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = alpha*y[i];
    }
  };

  template <class Scalar>
  struct ScaleOp {
    const Scalar *scale;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = x[i];
      x[i] = scale[i]*tmp;
    }
  };

  template <class Scalar>
  struct RecipScaleOp {
    const Scalar *scale;
    Scalar *x;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = x[i];
      x[i] = tmp / scale[i];
    }
  };

  template <class Scalar>
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

  template <class Scalar>
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

  template <class Scalar>
  struct SumAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(int i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <class Scalar>
  struct SumOp {
    const Scalar *x;
    typedef  Scalar ReductionType;
    inline static Scalar identity() {return Teuchos::ScalarTraits<Scalar>::zero();}
    Scalar reduce(Scalar x, Scalar y) {return x+y;}
    Scalar generate(int i) { return x[i]; }
  };

  template <class Scalar>
  struct MaxAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return std::max(x,y);}
    Magnitude generate(int i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <class Scalar>
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

  template <class Scalar>
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

  //REFACTOR// //! Class for providing GEMM for a particular Node
  //REFACTOR// template <class Scalar, class Node> 
  //REFACTOR// GEMM

  //! Class DefaultArithmetic, unimplemented
  template <class MV>
  class DefaultArithmetic {
    public:
      //! Multiply one MultiVector by another, element-wise: B *= A
      static void Multiply(const MV &A, MV &B) { 
        TEST_FOR_EXCEPTION(true,std::logic_error,"DefaultArithmetic<" << Teuchos::typeName(A) << ">: no specialization exists for given multivector type.");
      }

      //! Divide one MultiVector by another, element-wise: B /= A
      static void Divide(MV &A, const MV &B) {
        TEST_FOR_EXCEPTION(true,std::logic_error,"DefaultArithmetic<" << Teuchos::typeName(A) << ": no specialization exists for given multivector type.");
      }
  };

  //! class DefaultArithmetic, for Kokkos::MultiVector
  template <class Scalar, class Node>
  class DefaultArithmetic<MultiVector<Scalar,Node> > {

    public:
      //! Initialize multivector to constant value.
      inline static void Init(MultiVector<Scalar,Node> &A, Scalar alpha) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type stride = A.getStride();
        Node &node = A.getNode();
        Teuchos::ArrayRCP<Scalar> data = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addNonConstBuffer<Scalar>(data);
        rbh.end();
        // prepare op
        InitOp<Scalar> wdp;
        wdp.alpha = alpha;
        if (stride == nR) {
          // one kernel invocation for whole multivector
          wdp.x = data(0,nR*nC).getRawPtr();
          node.template parallel_for<InitOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.x = data(0,nR).getRawPtr();
            node.template parallel_for<InitOp<Scalar> >(0,nR,wdp);
            data += stride;
          }
        }
      }

      //! Multiply one MultiVector by another, element-wise: B *= A
      inline static void Multiply(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        TEST_FOR_EXCEPTION(((nC != B.getNumCols()) && B.getNumCols() != 1)  
                            || nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Multiply(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Adata);
        rbh.end();
        // prepare op
        ScaleOp<Scalar> wdp;
        if (B.getNumCols() == 1) {
          wdp.scale = Bdata(0,nR).getRawPtr();
          for (size_type j = 0; j < nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            node.template parallel_for<ScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
          }
        }
        else if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.scale = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node.template parallel_for<ScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.scale = Bdata(0,nR).getRawPtr();
            wdp.x = Adata(0,nR).getRawPtr();
            node.template parallel_for<ScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      //! Divide one MultiVector by another, element-wise: B /= A
      inline static void Divide(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        TEST_FOR_EXCEPTION(((nC != B.getNumCols()) && B.getNumCols() != 1)  
                           || nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Adata);
        rbh.end();
        RecipScaleOp<Scalar> wdp;
        if (B.getNumCols() == 1) {
          // one kernel invocation for each column
          wdp.scale = Bdata(0,nR).getRawPtr();
          for (size_type j=0; j<nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            node.template parallel_for<RecipScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
          }
        }
        else if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector
          wdp.scale = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node.template parallel_for<RecipScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            wdp.scale = Bdata(0,nR).getRawPtr();
            node.template parallel_for<RecipScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      //! Assign one MultiVector to another
      inline static void Assign(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        if (nC*nR == 0) return;
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Adata);
        rbh.end();
        // prepare op
        AssignOp<Scalar> wdp;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector assignment
          wdp.x = Adata(0,nR*nC).getRawPtr();
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          node.template parallel_for<AssignOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            wdp.y = Bdata(0,nR).getRawPtr();
            node.template parallel_for<AssignOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void Dot(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, const Teuchos::ArrayView<Scalar> &dots) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same dimensions.");
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_type>(dots.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): dots must have length as large as number of columns of A and B.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {
            dots[j] = Teuchos::ScalarTraits<Scalar>::zero();
          }
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues(),
                                        Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp2<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          op.y = Bdata(0,nR).getRawPtr();
          dots[j] = node.parallel_reduce(0,nR,op);
          Adata += Astride;
          Bdata += Bstride;
        }
      }

      inline static Scalar Dot(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same number of rows.");
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues(0),
                                        Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp2<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        op.y = Bdata(0,nR).getRawPtr();
        return node.parallel_reduce(0,nR,op);
      }

      inline static void GEMM(MultiVector<Scalar,Node> &C, Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, Scalar beta) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void GESUM(MultiVector<Scalar,Node> &B, Scalar alpha, const MultiVector<Scalar,Node> &A, Scalar beta) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        Teuchos::ArrayRCP<Scalar>       Bdata = B.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addNonConstBuffer<Scalar>(Bdata);
        rbh.end();
        GESUMOp<Scalar> wdp;
        wdp.alpha = alpha;
        wdp.beta  = beta;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node.template parallel_for<GESUMOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.y = Bdata(0,nR).getRawPtr();
            wdp.x = Adata(0,nR).getRawPtr();
            node.template parallel_for<GESUMOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void GESUM(MultiVector<Scalar,Node> &C, Scalar alpha, const MultiVector<Scalar,Node> &A, Scalar beta, const MultiVector<Scalar,Node> &B, Scalar gamma) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        const size_type Cstride = C.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(),
                                        Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Cdata = C.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Cdata);
        rbh.end();
        GESUMOp3<Scalar> wdp;
        wdp.alpha = alpha;
        wdp.beta  = beta;
        wdp.gamma = gamma;
        if (Astride == nR && Bstride == nR && Cstride == nR) {
          // one kernel invocation for whole multivector
          wdp.z = Cdata(0,nR*nC).getRawPtr();
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node.template parallel_for<GESUMOp3<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.z = Cdata(0,nR).getRawPtr();
            wdp.y = Bdata(0,nR).getRawPtr();
            wdp.x = Adata(0,nR).getRawPtr();
            node.template parallel_for<GESUMOp3<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
            Cdata += Cstride;
          }
        }
      }

      inline static void Norm1(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_type>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumAbsOp<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          norms[j] = node.parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType
      Norm1(const MultiVector<Scalar,Node> &A) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesNonConst(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumAbsOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node.parallel_reduce(0,nR,op);
      }

      inline static void Sum(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<Scalar> &sums) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > sums.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Sum(A,sums): sums must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_type j=0; j<nC; ++j) {sums[j] = zero;}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumOp<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          sums[j] = node.parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static Scalar Sum(const MultiVector<Scalar,Node> &A) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node.parallel_reduce(0,nR,op);
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormInf(const MultiVector<Scalar,Node> &A) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        MaxAbsOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node.parallel_reduce(0,nR,op);
      }

      inline static void NormInf(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_type>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::NormInf(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        MaxAbsOp<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          norms[j] = node.parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static void Norm2Squared(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_type>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm2Squared(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp1<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          norms[j] = node.parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
      Norm2Squared(const MultiVector<Scalar,Node> &A) {
        const size_type nR = A.getNumRows();
        if (nR == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp1<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node.parallel_reduce(0,nR,op);
      }

      inline static void WeightedNorm(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &weightVector, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Random(MultiVector<Scalar,Node> &A) {
        // TODO: consider adding rand() functionality to node
        // in the meantime, just generate random numbers via Teuchos and then copy to node
        typedef Teuchos::ScalarTraits<Scalar> SCT;
        const size_type stride = A.getStride();
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        if (nR*nC == 0) return;
        Node &node = A.getNode();
        Teuchos::ArrayRCP<Scalar> Adata = A.getValuesNonConst();
        // we'll overwrite all data covered by the multivector, but not off-stride data
        // therefore, we are write-only only in the case that stride=nR
        bool writeOnly = (stride == nR);
        Teuchos::ArrayRCP<Scalar> mvdata = node.template viewBufferNonConst<Scalar>(writeOnly,stride*(nC-1)+nR,Adata);
        for (size_type j=0; j<nC; ++j) {
          for (size_type i=0; i<nR; ++i) {
            mvdata[j*stride + i] = SCT::random();
          }
        }
        mvdata = Teuchos::null;
      }

      inline static void Abs(MultiVector<Scalar,Node> &B, const MultiVector<Scalar,Node> &A) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Scale(MultiVector<Scalar,Node> &B, Scalar alpha, const MultiVector<Scalar,Node> &A) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        const size_type Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        Teuchos::ArrayRCP<Scalar>       Bdata = B.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addNonConstBuffer<Scalar>(Bdata);
        rbh.end();
        MVScaleOp<Scalar> wdp;
        wdp.alpha = alpha;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.x = Bdata(0,nR*nC).getRawPtr();
          wdp.y = Adata(0,nR*nC).getRawPtr();
          node.template parallel_for<MVScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.x = Bdata(0,nR).getRawPtr();
            wdp.y = Adata(0,nR).getRawPtr();
            node.template parallel_for<MVScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void Scale(MultiVector<Scalar,Node> &A, Scalar alpha) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type stride = A.getStride();
        Node &node = A.getNode();
        Teuchos::ArrayRCP<Scalar> data = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addNonConstBuffer<Scalar>(data);
        rbh.end();
        // prepare op
        SingleScaleOp<Scalar> wdp;
        wdp.alpha = alpha;
        if (stride == nR) {
          // one kernel invocation for whole multivector
          wdp.x = data(0,nR*nC).getRawPtr();
          node.template parallel_for<SingleScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.x = data(0,nR).getRawPtr();
            node.template parallel_for<SingleScaleOp<Scalar> >(0,nR,wdp);
            data += stride;
          }
        }
      }

      inline static void Scale(MultiVector<Scalar,Node> &B, const Teuchos::ArrayView<Scalar> &alphas) {
        TEST_FOR_EXCEPT(true);
      }

  };

}

#endif
