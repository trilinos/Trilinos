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
  struct ScaleOp {
    const Scalar *x;
    Scalar *y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = y[i];
      y[i] = x[i]*tmp;
    }
  };

  template <class Scalar>
  struct RecipScaleOp {
    const Scalar *x;
    Scalar *y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = y[i];
      y[i] = tmp / x[i];
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
    private:
      typedef Teuchos::ArrayRCP<const char>  CBuf;
      typedef Teuchos::ArrayRCP<      char> NCBuf;

      template <class T>
      inline static CBuf toCBuf(const Teuchos::ArrayRCP<T> &arr) {
        return Teuchos::arcp_reinterpret_cast<const char>(arr);
      }

      template <class T>
      inline static NCBuf toNCBuf(const Teuchos::ArrayRCP<T> &arr) {
        return Teuchos::arcp_reinterpret_cast<char>(arr);
      }

    public:
      //! Initialize multivector to constant value.
      inline static void Init(MultiVector<Scalar,Node> &A, Scalar alpha) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type stride = A.getStride();
        Node &node = A.getNode();
        Teuchos::ArrayRCP<Scalar> data = A.getValues();
        node.readyBuffers( Teuchos::null, Teuchos::tuple<NCBuf>(toNCBuf(data)) );
        if (stride == nR) {
          // one kernel invocation for whole multivector
          InitOp<Scalar> wdp;
          wdp.x = data.get();
          wdp.alpha = alpha;
          node.template parallel_for<InitOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          InitOp<Scalar> wdp;
          wdp.alpha = alpha;
          for (size_type j=0; j<nC; ++j) {
            wdp.x = data.get();
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
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValuesConst();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValues();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Bdata)), Teuchos::tuple<NCBuf>(toNCBuf(Adata)) );
        ScaleOp<Scalar> wdp;
        if (B.getNumCols() == 1) {
          wdp.x = Bdata.get();
          for (size_type j = 0; j < nC; ++j) {
            wdp.y = Adata.get();
            node.template parallel_for<ScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
          }
        }
        else if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.y = Adata.get();
          wdp.x = Bdata.get();
          node.template parallel_for<ScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.y = Adata.get();
            wdp.x = Bdata.get();
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
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValuesConst();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValues();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Bdata)), Teuchos::tuple<NCBuf>(toNCBuf(Adata)) );
        RecipScaleOp<Scalar> wdp;
        if (B.getNumCols() == 1) {
          // one kernel invocation for each column
          wdp.x = Bdata.get();
          for (size_type j=0; j<nC; ++j) {
            wdp.y = Adata.get();
            node.template parallel_for<RecipScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
          }
        }
        else if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector
          wdp.x = Bdata.get();
          wdp.y = Adata.get();
          node.template parallel_for<RecipScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.y = Adata.get();
            wdp.x = Bdata.get();
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
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValuesConst();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValues();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Bdata)), Teuchos::tuple<NCBuf>(toNCBuf(Adata)) );
        AssignOp<Scalar> wdp;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector assignment
          wdp.x = Adata.get();
          wdp.y = Bdata.get();
          node.template parallel_for<AssignOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.x = Adata.get();
            wdp.y = Bdata.get();
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
        TEST_FOR_EXCEPTION(nC > dots.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): dots must have length as large as number of columns of A and B.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {
            dots[j] = Teuchos::ScalarTraits<Scalar>::zero();
          }
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValuesConst(),
                               Adata = A.getValuesConst();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata),toCBuf(Bdata)), Teuchos::null );
        DotOp2<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata.get();
          op.y = Bdata.get();
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
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValuesConst(0),
                               Adata = A.getValuesConst(0);
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata),toCBuf(Bdata)), Teuchos::null );
        DotOp2<Scalar> op;
        op.x = Adata.get();
        op.y = Bdata.get();
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
        TEST_FOR_EXCEPTION(((nC != B.getNumCols()) && B.getNumCols() != 1)  
                           || nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst();
        Teuchos::ArrayRCP<Scalar>       Bdata = B.getValues();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::tuple<NCBuf>(toNCBuf(Bdata)) );
        GESUMOp<Scalar> wdp;
        wdp.alpha = alpha;
        wdp.beta  = beta;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.y = Bdata.get();
          wdp.x = Adata.get();
          node.template parallel_for<GESUMOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_type j=0; j<nC; ++j) {
            wdp.y = Bdata.get();
            wdp.x = Adata.get();
            node.template parallel_for<GESUMOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void GESUM(MultiVector<Scalar,Node> &C, Scalar alpha, const MultiVector<Scalar,Node> &A, Scalar beta, const MultiVector<Scalar,Node> &B, Scalar gamma) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Norm1(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > norms.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        SumAbsOp<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata.get();
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
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        SumAbsOp<Scalar> op;
        op.x = Adata.get();
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
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        SumOp<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata.get();
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
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst(0);
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        SumOp<Scalar> op;
        op.x = Adata.get();
        return node.parallel_reduce(0,nR,op);
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormInf(const MultiVector<Scalar,Node> &A) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst(0);
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        MaxAbsOp<Scalar> op;
        op.x = Adata.get();
        return node.parallel_reduce(0,nR,op);
      }

      inline static void NormInf(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > norms.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::NormInf(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        MaxAbsOp<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata.get();
          norms[j] = node.parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static void Norm2Squared(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_type nR = A.getNumRows();
        const size_type nC = A.getNumCols();
        const size_type Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > norms.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm2Squared(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (size_type j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst();
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        DotOp1<Scalar> op;
        for (size_type j=0; j<nC; ++j) {
          op.x = Adata.get();
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
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValuesConst(0);
        node.readyBuffers( Teuchos::tuple<CBuf>(toCBuf(Adata)), Teuchos::null );
        DotOp1<Scalar> op;
        op.x = Adata.get();
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
        Teuchos::ArrayRCP<Scalar> Adata = A.getValues();
        // we'll overwrite all data covered by the multivector, but not off-stride data
        // therefore, we are write-only only in the case that stride=nR
        bool writeOnly = (stride == nR);
        Teuchos::ArrayRCP<Scalar> mvdata = node.template viewBuffer<Scalar>(writeOnly,stride*(nC-1)+nR,Adata);
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
        TEST_FOR_EXCEPT(true);
      }

      inline static void Scale(MultiVector<Scalar,Node> &B, Scalar alpha) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Scale(MultiVector<Scalar,Node> &B, const Teuchos::ArrayView<Scalar> &alphas) {
        TEST_FOR_EXCEPT(true);
      }

  };

}

#endif
