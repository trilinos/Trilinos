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
#include <stdexcept>

#include "Kokkos_MultiVector.hpp"

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar, class Node>
  struct InitOp {
    Teuchos::ArrayRCP<Scalar> x;
    Scalar alpha;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = alpha;
    }
  };

  template <class Scalar, class Node>
  struct AssignOp {
    Teuchos::ArrayRCP<Scalar> x;
    Teuchos::ArrayRCP<const Scalar> y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      x[i] = y[i];
    }
  };

  template <class Scalar, class Node>
  struct ScaleOp {
    Teuchos::ArrayRCP<const Scalar> x;
    Teuchos::ArrayRCP<Scalar> y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = y[i];
      y[i] = x[i]*tmp;
    }
  };

  template <class Scalar, class Node>
  struct RecipScaleOp {
    Teuchos::ArrayRCP<const Scalar> x;
    Teuchos::ArrayRCP<Scalar> y;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = y[i];
      y[i] = tmp / x[i];
    }
  };

  template <class Scalar, class Node>
  struct GESUMOp {
    Teuchos::ArrayRCP<const Scalar> x;
    Teuchos::ArrayRCP<Scalar> y;
    Scalar alpha, beta;
    inline KERNEL_PREFIX void execute(int i) const
    {
      Scalar tmp = y[i];
      y[i] = alpha * x[i] + beta * tmp;
    }
  };

  template <class Scalar, class Node>
  struct SumAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    Teuchos::ArrayRCP<const Scalar> x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(int i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <class Scalar, class Node>
  struct SumOp {
    Teuchos::ArrayRCP<const Scalar> x;
    typedef  Scalar ReductionType;
    inline static Scalar identity() {return Teuchos::ScalarTraits<Scalar>::zero();}
    Scalar reduce(Scalar x, Scalar y) {return x+y;}
    Scalar generate(int i) { return x[i]; }
  };

  template <class Scalar, class Node>
  struct MaxAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    Teuchos::ArrayRCP<const Scalar> x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return std::max(x,y);}
    Magnitude generate(int i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <class Scalar, class Node>
  struct DotOp1 {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    Teuchos::ArrayRCP<const Scalar> x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(int i) {
      Scalar xi = x[i]; 
      return SCT::real( SCT::conjugate(xi)*xi );
    }
  };

  template <class Scalar, class Node>
  struct DotOp2 {
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    Teuchos::ArrayRCP<const Scalar> x, y;
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
  template <class Scalar, class Ordinal, class Node>
  class DefaultArithmetic<MultiVector<Scalar,Ordinal,Node> > {
    public:

      //! Initialize multivector to constant value.
      inline static void Init(MultiVector<Scalar,Ordinal,Node> &A, Scalar alpha) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        Node &node = A.getNode();
        if (A.getStride() == nR) {
          // one kernel invocation for whole multivector
          InitOp<Scalar,Node> wdp;
          wdp.x = A.getValues(0);
          wdp.alpha = alpha;
          node.template parallel_for<InitOp<Scalar,Node> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          InitOp<Scalar,Node> wdp;
          wdp.alpha = alpha;
          for (Ordinal j=0; j<nC; ++j) {
            wdp.x = A.getValues(j);
            node.template parallel_for<InitOp<Scalar,Node> >(0,nR,wdp);
          }
        }
      }

      //! Multiply one MultiVector by another, element-wise: B *= A
      inline static void Multiply(MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &B) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(((nC != B.getNumCols()) && B.getNumCols() != 1)  ||
                           nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Multiply(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        if (B.getNumCols() == 1) {
          ScaleOp<Scalar,Node> wdp;
          wdp.x = B.getValues(0);
          for (Ordinal j=0; j<nC; ++j) {
            wdp.y = A.getValues(j);
            node.template parallel_for<ScaleOp<Scalar,Node> >(0,nR,wdp);
          }
        }
        else if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector
          ScaleOp<Scalar,Node> wdp;
          wdp.y = A.getValues(0);
          wdp.x = B.getValues(0);
          node.template parallel_for<ScaleOp<Scalar,Node> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          ScaleOp<Scalar,Node> wdp;
          for (Ordinal j=0; j<nC; ++j) {
            wdp.y = A.getValues(j);
            wdp.x = B.getValues(j);
            node.template parallel_for<ScaleOp<Scalar,Node> >(0,nR,wdp);
          }
        }
      }

      //! Divide one MultiVector by another, element-wise: B /= A
      inline static void Divide(MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &B) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(((nC != B.getNumCols()) && B.getNumCols() != 1)  ||
                           nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        RecipScaleOp<Scalar,Node> wdp;
        if (B.getNumCols() == 1) {
          // one kernel invocation for each column
          wdp.x = B.getValues(0);
          for (Ordinal j=0; j<nC; ++j) {
            wdp.y = A.getValues(j);
            node.template parallel_for<RecipScaleOp<Scalar,Node> >(0,nR,wdp);
          }
        }
        else if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector
          wdp.x = B.getValues(0);
          wdp.y = A.getValues(0);
          node.template parallel_for<RecipScaleOp<Scalar,Node> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (Ordinal j=0; j<nC; ++j) {
            wdp.y = A.getValues(j);
            wdp.x = B.getValues(j);
            node.template parallel_for<RecipScaleOp<Scalar,Node> >(0,nR,wdp);
          }
        }
      }

      //! Assign one MultiVector to another
      inline static void Assign(MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &B) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() ||
                           nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        if (nC*nR == 0) return;
        Node &node = A.getNode();
        AssignOp<Scalar,Node> wdp;
        if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector assignment
          wdp.x = A.getValues(0);
          wdp.y = B.getValues(0);
          node.template parallel_for<AssignOp<Scalar,Node> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (Ordinal j=0; j<nC; ++j) {
            wdp.x = A.getValues(j);
            wdp.y = B.getValues(j);
            node.template parallel_for<AssignOp<Scalar,Node> >(0,nR,wdp);
          }
        }
      }

      inline static void Dot(const MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &B, const Teuchos::ArrayView<Scalar> &dots) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() ||
                           nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same dimensions.");
        TEST_FOR_EXCEPTION(nC > dots.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): dots must have length as large as number of columns of A and B.");
        if (nR*nC == 0) {
          for (int j=0; j<nC; ++j) {
            dots[j] = Teuchos::ScalarTraits<Scalar>::zero();
          }
          return;
        }
        Node &node = A.getNode();
        DotOp2<Scalar,Node> op;
        for (int j=0; j<nC; ++j) {
          op.x = A.getValues(j);
          op.y = B.getValues(j);
          dots[j] = node.parallel_reduce(0,nR,op);
        }
      }

      inline static Scalar Dot(const MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &B) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same number of rows.");
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        Node &node = A.getNode();
        DotOp2<Scalar,Node> op;
        op.x = A.getValues(0);
        op.y = B.getValues(0);
        return node.parallel_reduce(0,nR,op);
      }

      inline static void GEMM(MultiVector<Scalar,Ordinal,Node> &C, Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &B, Scalar beta) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void GESUM(MultiVector<Scalar,Ordinal,Node> &B, Scalar alpha, const MultiVector<Scalar,Ordinal,Node> &A, Scalar beta) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(((nC != B.getNumCols()) && B.getNumCols() != 1)  ||
                           nR != B.getNumRows(), 
                           std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Divide(A,B): A and B must have the same dimensions.");
        Node &node = B.getNode();
        GESUMOp<Scalar,Node> wdp;
        wdp.alpha = alpha;
        wdp.beta  = beta;
        if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector
          wdp.y = B.getValues(0);
          wdp.x = ((const MultiVector<Scalar,Ordinal,Node> &)A).getValues(0);
          node.template parallel_for<GESUMOp<Scalar,Node> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (Ordinal j=0; j<nC; ++j) {
            wdp.y = B.getValues(j);
            wdp.x = ((const MultiVector<Scalar,Ordinal,Node> &)A).getValues(j);
            node.template parallel_for<GESUMOp<Scalar,Node> >(0,nR,wdp);
          }
        }
      }

      inline static void GESUM(MultiVector<Scalar,Ordinal,Node> &C, Scalar alpha, const MultiVector<Scalar,Ordinal,Node> &A, Scalar beta, const MultiVector<Scalar,Ordinal,Node> &B, Scalar gamma) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Norm1(const MultiVector<Scalar,Ordinal,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nC > norms.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (int j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        SumAbsOp<Scalar,Node> op;
        for (int j=0; j<nC; ++j) {
          op.x = A.getValues(j);
          norms[j] = node.parallel_reduce(0,nR,op);
        }
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType
      Norm1(const MultiVector<Scalar,Ordinal,Node> &A) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        SumAbsOp<Scalar,Node> op;
        op.x = A.getValues(0);
        return node.parallel_reduce(0,nR,op);
      }

      inline static void Sum(const MultiVector<Scalar,Ordinal,Node> &A, const Teuchos::ArrayView<Scalar> &sums) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nC > sums.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Sum(A,sums): sums must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
          for (int j=0; j<nC; ++j) {sums[j] = zero;}
          return;
        }
        Node &node = A.getNode();
        SumOp<Scalar,Node> op;
        for (int j=0; j<nC; ++j) {
          op.x = A.getValues(j);
          sums[j] = node.parallel_reduce(0,nR,op);
        }
      }

      inline static Scalar Sum(const MultiVector<Scalar,Ordinal,Node> &A) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        Node &node = A.getNode();
        SumOp<Scalar,Node> op;
        op.x = A.getValues(0);
        return node.parallel_reduce(0,nR,op);
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormInf(const MultiVector<Scalar,Ordinal,Node> &A) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        MaxAbsOp<Scalar,Node> op;
        op.x = A.getValues(0);
        return node.parallel_reduce(0,nR,op);
      }

      inline static void NormInf(const MultiVector<Scalar,Ordinal,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nC > norms.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::NormInf(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (int j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        MaxAbsOp<Scalar,Node> op;
        for (int j=0; j<nC; ++j) {
          op.x = A.getValues(j);
          norms[j] = node.parallel_reduce(0,nR,op);
        }
      }

      inline static void Norm2Squared(const MultiVector<Scalar,Ordinal,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        TEST_FOR_EXCEPTION(nC > norms.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm2Squared(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          for (int j=0; j<nC; ++j) {norms[j] = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();}
          return;
        }
        Node &node = A.getNode();
        DotOp1<Scalar,Node> op;
        for (int j=0; j<nC; ++j) {
          op.x = A.getValues(j);
          norms[j] = node.parallel_reduce(0,nR,op);
        }
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
      Norm2Squared(const MultiVector<Scalar,Ordinal,Node> &A) {
        const Ordinal nR = A.getNumRows();
        if (nR == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        Node &node = A.getNode();
        DotOp1<Scalar,Node> op;
        op.x = A.getValues(0);
        return node.parallel_reduce(0,nR,op);
      }

      inline static void WeightedNorm(const MultiVector<Scalar,Ordinal,Node> &A, const MultiVector<Scalar,Ordinal,Node> &weightVector, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Random(MultiVector<Scalar,Ordinal,Node> &A) {
        // TODO: consider adding rand() functionality to node
        // in the meantime, just generate random numbers via Teuchos and then copy to node
        typedef Teuchos::ScalarTraits<Scalar> SCT;
        const Ordinal stride = A.getStride();
        const Ordinal nR = A.getNumRows();
        const Ordinal nC = A.getNumCols();
        if (nR*nC == 0) return;
        Node &node = A.getNode();
        // we'll overwrite all data covered by the multivector, but not off-stride data
        // therefore, we are write-only only in the case that stride=nR
        bool writeOnly = (stride == nR);
        Scalar *mvdata = node.template viewBuffer<Scalar>(writeOnly,stride*(nC-1)+nR,A.getValues(0),0);
        for (int j=0; j<nC; ++j) {
          for (int i=0; i<nR; ++i) {
            mvdata[j*stride + i] = SCT::random();
          }
        }
        mvdata = Teuchos::null;
      }

      inline static void Abs(MultiVector<Scalar,Ordinal,Node> &B, const MultiVector<Scalar,Ordinal,Node> &A) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Scale(MultiVector<Scalar,Ordinal,Node> &B, Scalar alpha, const MultiVector<Scalar,Ordinal,Node> &A) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Scale(MultiVector<Scalar,Ordinal,Node> &B, Scalar alpha) {
        TEST_FOR_EXCEPT(true);
      }

      inline static void Scale(MultiVector<Scalar,Ordinal,Node> &B, const Teuchos::ArrayView<Scalar> &alphas) {
        TEST_FOR_EXCEPT(true);
      }

  };

}

#endif
