//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// Lesser General Public License for more detail.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef TPETRA_RTI_detail_HPP
#define TPETRA_RTI_detail_HPP

#include <Teuchos_TestForException.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Kokkos_NodeHelpers.hpp>

#include "Tpetra_Vector.hpp"

namespace Tpetra {

  namespace RTI {

    //! Internal detail for Tpetra::RTI. Methods and class here are not guaranteed to be backwards compatible.
    namespace detail {

      //! Utility base class for kernels used to define Tpetra::Operator objects.
      template <class S>
      class StdOpKernel
      {
        protected:
          S _alpha, _beta;
          S       * _vec_inout;
          const S * _vec_in2;
        public:
          inline StdOpKernel() : _alpha(ScalarTraits<S>::one()), _beta(ScalarTraits<S>::zero()) {}
          inline void setData(S * vec_inout, const S * vec_in2)   { _vec_inout = vec_inout; _vec_in2 = vec_in2; }
          inline void setAlphaBeta(const S &alpha, const S &beta) { _alpha = alpha; _beta = beta; }
      };

      //! adapter class between kernels for Tpetra::RTI::unary_transform and Tpetra::RTI::detail::unary_transform
      template <class OP, class S>
      class UnaryFunctorAdapter {
        protected:
          OP   _op;
          S  * _vec;
        public:
          UnaryFunctorAdapter(OP op) : _op(op) {}
          inline void setData(S *vec)                    { _vec = vec;             }
          inline void execute(int i)                     { _vec[i] = _op(_vec[i]); }
      };

      //! adapter class between kernels for Tpetra::RTI::binary_transform and Tpetra::RTI::detail::binary_transform
      template <class OP, class S1, class S2>
      class BinaryFunctorAdapter {
        protected:
          OP         _op;
          S1       * _vec_inout;
          const S2 * _vec_in2;
        public:
          BinaryFunctorAdapter(OP op) : _op(op)    {}
          inline void setData(S1 *vec_inout, const S2 *vec_in2) { _vec_inout = vec_inout; _vec_in2 = vec_in2;   }
          inline void execute(int i)                            { _vec_inout[i] = _op(_vec_inout[i], _vec_in2[i]); }
      };

      //! adapter class between binary functors and BinaryOp
      template <class OP, class S>
      class BinaryFunctorAdapterWithAlphaBeta : public StdOpKernel<S> {
        protected:
          OP        _op;
          S       * _vec_inout;
          const S * _vec_in2;
        public:
          BinaryFunctorAdapterWithAlphaBeta(OP op) : _op(op)  {}
          inline void setData(S *vec_inout, const S *vec_in2) { _vec_inout = vec_inout; _vec_in2 = vec_in2;   }
          inline void execute(int i) 
          { 
            S res = _op(_vec_inout[i], _vec_in2[i]); 
            _vec_inout[i] = this->_alpha * res + this->_beta * _vec_inout[i];
          }
      };

      //! adapter class between kernels for Tpetra::RTI::binary_transform and Tpetra::RTI::detail::binary_transform
      template <class Glob, class S1, class S2>
      class RTIReductionAdapter {
        public:
          typedef typename Glob::GenOP                GenOP;
          typedef typename Glob::RedOP                RedOP;
          typedef typename Glob::IdOP                  IdOP;
          typedef typename RedOP::result_type ReductionType;
        protected:
          GenOP      _genop;
          RedOP      _redop;
          const S1 * _vec_in1;
          const S2 * _vec_in2;
        public:
          RTIReductionAdapter(Glob glob)                                : _genop(glob.genop), _redop(glob.redop)  {}
          inline void setData(const S1 *vec_in1, const S2 *vec_in2)     { _vec_in1 = vec_in1; _vec_in2 = vec_in2;  }
          static inline ReductionType identity()                        { return IdOP::identity();                 }
          inline ReductionType generate(int i)                          { return _genop(_vec_in1[i], _vec_in2[i]); }
          inline ReductionType reduce(ReductionType a, ReductionType b) { return _redop(a, b);                     }
      };

      //! adapter class between kernels for Tpetra::RTI::binary_transform and Tpetra::RTI::detail::binary_transform for three vectors
      template <class Glob, class S1, class S2, class S3>
      class RTIReductionAdapter3 {
        public:
          typedef typename Glob::GenOP                GenOP;
          typedef typename Glob::RedOP                RedOP;
          typedef typename Glob::IdOP                  IdOP;
          typedef typename RedOP::result_type ReductionType;
        protected:
          GenOP      _genop;
          RedOP      _redop;
          const S1 * _vec_in1;
          const S2 * _vec_in2;
          const S3 * _vec_in3;
        public:
          RTIReductionAdapter3(Glob glob)                                              : _genop(glob.genop), _redop(glob.redop)                     {}
          inline void setData(const S1 *vec_in1, const S2 *vec_in2, const S3 *vec_in3) { _vec_in1 = vec_in1; _vec_in2 = vec_in2; _vec_in3 = vec_in3; }
          static inline ReductionType identity()                                       { return IdOP::identity();                                    }
          inline ReductionType generate(int i)                                         { return _genop(_vec_in1[i], _vec_in2[i], _vec_in3[i]);       }
          inline ReductionType reduce(ReductionType a, ReductionType b)                { return _redop(a, b);                                        }
      };

      //! adapter class between kernels for Tpetra::RTI::binary_pre_transform_reduce and Tpetra::RTI::detail::binary_transform
      template <class Glob, class S1, class S2>
      class RTIPreTransformReductionAdapter {
        public:
          typedef typename Glob::TOP                    TOP;
          typedef typename Glob::GenOP                GenOP;
          typedef typename Glob::RedOP                RedOP;
          typedef typename Glob::IdOP                  IdOP;
          typedef typename RedOP::result_type ReductionType;
        protected:
          TOP        _top;
          GenOP      _genop;
          RedOP      _redop;
          S1       * _vec_inout;
          const S2 * _vec_in2;
        public:
          RTIPreTransformReductionAdapter(Glob glob)                    : _top(glob.top), _genop(glob.genop), _redop(glob.redop) {}
          inline void setData(S1 *vec_inout, const S2 *vec_in2)         { _vec_inout = vec_inout; _vec_in2 = vec_in2; }
          static inline ReductionType identity()                        { return IdOP::identity();                    }
          inline ReductionType reduce(ReductionType a, ReductionType b) { return _redop(a, b);                        }
          inline ReductionType generate(int i)
          {
            _vec_inout[i] = _top( _vec_inout[i], _vec_in2[i] );
            return _genop( _vec_inout[i], _vec_in2[i] ); 
          }
      };

      //! decorator for Kokkos reduction kernels to satisfy requirements for Teuchos::ValueTypeReductionOp
      template <class OP>
      class TeuchosValueTypeReductionOpAdapter : public Teuchos::ValueTypeReductionOp<int,typename OP::ReductionType> 
      {
        protected:
          mutable OP _op;
        public:
          typedef typename OP::ReductionType Packet;
          TeuchosValueTypeReductionOpAdapter(OP op) : _op(op) {}
          void reduce(const int count, const Packet inBuffer[], Packet inoutBuffer []) const 
          {
            for (int i=0; i != count; ++i) {
              inoutBuffer[i] = _op.reduce( inoutBuffer[i], inBuffer[i] );
            }
          }
      };

      //! pass \c vec data pointer to \c op, then execute via node parallel_for
      template <class S, class LO, class GO, class Node, class OP>
      void unary_transform(Vector<S,LO,GO,Node> &vec, OP op) 
      {
        Kokkos::MultiVector<S,Node> &mv = vec.getLocalMVNonConst();
        const RCP<Node> node = mv.getNode();
        // ready data
        Kokkos::ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        S * out_ptr = rbh.addNonConstBuffer(mv.getValuesNonConst());
        rbh.end();
        op.setData(out_ptr);
        const size_t N = mv.getNumRows();
        node->template parallel_for(0, N, op);
      }

      //! pass \c vec_inout and \c vec_in2 data pointers to \c op, then execute via node parallel_for
      template <class S1, class S2, class LO, class GO, class Node, class OP>
      void binary_transform(Vector<S1,LO,GO,Node> &vec_inout, const Vector<S2,LO,GO,Node> &vec_in2, OP op) 
      {
        Kokkos::MultiVector<S1,Node>       &mv_inout = vec_inout.getLocalMVNonConst();
        const Kokkos::MultiVector<S2,Node> &mv_in2   = vec_in2.getLocalMV();
        const RCP<Node> node = mv_inout.getNode();
        // ready data
        Kokkos::ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        S1       * out_ptr = rbh.addNonConstBuffer(mv_inout.getValuesNonConst());
        const S2 * in_ptr  = rbh.addConstBuffer(mv_in2.getValues());
        rbh.end();
        op.setData(out_ptr, in_ptr);
        const size_t N = mv_inout.getNumRows();
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION( mv_in2.getNode() != mv_inout.getNode(), std::runtime_error, 
            "Tpetra::RTI::detail::binary_transform(): multivectors must share the same node.");
#endif
        node->template parallel_for(0, N, op);
      }

      //! pass \c vec_in1 and \c vec_in2 data pointers to \ op, then execute via node parallel_reduce.
      template <class S1, class S2, class LO, class GO, class Node, class OP>
      typename OP::ReductionType 
      reduce(const Vector<S1,LO,GO,Node> &vec_in1, const Vector<S2,LO,GO,Node> &vec_in2, OP op) 
      {
        const Kokkos::MultiVector<S1,Node> &mv_in1 = vec_in1.getLocalMV(),
                                           &mv_in2 = vec_in2.getLocalMV();
        const RCP<Node> node = mv_in1.getNode();
        const RCP<const Teuchos::Comm<int> > comm = vec_in1.getMap()->getComm();
        // ready data
        Kokkos::ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        const S1 * in_ptr1 = rbh.addConstBuffer(mv_in1.getValues());
        const S2 * in_ptr2 = rbh.addConstBuffer(mv_in2.getValues());
        rbh.end();
        op.setData( in_ptr1, in_ptr2 );
        const size_t N = mv_in1.getNumRows();
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION( mv_in1.getNode() != mv_in2.getNode(), std::runtime_error, 
            "Tpetra::RTI::detail::reduce(): multivectors must share the same node.");
#endif
        // compute local reduction
        typename OP::ReductionType gbl_res, lcl_res;
        lcl_res = node->template parallel_reduce(0, N, op);
        // compute global reduction
        TeuchosValueTypeReductionOpAdapter<OP> vtrop(op);
        Teuchos::reduceAll(*comm, vtrop, 1, &lcl_res, &gbl_res);
        return gbl_res;
      }

      //! pass \c vec_in1, \c vec_in2 and \c vec_in3 data pointers to \ op, then execute via node parallel_reduce.
      template <class S1, class S2, class S3, class LO, class GO, class Node, class OP>
      typename OP::ReductionType 
      reduce(const Vector<S1,LO,GO,Node> &vec_in1, const Vector<S2,LO,GO,Node> &vec_in2, const Vector<S3,LO,GO,Node> &vec_in3, OP op) 
      {
        const Kokkos::MultiVector<S1,Node> &mv_in1 = vec_in1.getLocalMV();
        const Kokkos::MultiVector<S2,Node> &mv_in2 = vec_in2.getLocalMV();
        const Kokkos::MultiVector<S3,Node> &mv_in3 = vec_in3.getLocalMV();
        const RCP<Node> node = mv_in1.getNode();
        const RCP<const Teuchos::Comm<int> > comm = vec_in1.getMap()->getComm();
        // ready data
        Kokkos::ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        const S1 * in_ptr1 = rbh.addConstBuffer(mv_in1.getValues());
        const S2 * in_ptr2 = rbh.addConstBuffer(mv_in2.getValues());
        const S3 * in_ptr3 = rbh.addConstBuffer(mv_in3.getValues());
        rbh.end();
        op.setData( in_ptr1, in_ptr2, in_ptr3 );
        const size_t N = mv_in1.getNumRows();
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION( mv_in1.getNode() != mv_in2.getNode() || mv_in2.getNode() != mv_in3.getNode(), std::runtime_error, 
            "Tpetra::RTI::detail::reduce(): multivectors must share the same node.");
#endif
        // compute local reduction
        typename OP::ReductionType gbl_res, lcl_res;
        lcl_res = node->template parallel_reduce(0, N, op);
        // compute global reduction
        TeuchosValueTypeReductionOpAdapter<OP> vtrop(op);
        Teuchos::reduceAll(*comm, vtrop, 1, &lcl_res, &gbl_res);
        return gbl_res;
      }

      //! pass \c vec_inout and \c vec_in2 data pointers to \ op, then execute via node parallel_reduce.
      template <class S1, class S2, class LO, class GO, class Node, class OP>
      typename OP::ReductionType 
      transform_reduce(Vector<S1,LO,GO,Node> &vec_inout, const Vector<S2,LO,GO,Node> &vec_in2, OP op) 
      {
        Kokkos::MultiVector<S1,Node>       &mv_inout = vec_inout.getLocalMVNonConst();
        const Kokkos::MultiVector<S2,Node> &mv_in2   = vec_in2.getLocalMV();
        const RCP<Node> node = mv_inout.getNode();
        const RCP<const Teuchos::Comm<int> > comm = vec_inout.getMap()->getComm();
        // ready data
        Kokkos::ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        S1 * in_ptr1 = rbh.addNonConstBuffer(mv_inout.getValuesNonConst());
        const S2 * in_ptr2 = rbh.addConstBuffer(mv_in2.getValues());
        rbh.end();
        op.setData( in_ptr1, in_ptr2 );
        const size_t N = mv_inout.getNumRows();
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION( mv_inout.getNode() != mv_in2.getNode(), std::runtime_error, 
            "Tpetra::RTI::detail::binary_transform(): multivectors must share the same node.");
#endif
        // compute local reduction
        typename OP::ReductionType gbl_res, lcl_res;
        lcl_res = node->template parallel_reduce(0, N, op);
        // compute global reduction
        TeuchosValueTypeReductionOpAdapter<OP> vtrop(op);
        Teuchos::reduceAll(*comm, vtrop, 1, &lcl_res, &gbl_res);
        return gbl_res;
      }

    } // end of namespace Tpetra::RTI::detail

  } // end of namespace Tpetra::RTI

} // end of namespace Tpetra

#endif // TPETRA_RTI_detail_HPP
