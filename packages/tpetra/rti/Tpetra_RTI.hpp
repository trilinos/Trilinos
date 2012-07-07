// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_RTI_HPP
#define TPETRA_RTI_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_Assert.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_Operator.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RTI_detail.hpp"

namespace Tpetra {

  namespace RTI {
    
    /// \class ZeroOp
    /// \brief A static identity functor, providing a static method identity() that returns zero.
    ///
    /// This class is useful for providing zero as the initial value
    /// of a reduction.  It may be used as the IOP template parameter
    /// of a ReductionGlob.
    template <class T>
    class ZeroOp {
      public:
      static inline T identity() {return Teuchos::ScalarTraits<T>::zero();}
    };

    /// \class OneOp
    /// \brief A static identity functor, providing a static method identity() that returns one.
    /// 
    /// This class is useful for providing one as the initial value of
    /// a reduction.  It may be used as the IOP template parameter of
    /// a ReductionGlob.
    template <class T>
    class OneOp {
      public:
      static inline T identity() {return Teuchos::ScalarTraits<T>::one();}
    };

    /// \class ReductionGlob
    /// \brief A type glob containing the types needed for calling Tpetra::RTI::reduce() with individual functors.
    ///
    /// \tparam GOP Type of the operator genop, that generates
    ///   successive new inputs of the reduction.
    /// 
    /// \tparam ROP Type of the operator that performs the pairwise
    ///   reduction operations.
    ///
    /// \tparam IOP Type of the operator that provides (via a
    ///   zero-argument static function) the initial value of the
    ///   reduction.
    ///
    /// For reducing a pair of vectors v1, v2, successive reduction
    /// elements are generated in a way equivalent to <tt>genop(v1[i],
    /// v2[i])</tt> for all indices i of the vector.
    ///
    /// For reducing a triple of vectors v1, v2, v3, successive
    /// reduction elements are generated in a way equivalent to
    /// <tt>genop(v1[i], v2[i], v3[i])</tt> for all indices i of the
    /// vector.
    ///
    /// Regardless, each genop invocation generates a single value,
    /// and the sequence of these values is reduced using the binary
    /// operator redop.  The initial value of this sequence comes from
    /// the static <tt>identity()</tt> method of IOP.
    template <class GOP, class ROP, class IOP> 
    class ReductionGlob {
      public:
        typedef GOP GenOP;
        typedef ROP RedOP;
        typedef IOP  IdOP;
        GenOP genop;
        RedOP redop;
        ReductionGlob(GenOP gop, RedOP rop) : genop(gop), redop(rop) {}
    };

    //! A type glob containing the types needed for calling Tpetra::RTI::binary_pre_transform_reduce() with individual functors.
    template <class TxOP, class GOP, class ROP, class IOP> 
    class TransformReductionGlob {
      public:
        typedef TxOP   TOP;
        typedef GOP  GenOP;
        typedef ROP  RedOP;
        typedef IOP   IdOP;
        TOP     top;
        GenOP genop;
        RedOP redop;
        TransformReductionGlob(TOP txop, GenOP gop, RedOP rop) : top(txop), genop(gop), redop(rop) {}
    };

    //! Non-member constructor to instantiate a type glob of a static identity functor and generation and reduction functor objects.
    template <class IOP, class GOP, class ROP> 
    inline ReductionGlob<GOP,ROP,IOP> reductionGlob(GOP gop, ROP rop) 
    {
      return ReductionGlob<GOP,ROP,IOP>(gop,rop);
    }

    //! Non-member constructor to instantiate a type glob of a static identity functor and transform, generation and reduction functor objects.
    template <class IOP, class TOP, class GOP, class ROP> 
    inline TransformReductionGlob<TOP,GOP,ROP,IOP> reductionGlob(TOP top, GOP gop, ROP rop) 
    {
      return TransformReductionGlob<TOP,GOP,ROP,IOP>(top,gop,rop);
    }

    //! \brief Transform values of \c vec_inout using via operator \c op.
    /** For each element <tt>vec_inout[i]</tt>, assign <tt>vec_inout[i] = op( vec_inout[i] )</tt>
        
        Calls Tpetra::RTI::detail::unary_transform via the Tpetra::RTI::detail::UnaryFunctorAdapter.
      */
    template <class S, class LO, class GO, class Node, class OP>
    void unary_transform(Vector<S,LO,GO,Node> &vec_inout, OP op) 
    {
      Tpetra::RTI::detail::UnaryFunctorAdapter<OP,S> Adapter_op(op);
      Tpetra::RTI::detail::unary_transform(vec_inout, Adapter_op);
    }

    //! \brief Transform values of \c vec_inout using \c vec_inout, \c vec_in2 and operator \c op.
    /** For each element <tt>vec_inout[i]</tt>, assign <tt>vec_inout[i] = op( vec_inout[i], vec_in2[i] )</tt>
        
        Calls Tpetra::RTI::detail::binary_transform via the Tpetra::RTI::detail::BinaryFunctorAdapter.
      */
    template <class S1, class S2, class LO, class GO, class Node, class OP>
    void binary_transform(Vector<S1,LO,GO,Node> &vec_inout, const Vector<S2,LO,GO,Node> &vec_in2, OP op) 
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( vec_inout.getLocalLength() != vec_in2.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::binary_transform(vec_inout,vec_in2): vec_in2 and vec_inout must have the same local length.");
#endif
      Tpetra::RTI::detail::BinaryFunctorAdapter<OP,S1,S2> adapter_op(op);
      Tpetra::RTI::detail::binary_transform(vec_inout, vec_in2, adapter_op);
    }

    //! \brief Transform values of \c vec_inout using \c vec_inout, \c vec_in2, \c vec_in3 and operator \c op.
    /** For each element <tt>vec_inout[i]</tt>, assign <tt>vec_inout[i] = op( vec_inout[i], vec_in2[i], vec_in3[i] )</tt>
        
        Calls Tpetra::RTI::detail::tertiary_transform via the Tpetra::RTI::detail::TertiaryFunctorAdapter.
      */
    template <class S1, class S2, class S3, class LO, class GO, class Node, class OP>
    void tertiary_transform(Vector<S1,LO,GO,Node> &vec_inout, const Vector<S2,LO,GO,Node> &vec_in2, const Vector<S3,LO,GO,Node> &vec_in3, OP op) 
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( vec_inout.getLocalLength() != vec_in2.getLocalLength() || vec_in2.getLocalLength() != vec_in3.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::tertiary_transform(vec_inout,vec_in2,vec_in3): vec_inout, vec_in2 and vec_in3 must have the same local length.");
#endif
      Tpetra::RTI::detail::TertiaryFunctorAdapter<OP,S1,S2,S3> adapter_op(op);
      Tpetra::RTI::detail::tertiary_transform(vec_inout, vec_in2, vec_in3, adapter_op);
    }

    //! \brief Reduce values of \c vec_in using the operators instantiated in \c glob.
    /** For each element <tt>vec_in[i]</tt>, generates reduction elements via <tt>glob.genop( vec_in[i] )</tt> and reduces them via 
        the <tt>glob.redop</tt> binary functor. 

        Calls Tpetra::RTI::detail::reduce via the Tpetra::RTI::detail::RTIReductionAdapter1.
      */
    template <class S, class LO, class GO, class Node, class Glob>
    typename Glob::RedOP::result_type 
    reduce( const Vector<S,LO,GO,Node> &vec_in, Glob glob)
    {
      Tpetra::RTI::detail::RTIReductionAdapter1<Glob,S> adapter_op(glob);
      return Tpetra::RTI::detail::reduce(vec_in, adapter_op);
    }

    //! \brief Reduce values of \c vec_in1 and \c vec_in2 using the operators instantiated in \c glob.
    /** For each element pair <tt>vec_in1[i]</tt> and <tt>vec_in2[i]</tt>, generates reduction elements via <tt>glob.genop( vec_in1[i], vec_in2[i] )</tt> and reduces them via 
        the <tt>glob.redop</tt> binary functor. 

        Calls Tpetra::RTI::detail::reduce via the Tpetra::RTI::detail::RTIReductionAdapter2.
      */
    template <class S1, class S2, class LO, class GO, class Node, class Glob>
    typename Glob::RedOP::result_type 
    reduce( const Vector<S1,LO,GO,Node> &vec_in1, const Vector<S2,LO,GO,Node> &vec_in2, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( vec_in1.getLocalLength() != vec_in2.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::reduce(vec_in1,vec_in2): vec_in1 and vec_in2 must have the same local length.");
#endif
      Tpetra::RTI::detail::RTIReductionAdapter2<Glob,S1,S2> adapter_op(glob);
      return Tpetra::RTI::detail::reduce(vec_in1, vec_in2, adapter_op);
    }

    //! \brief Reduce values of \c vec_in1, \c vec_in2  and \c vec_in3 using the operators instantiated in \c glob.
    /** For each element triplet <tt>vec_in1[i]</tt>, <tt>vec_in2[i]</tt> and <tt>vec_in3[i]</tt>, generates reduction elements via 
        <tt>glob.genop( vec_in1[i], vec_in2[i], vec_in3[i] )</tt> and reduces them via 
        the <tt>glob.redop</tt> binary functor. 

        Calls Tpetra::RTI::detail::reduce via the Tpetra::RTI::detail::RTIReductionAdapter3.
      */
    template <class S1, class S2, class S3, class LO, class GO, class Node, class Glob>
    typename Glob::RedOP::result_type 
    reduce(const Vector<S1,LO,GO,Node> &vec_in1, const Vector<S2,LO,GO,Node> &vec_in2, const Vector<S3,LO,GO,Node> &vec_in3, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( vec_in1.getLocalLength() != vec_in2.getLocalLength() || vec_in2.getLocalLength() != vec_in3.getLocalLength(), 
          std::runtime_error, "Tpetra::RTI::reduce(vec_in1,vec_in2): vec_in1 and vec_in2 must have the same local length.");
#endif
      Tpetra::RTI::detail::RTIReductionAdapter3<Glob,S1,S2,S3> adapter_op(glob);
      return Tpetra::RTI::detail::reduce(vec_in1, vec_in2, vec_in3, adapter_op);
    }

    //! \brief Transforms values of \c vec_inout while simultaneously performing a parallel reduction.
    /** For each element pair <tt>vec_inout[i]</tt> and <tt>vec_in2[i]</tt>, 
        assigns <tt>vec_inout[i] = glob.top( vec_inout[i], vec_in2[i] )</tt>. Simultaneously, generates reduction elements via <tt>glob.genop( vec_inout[i], vec_in2[i] )</tt> (using the transformed values) and reduces them via 
        the <tt>glob.redop</tt> binary functor. 
        
        Calls Tpetra::RTI::detail::transform_reduce via the Tpetra::RTI::detail::RTIPreTransformReductionAdapter.
      */
    template <class S1, class S2, class LO, class GO, class Node,class Glob>
    typename Glob::RedOP::result_type 
    binary_pre_transform_reduce(Vector<S1,LO,GO,Node> &vec_inout, const Vector<S2,LO,GO,Node> &vec_in2, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( vec_inout.getLocalLength() != vec_in2.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::binary_pre_transform_reduce(vec_in1,vec_in2): vec_in1 and vec_in2 must have the same local length.");
#endif
      Tpetra::RTI::detail::RTIPreTransformReductionAdapter<Glob,S1,S2> adapter_op(glob);
      return Tpetra::RTI::detail::transform_reduce(vec_inout, vec_in2, adapter_op);
    }

    //! \brief Transforms values of \c vec_inout while simultaneously performing a parallel reduction.
    /** For each element triplet <tt>vec_inout[i]</tt> and <tt>vec_in2[i]</tt> and <tt>vec_in3[i]</tt>, 
        assigns <tt>vec_inout[i] = glob.top( vec_inout[i], vec_in2[i], vec_in3[i] )</tt>. Simultaneously, generates reduction 
        elements via <tt>glob.genop( vec_inout[i], vec_in2[i], vec_in3[i] )</tt> (using the transformed values) and reduces them via 
        the <tt>glob.redop</tt> tertiary functor. 
        
        Calls Tpetra::RTI::detail::transform_reduce via the Tpetra::RTI::detail::RTIPreTransformReductionAdapter3.
      */
    template <class S1, class S2, class S3, class LO, class GO, class Node,class Glob>
    typename Glob::RedOP::result_type 
    tertiary_pre_transform_reduce(Vector<S1,LO,GO,Node> &vec_inout, const Vector<S2,LO,GO,Node> &vec_in2, const Vector<S3,LO,GO,Node> &vec_in3, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( vec_inout.getLocalLength() != vec_in2.getLocalLength() && vec_in2.getLocalLength() != vec_in3.getLocalLength(), 
          std::runtime_error, "Tpetra::RTI::tertiary_pre_transform_reduce(vec_in1,vec_in2,vec_in3): vec_in1, vec_in2 and vec_in3 must have the same local length.");
#endif
      Tpetra::RTI::detail::RTIPreTransformReductionAdapter3<Glob,S1,S2,S3> adapter_op(glob);
      return Tpetra::RTI::detail::transform_reduce(vec_inout, vec_in2, vec_in3, adapter_op);
    }

  } // end of namespace RTI

} // end of namespace Tpetra

#define TPETRA_UNARY_TRANSFORM(out,expr) \
  Tpetra::RTI::unary_transform( *out, [=](decltype((out)->meanValue()) out) \
                                         {return expr;})

#define TPETRA_BINARY_TRANSFORM(out,in,expr) \
  Tpetra::RTI::binary_transform( *out, *in, [=](decltype((out)->meanValue()) out, \
                                                    decltype((in)->meanValue()) in) \
                                                    {return expr;})

#define TPETRA_TERTIARY_TRANSFORM(out,in2,in3,expr) \
  Tpetra::RTI::tertiary_transform( *out, *in2, *in3, \
                                   [=](decltype((out)->meanValue()) out, \
                                       decltype((in2)->meanValue()) in2, \
                                       decltype((in3)->meanValue()) in3) \
                                       {return expr;})


#define TPETRA_REDUCE1(in, gexp, id, robj ) \
  Tpetra::RTI::reduce( *in,                                                      \
    Tpetra::RTI::reductionGlob<id>( [=]( decltype((in)->meanValue()) in )        \
                                       { return gexp; },                         \
                                       robj ) )

#define TPETRA_REDUCE2(in1,in2, gexp, id, robj ) \
  Tpetra::RTI::reduce( *in1, *in2,                                                \
    Tpetra::RTI::reductionGlob<id>( [=]( decltype((in1)->meanValue()) in1,       \
                                         decltype((in2)->meanValue()) in2 )        \
                                       { return gexp; },                         \
                                       robj ) )

#define TPETRA_REDUCE3(in1,in2,in3, gexp, id, robj ) \
  Tpetra::RTI::reduce( *in1, *in2, *in3,                  \
    Tpetra::RTI::reductionGlob<id>( [=]( decltype((in1)->meanValue()) in1,       \
                                         decltype((in2)->meanValue()) in2,       \
                                         decltype((in3)->meanValue()) in3 )      \
                                       { return gexp; },                         \
                                       robj ) )

#define TPETRA_BINARY_PRETRANSFORM_REDUCE(out,in,  texp, gexp, id, robj ) \
  Tpetra::RTI::binary_pre_transform_reduce( *out, *in,                     \
    Tpetra::RTI::reductionGlob<id>( [=]( decltype((out)->meanValue()) out,  \
                                         decltype((in)->meanValue()) in )   \
                                       { return texp; },                    \
                                    [=]( decltype((out)->meanValue()) out,  \
                                         decltype((in)->meanValue()) in )   \
                                       { return gexp; },                    \
                                       robj ) )

#define TPETRA_TERTIARY_PRETRANSFORM_REDUCE(out,in2,in3, texp, gexp, id, robj ) \
  Tpetra::RTI::tertiary_pre_transform_reduce( *out, *in2, *in3,                  \
    Tpetra::RTI::reductionGlob<id>( [=]( decltype((out)->meanValue()) out,       \
                                         decltype((in2)->meanValue()) in2,       \
                                         decltype((in3)->meanValue()) in3 )      \
                                       { return texp; },                         \
                                    [=]( decltype((out)->meanValue()) out,       \
                                         decltype((in2)->meanValue()) in2,       \
                                         decltype((in3)->meanValue()) in3 )      \
                                       { return gexp; },                         \
                                       robj ) )

#endif // TPETRA_RTI_HPP
