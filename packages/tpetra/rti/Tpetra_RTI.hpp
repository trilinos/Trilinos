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

#ifndef TPETRA_RTI_HPP
#define TPETRA_RTI_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_TestForException.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_Operator.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RTI_detail.hpp"

namespace Tpetra {

  namespace RTI {
    
    //! A static identity functor, providing a static method identity() that returns zero.
    template <class T>
    class ZeroOp {
      public:
      static inline T identity() {return Teuchos::ScalarTraits<T>::zero();}
    };

    //! A static identity functor, providing a static method identity() that returns one.
    template <class T>
    class OneOp {
      public:
      static inline T identity() {return Teuchos::ScalarTraits<T>::one();}
    };

    //! A type glob containing the types needed for calling Tpetra::RTI::binary_pre_transform_reduce() with individual functors.
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

    //! A type glob containing the types needed for calling Tpetra::RTI::reduce() with individual functors.
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
      TEST_FOR_EXCEPTION( vec_inout.getLocalLength() != vec_in2.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::binary_transform(vec_inout,vec_in2): vec_in2 and vec_inout must have the same local length.");
#endif
      Tpetra::RTI::detail::BinaryFunctorAdapter<OP,S1,S2> adapter_op(op);
      Tpetra::RTI::detail::binary_transform(vec_inout, vec_in2, adapter_op);
    }

    //! \brief Reduce values of \c vec_in1 and \c vec_in2 using the operators instantiated in \c glob.
    /** For each element pair <tt>vec_in1[i]</tt> and <tt>vec_in2[i]</tt>, generates reduction elements via <tt>glob.genop( vec_in1[i], vec_in2[i] )</tt> and reduces them via 
        the <tt>glob.redop</tt> binary functor. 

        Calls Tpetra::RTI::detail::reduce via the Tpetra::RTI::detail::RTIReductionAdapter.
      */
    template <class S1, class S2, class LO, class GO, class Node, class Glob>
    typename Glob::RedOP::result_type 
    reduce( const Vector<S1,LO,GO,Node> &vec_in1, const Vector<S2,LO,GO,Node> &vec_in2, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( vec_in1.getLocalLength() != vec_in2.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::binary_transform(vec_in1,vec_in2): vec_in1 and vec_in2 must have the same local length.");
#endif
      Tpetra::RTI::detail::RTIReductionAdapter<Glob,S1,S2> adapter_op(glob);
      return Tpetra::RTI::detail::reduce(vec_in1, vec_in2, adapter_op);
    }

    //! \brief Reduce values of \c vec_in1, \c vec_in2  and \c vec_in3 using the operators instantiated in \c glob.
    /** For each element triplet <tt>vec_in1[i]</tt>, <tt>vec_in2[i]</tt> and <tt>vec_in3[i]</tt>, generates reduction elements via 
        <tt>glob.genop( vec_in1[i], vec_in2[i], vec_in3[i] )</tt> and reduces them via 
        the <tt>glob.redop</tt> binary functor. 

        Calls Tpetra::RTI::detail::reduce via the Tpetra::RTI::detail::RTIReductionAdapter.
      */
    template <class S1, class S2, class S3, class LO, class GO, class Node, class Glob>
    typename Glob::RedOP::result_type 
    reduce(const Vector<S1,LO,GO,Node> &vec_in1, const Vector<S2,LO,GO,Node> &vec_in2, const Vector<S3,LO,GO,Node> &vec_in3, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( vec_in1.getLocalLength() != vec_in2.getLocalLength() || vec_in2.getLocalLength() != vec_in3.getLocalLength(), 
          std::runtime_error, "Tpetra::RTI::binary_transform(vec_in1,vec_in2): vec_in1 and vec_in2 must have the same local length.");
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
    template <class S, class LO, class GO, class Node,class Glob>
    typename Glob::RedOP::result_type 
    binary_pre_transform_reduce(Vector<S,LO,GO,Node> &vec_inout, const Vector<S,LO,GO,Node> &vec_in2, Glob glob)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( vec_inout.getLocalLength() != vec_in2.getLocalLength(), std::runtime_error,
          "Tpetra::RTI::binary_transform_reduce(vec_in1,vec_in2): vec_in1 and vec_in2 must have the same local length.");
#endif
      Tpetra::RTI::detail::RTIPreTransformReductionAdapter<Glob,S,S> adapter_op(glob);
      return Tpetra::RTI::detail::transform_reduce(vec_inout, vec_in2, adapter_op);
    }

  } // end of namespace RTI

} // end of namespace Tpetra

#endif // TPETRA_RTI_HPP
