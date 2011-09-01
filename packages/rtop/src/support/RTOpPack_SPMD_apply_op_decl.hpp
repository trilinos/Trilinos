// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_SPMD_APPLY_OP_DECL_HPP
#define RTOPPACK_SPMD_APPLY_OP_DECL_HPP

#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_Serializer.hpp"
#include "Teuchos_ReductionOp.hpp"


namespace Teuchos { template<typename Ordinal> class Comm; }


// Enable this by hand to enable showing the dump of the RTOp
#define RTOPPACK_ENABLE_SHOW_DUMP


#ifdef RTOP_DEBUG
#  define RTOPPACK_ENABLE_SHOW_DUMP
#endif


namespace RTOpPack {


#ifdef RTOPPACK_ENABLE_SHOW_DUMP
extern bool show_spmd_apply_op_dump;
#endif // RTOPPACK_ENABLE_SHOW_DUMP


/** \brief Return the size in bytes of an external representation of a
 * <tt>ReductTarget</tt> object.
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class PrimitiveScalar>
int serializedSize(
  int   num_values
  ,int  num_indexes
  ,int  num_chars
  );


/** \brief Serialize a <tt>ReductTarget</tt> object.
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void serialize(
  const RTOpT<Scalar> &op,
  Ordinal num_values,
  Ordinal num_indexes,
  Ordinal num_chars,
  const ReductTarget &reduct_obj,
  char reduct_obj_ext[]
  );


/** \brief Deserialize a <tt>ReductTarget</tt> object.
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void deserialize(
  const RTOpT<Scalar> &op,
  int num_values,
  int num_indexes,
  int num_chars,
  const char reduct_obj_ext[],
  ReductTarget *reduct_obj
  );


/** \brief Serializer subclass for <tt>ReductTarget</tt> objects.
 *
 * The copy constructor is allowed and has shallow copy semantics.
 */
template<class Scalar>
class ReductTargetSerializer : public Teuchos::Serializer<index_type,ReductTarget> {
public:
  /** \brief . */
  ReductTargetSerializer(
    const Teuchos::RCP<const RTOpT<Scalar> > &op
    );
  /** \name Public functions overridden from Teuchos::Serializer */
  //@{
  /** \brief . */
  index_type getBufferSize(const index_type count) const;
  /** \brief . */
  void serialize(
    const index_type            count
    ,const ReductTarget* const  reduct_objs[]
    ,const index_type           bytes
    ,char                       charBuffer[]
    ) const;
  /** \brief . */
  Teuchos::RCP<ReductTarget> createObj() const;
  /** \brief . */
  void deserialize(
    const index_type       bytes
    ,const char            charBuffer[]
    ,const index_type      count
    ,ReductTarget* const   reduct_objs[]
    ) const;
  //@}
private:
  Teuchos::RCP<const RTOpT<Scalar> >   op_;
  int                                          num_values_;
  int                                          num_indexes_;
  int                                          num_chars_;
  int                                          reduct_obj_ext_size_;
  // Not defined and not to be called!
  ReductTargetSerializer();
  ReductTargetSerializer& operator=(const ReductTargetSerializer&);
};


/** \brief ReductionOp subclass for <tt>ReductTarget</tt> objects.
 *
 * The copy constructor is allowed and has shallow copy semantics.
 */
template<class Scalar>
class ReductTargetReductionOp
  : public Teuchos::ReferenceTypeReductionOp<Teuchos_Index,ReductTarget>
{
public:
  /** \brief . */
  typedef Teuchos_Index Ordinal;
  /** \brief . */
  ReductTargetReductionOp(
    const Teuchos::RCP<const RTOpT<Scalar> >  &op
    );
  /** \name Overridden from Teuchos::ReferenceTypeReductionOp */
  //@{
  /** \brief . */
  void reduce(
    const Ordinal              count
    ,const ReductTarget*const  inBuffer[]
    ,ReductTarget*const        inoutBuffer[]
    ) const;
  //@}
private:
  Teuchos::RCP<const RTOpT<Scalar> >  op_;
  // Not defined and not to be called!
  ReductTargetReductionOp();
  ReductTargetReductionOp<Scalar>(const ReductTargetReductionOp<Scalar>&);
  ReductTargetReductionOp<Scalar>& operator=(const ReductTargetReductionOp<Scalar>&);
};


/** \brief Reduce a set of reduction objects.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void SPMD_all_reduce(
  const Teuchos::Comm<index_type>     *comm
  ,const RTOpT<Scalar>                &op
  ,const int                          num_cols
  ,const ReductTarget*const           i_reduct_objs[]
  ,ReductTarget*const                 reduct_objs[]
  );


/** \brief Apply an RTOp in SMPD mode to a set of vectors with contiguous
 * storage per process.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void SPMD_apply_op(
  const Teuchos::Comm<index_type>     *comm
  ,const RTOpT<Scalar>                &op
  ,const int                          num_vecs
  ,const ConstSubVectorView<Scalar>   sub_vecs[]
  ,const int                          num_targ_vecs
  ,const SubVectorView<Scalar>        targ_sub_vecs[]
  ,ReductTarget                       *reduct_obj
  );


/** \brief Apply an RTOp in SMPD mode to a set of columns to a set of
 * multi-vectors with contiguous storage per process.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void SPMD_apply_op(
  const Teuchos::Comm<index_type>          *comm
  ,const RTOpT<Scalar>                     &op
  ,const int                               num_cols
  ,const int                               num_multi_vecs
  ,const ConstSubMultiVectorView<Scalar>   sub_multi_vecs[]
  ,const int                               num_targ_multi_vecs
  ,const SubMultiVectorView<Scalar>        targ_sub_multi_vecs[]
  ,ReductTarget*const                      reduct_objs[]
  );


/** \brief Apply an RTOp in SMPD mode to a set of columns to a set of
 * multi-vectors with contiguous storage per process.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
template<class Scalar>
void  SPMD_apply_op(
  const Teuchos::Comm<index_type>           *comm
  ,const RTOpT<Scalar>                      &op
  ,const int                                num_cols
  ,const int                                num_vecs
  ,const ConstSubVectorView<Scalar>         sub_vecs[]
  ,const int                                num_targ_vecs
  ,const SubVectorView<Scalar>              sub_targ_vecs[]
  ,ReductTarget*const                       reduct_objs[]
  );


} // end namespace RTOpPack


#endif // RTOPPACK_SPMD_APPLY_OP_DECL_HPP
