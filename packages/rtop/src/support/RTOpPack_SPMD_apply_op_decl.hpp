// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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


/** \brief Set up to show a dump of RTOps applied through SPMD_apply_op().
 *
 * \param dumOut [in] RCP to output stream.  If non-null, output will be
 * dumped to this stream.  If null, then no output will be dumped.
 *
 * \ingroup RTOpPack_parallel_helpers_grp
 */
void set_SPMD_apply_op_dump_out(const RCP<FancyOStream> &dumpOut);


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
  : public Teuchos::ReferenceTypeReductionOp<Teuchos_Ordinal,ReductTarget>
{
public:
  /** \brief . */
  typedef Teuchos_Ordinal Ordinal;
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
  ReductTargetReductionOp(const ReductTargetReductionOp<Scalar>&);
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
