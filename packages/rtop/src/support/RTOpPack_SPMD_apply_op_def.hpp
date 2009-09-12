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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_SPMD_APPLY_OP_DEF_HPP
#define RTOPPACK_SPMD_APPLY_OP_DEF_HPP

#include "RTOpPack_SPMD_apply_op_decl.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_CommHelpers.hpp"

#ifdef RTOPPACK_DEBUG
#  include "Teuchos_VerboseObject.hpp"
#endif // RTOPPACK_DEBUG


namespace RTOpPack {


#ifdef RTOPPACK_DEBUG


template<class Scalar>
void print( const ConstSubVectorView<Scalar> &v, Teuchos::FancyOStream &out_arg )
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(&out_arg,false);
  Teuchos::OSTab tab(out);
  *out << "globalOffset="<<v.globalOffset()<<"\n";
  *out << "subDim="<<v.subDim()<<"\n";
  *out << "values:\n";
  tab.incrTab();
  for( int i = 0; i < v.subDim(); ++i )
    *out << " " << v(i) << ":" << (v.globalOffset()+i);
  *out << "\n";
}

# include "Teuchos_VerboseObject.hpp"


#endif // RTOPPACK_DEBUG


} // namespace RTOpPack


// ///////////////////////////
// Template implementations


//
// Misc Helper functions
//


template<class PrimitiveScalar>
int RTOpPack::serializedSize(
  int num_values,
  int num_indexes,
  int num_chars
  )
{
  return 3 * sizeof(index_type)
    + num_values * sizeof(PrimitiveScalar)
    + num_indexes * sizeof(index_type)
    + num_chars * sizeof(char_type);
}


template<class Scalar>
void RTOpPack::serialize(
  const RTOpT<Scalar> &op,
  int num_values,
  int num_indexes,
  int num_chars,
  const ReductTarget &reduct_obj,
  char reduct_obj_ext[]
  )
{
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  typedef Teuchos::SerializationTraits<int,primitive_value_type> PVTST;
  typedef Teuchos::SerializationTraits<int,index_type> ITST;
  typedef Teuchos::SerializationTraits<int,char_type> CTST;
  const int
    prim_value_type_size = PVTST::fromCountToIndirectBytes(1),
    index_type_size = ITST::fromCountToIndirectBytes(1);
  //char_type_size = CTST::fromCountToIndirectBytes(1);
  const int
    num_values_off = 0,
    num_indexes_off = num_values_off + index_type_size,
    num_chars_off = num_indexes_off + index_type_size,
    values_off = num_chars_off + index_type_size,
    indexes_off = values_off + num_values * prim_value_type_size,
    chars_off = indexes_off + num_indexes * index_type_size;
  ITST::serialize(1,&num_values,index_type_size,&reduct_obj_ext[num_values_off]);
  ITST::serialize(1,&num_indexes,index_type_size,&reduct_obj_ext[num_indexes_off]);
  ITST::serialize(1,&num_chars,index_type_size,&reduct_obj_ext[num_chars_off]);
  op.extract_reduct_obj_state(
    reduct_obj
    ,num_values, num_values ? PVTST::convertFromCharPtr(&reduct_obj_ext[values_off]) : 0
    ,num_indexes, num_indexes ? ITST::convertFromCharPtr(&reduct_obj_ext[indexes_off]) : 0
    ,num_chars, num_chars ? CTST::convertFromCharPtr(&reduct_obj_ext[chars_off]) : 0
    );
  // ToDo: Change above implementation to only require indirect serialization!
}


template<class Scalar>
void RTOpPack::deserialize(
  const RTOpT<Scalar> &op,
  int num_values_in,
  int num_indexes_in,
  int num_chars_in,
  const char reduct_obj_ext[],
  ReductTarget *reduct_obj
  )
{
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  typedef Teuchos::SerializationTraits<int,primitive_value_type> PVTST;
  typedef Teuchos::SerializationTraits<int,index_type> ITST;
  typedef Teuchos::SerializationTraits<int,char_type> CTST;
  const int
    prim_value_type_size = PVTST::fromCountToIndirectBytes(1),
    index_type_size = ITST::fromCountToIndirectBytes(1);
  //char_type_size = CTST::fromCountToIndirectBytes(1);
  const int
    num_values_off = 0,
    num_indexes_off = num_values_off + index_type_size,
    num_chars_off = num_indexes_off + index_type_size,
    values_off = num_chars_off + index_type_size,
    indexes_off = values_off + num_values_in * prim_value_type_size,
    chars_off = indexes_off + num_indexes_in * index_type_size;
#ifdef RTOPPACK_DEBUG
  int num_values = -1, num_indexes = -1, num_chars = -1;
  ITST::deserialize(index_type_size,&reduct_obj_ext[num_values_off],1,&num_values);
  ITST::deserialize(index_type_size,&reduct_obj_ext[num_indexes_off],1,&num_indexes);
  ITST::deserialize(index_type_size,&reduct_obj_ext[num_chars_off],1,&num_chars);
  TEST_FOR_EXCEPT(
    !(
      num_values==num_values_in && num_indexes==num_indexes_in
      && num_chars==num_chars_in )
    );
#endif
  op.load_reduct_obj_state(
    num_values_in, 
    num_values_in ? PVTST::convertFromCharPtr(&reduct_obj_ext[values_off]) : 0
    ,num_indexes_in, num_indexes_in ? ITST::convertFromCharPtr(&reduct_obj_ext[indexes_off]) : 0
    ,num_chars_in, num_chars_in ? CTST::convertFromCharPtr(&reduct_obj_ext[chars_off]) : 0
    ,reduct_obj
    );
  // ToDo: Change above implementation to only require indirect serialization!
}


namespace RTOpPack {


//
// ReductTargetSerializer
//


template<class Scalar>
ReductTargetSerializer<Scalar>::ReductTargetSerializer(
  const Teuchos::RCP<const RTOpT<Scalar> > &op
  )
  :op_(op.assert_not_null())
{
  typedef typename RTOpT<Scalar>::primitive_value_type PrimitiveScalar;
  op_->get_reduct_type_num_entries(
    &num_values_,&num_indexes_,&num_chars_
    );
  reduct_obj_ext_size_
    = serializedSize<PrimitiveScalar>(num_values_,num_indexes_,num_chars_);
}


template<class Scalar>
index_type
ReductTargetSerializer<Scalar>::getBufferSize(const index_type count) const
{
  return reduct_obj_ext_size_ * count;
}


template<class Scalar>
void ReductTargetSerializer<Scalar>::serialize(
  const index_type count
  ,const ReductTarget * const reduct_objs[]
  ,const index_type bytes
  ,char charBuffer[]
  ) const
{
#ifdef RTOPPACK_DEBUG
  TEST_FOR_EXCEPT( !(count > 0) );
  TEST_FOR_EXCEPT( !reduct_objs );
  TEST_FOR_EXCEPT( !(bytes==this->getBufferSize(count)) );
  TEST_FOR_EXCEPT( !charBuffer );
#endif
  int offset = 0;
  for( int i = 0; i < count; ++i, offset += reduct_obj_ext_size_ ) {
    RTOpPack::serialize(
      *op_,num_values_,num_indexes_,num_chars_
      ,*reduct_objs[i],&charBuffer[offset]
      );
  }
}


template<class Scalar>
Teuchos::RCP<ReductTarget>
ReductTargetSerializer<Scalar>::createObj() const
{
  return op_->reduct_obj_create();
}

template<class Scalar>
void ReductTargetSerializer<Scalar>::deserialize(
  const index_type bytes
  ,const char charBuffer[]
  ,const index_type count
  ,ReductTarget * const reduct_objs[]
  ) const
{
#ifdef RTOPPACK_DEBUG
  TEST_FOR_EXCEPT( !(bytes > 0) );
  TEST_FOR_EXCEPT( !charBuffer );
  TEST_FOR_EXCEPT( !(bytes==getBufferSize(count)) );
  TEST_FOR_EXCEPT( !reduct_objs );
#endif
  int offset = 0;
  for( int i = 0; i < count; ++i, offset += reduct_obj_ext_size_ ) {
    RTOpPack::deserialize(
      *op_,num_values_,num_indexes_,num_chars_
      ,&charBuffer[offset],reduct_objs[i]
      );
  }
}


//
// ReductTargetReductionOp
//


template<class Scalar>
ReductTargetReductionOp<Scalar>::ReductTargetReductionOp(
  const Teuchos::RCP<const RTOpT<Scalar> > &op
  )
  :op_(op)
{}

 
template<class Scalar>
void ReductTargetReductionOp<Scalar>::reduce(
  const Ordinal count
  ,const ReductTarget*const inBuffer[]
  ,ReductTarget*const inoutBuffer[]
  ) const
{
  for( int i = 0; i < count; ++i )
    op_->reduce_reduct_objs( *inBuffer[i], inoutBuffer[i] );
}


} // namespace RTOpPack


template<class Scalar>
void RTOpPack::SPMD_all_reduce(
  const Teuchos::Comm<index_type> *comm,
  const RTOpT<Scalar> &op,
  const int num_cols,
  const ReductTarget*const i_reduct_objs[],
  ReductTarget*const reduct_objs[]
  )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  Workspace<Teuchos::RCP<ReductTarget> >
    i_i_reduct_objs( wss, num_cols );
  Workspace<ReductTarget*>
    _i_i_reduct_objs( wss, num_cols );
  for( int kc = 0; kc < num_cols; ++kc ) {
    i_i_reduct_objs[kc] = op.reduct_obj_create();
    _i_i_reduct_objs[kc] = &*i_i_reduct_objs[kc];
  }
  ReductTargetSerializer<Scalar>
    serializer(Teuchos::rcp(&op,false));
  ReductTargetReductionOp<Scalar>
    reductOp(Teuchos::rcp(&op,false));
  reduceAll(
    *comm,serializer,reductOp
    ,num_cols,&i_reduct_objs[0],&_i_i_reduct_objs[0]
    );
  for( int kc = 0; kc < num_cols; ++kc ) {
    op.reduce_reduct_objs(*_i_i_reduct_objs[kc],reduct_objs[kc]);
  }
}


template<class Scalar>
void RTOpPack::SPMD_apply_op(
  const Teuchos::Comm<index_type> *comm,
  const RTOpT<Scalar> &op,
  const int num_vecs,
  const RTOpPack::ConstSubVectorView<Scalar> sub_vecs[],
  const int num_targ_vecs,
  const RTOpPack::SubVectorView<Scalar> targ_sub_vecs[],
  ReductTarget *reduct_obj
  )
{
  ReductTarget* reduct_objs[] = { reduct_obj };
  SPMD_apply_op(
    comm,op,1,num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs
    ,reduct_obj ? reduct_objs : NULL
    );
}


/** \brief . */
template<class Scalar>
void RTOpPack::SPMD_apply_op(
  const Teuchos::Comm<index_type> *comm,
  const RTOpT<Scalar> &op,
  const int num_cols,
  const int num_multi_vecs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> sub_multi_vecs[],
  const int num_targ_multi_vecs,
  const RTOpPack::SubMultiVectorView<Scalar> targ_sub_multi_vecs[],
  RTOpPack::ReductTarget*const reduct_objs[]
  )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  int k, j, off;
  Workspace<ConstSubVectorView<Scalar> > c_sub_vecs(wss,num_multi_vecs*num_cols);
  if(sub_multi_vecs) {
    for( off = 0, j = 0; j < num_cols; ++j ) {
      for( k = 0; k < num_multi_vecs; ++k ) {
        const ConstSubMultiVectorView<Scalar> &mv = sub_multi_vecs[k];
        c_sub_vecs[off++].initialize(mv.globalOffset(),mv.subDim(),&mv(0,j),1);
      }
    }
  }
  Workspace<SubVectorView<Scalar> > c_targ_sub_vecs(wss,num_targ_multi_vecs*num_cols);
  if(targ_sub_multi_vecs) {
    for( off = 0, j = 0; j < num_cols; ++j ) {
      for( k = 0; k < num_targ_multi_vecs; ++k ) {
        const SubMultiVectorView<Scalar> &mv = targ_sub_multi_vecs[k];
        c_targ_sub_vecs[off++].initialize(mv.globalOffset(),mv.subDim(),&mv(0,j),1);
      }
    }
  }
  SPMD_apply_op(
    comm,op,num_cols
    ,num_multi_vecs, num_multi_vecs && sub_multi_vecs ? &c_sub_vecs[0] : NULL
    ,num_targ_multi_vecs, num_targ_multi_vecs && targ_sub_multi_vecs ? &c_targ_sub_vecs[0] : NULL
    ,reduct_objs
    );
}


template<class Scalar>
void RTOpPack::SPMD_apply_op(
  const Teuchos::Comm<index_type> *comm,
  const RTOpT<Scalar> &op,
  const int num_cols,
  const int num_vecs,
  const ConstSubVectorView<Scalar> sub_vecs[],
  const int num_targ_vecs,
  const SubVectorView<Scalar> sub_targ_vecs[],
  ReductTarget*const reduct_objs[]
  )
{
#ifdef RTOPPACK_DEBUG
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(show_spmd_apply_op_dump) {
    *out << "\nEntering RTOpPack::SPMD_apply_op(...) ...\n";
    *out
      << "\ncomm = " << (comm?comm->description():"NULL")
      << "\nop = " << op.description()
      << "\nnum_cols = " << num_cols
      << "\nnum_vecs = " << num_vecs
      << "\nnum_targ_vecs = " << num_targ_vecs
      << "\n";
    if( num_vecs && sub_vecs ) {
      *out << "\nInput vectors:\n";
      Teuchos::OSTab tab(out);
      for( int kc = 0; kc < num_cols; ++kc ) {
        for( int k = 0; k < num_vecs; ++k ) {
          *out << "\nvecs["<<kc<<","<<k<<"] =\n";
          print(sub_vecs[kc*num_vecs+k],*out);
        }
      }
    }
    if( num_targ_vecs && sub_targ_vecs ) {
      *out << "\nInput/output vectors *before* transforamtion:\n";
      Teuchos::OSTab tab(out);
      for( int kc = 0; kc < num_cols; ++kc ) {
        for( int k = 0; k < num_targ_vecs; ++k ) {
          *out << "\nvecs["<<kc<<","<<k<<"] =\n";
          print(sub_targ_vecs[kc*num_targ_vecs+k],*out);
        }
      }
    }
    if(reduct_objs) {
      *out << "\nInput/output reduction objects *before* reduction:\n";
      Teuchos::OSTab tab(out);
      for( int kc = 0; kc < num_cols; ++kc ) {
        *out
          << "\nreduct_objs["<<kc<<"] =\n"
          << describe(*reduct_objs[kc],Teuchos::VERB_EXTREME);
      }
    }
  }
#endif // RTOPPACK_DEBUG
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  if( reduct_objs == NULL && sub_vecs == NULL && sub_targ_vecs == NULL ) {
    // This is a transformation operation with no data on this processor.
    // Therefore, we can just exist!
  }
  else {
    const int localSubDim =
      ( num_vecs
        ? ( sub_vecs ? sub_vecs[0].subDim() : 0 )
        : ( sub_targ_vecs ? sub_targ_vecs[0].subDim() : 0 )
        );
    // See if we need to do any global communication at all?
    if( comm==NULL || reduct_objs == NULL ) {
      if( ( sub_vecs || sub_targ_vecs ) && localSubDim ) {
        for( int kc = 0; kc < num_cols; ++kc ) {
          op.apply_op(
            num_vecs,sub_vecs+kc*num_vecs,num_targ_vecs,sub_targ_vecs+kc*num_targ_vecs
            ,reduct_objs ? reduct_objs[kc] : NULL
            );
        }
      }
    }
    else {
      // Check the preconditions for excluding empty target vectors.
      TEST_FOR_EXCEPTION(
        ( ( num_vecs && !sub_vecs) || ( num_targ_vecs && !sub_targ_vecs) ) && !( !sub_vecs && !sub_targ_vecs )
        ,std::logic_error
        ,"SPMD_apply_op(...): Error, invalid arguments num_vecs = " << num_vecs
        << ", sub_vecs = " << sub_vecs << ", num_targ_vecs = " << num_targ_vecs
        << ", sub_targ_vecs = " << sub_targ_vecs
        );
      //
      // There is a non-null reduction target object and we are using
      // SPMD so we need to reduce it across processors
      //
      // Allocate the intermediate target object and perform the
      // reduction for the vector elements on this processor.
      //
      Workspace<Teuchos::RCP<ReductTarget> >
        i_reduct_objs( wss, num_cols );
      for( int kc = 0; kc < num_cols; ++kc ) {
        i_reduct_objs[kc] = op.reduct_obj_create();
        if( ( sub_vecs || sub_targ_vecs ) && localSubDim ) {
          op.apply_op(
            num_vecs, sub_vecs+kc*num_vecs, num_targ_vecs, sub_targ_vecs+kc*num_targ_vecs
            ,&*i_reduct_objs[kc]
            );
        }
      }
#ifdef RTOPPACK_DEBUG
      if(show_spmd_apply_op_dump) {
        if(reduct_objs) {
          *out << "\nIntermediate reduction objects in this process before global reduction:\n";
          Teuchos::OSTab tab(out);
          for( int kc = 0; kc < num_cols; ++kc ) {
            *out
              << "\ni_reduct_objs["<<kc<<"] =\n"
              << describe(*i_reduct_objs[kc],Teuchos::VERB_EXTREME);
          }
        }
      }
#endif // RTOPPACK_DEBUG
      //
      // Reduce the local intermediate reduction objects into the global reduction objects
      //
      Workspace<const ReductTarget*>
        _i_reduct_objs( wss, num_cols );
      for( int kc = 0; kc < num_cols; ++kc ) {
        _i_reduct_objs[kc] = &*i_reduct_objs[kc];
      }
#ifdef RTOPPACK_DEBUG
      if(show_spmd_apply_op_dump) {
        if(reduct_objs) {
          *out << "\nPerforming global reduction ...\n";
        }
      }
#endif // RTOPPACK_DEBUG
      SPMD_all_reduce(comm,op,num_cols,&_i_reduct_objs[0],reduct_objs);
    }
  }
#ifdef RTOPPACK_DEBUG
  if(show_spmd_apply_op_dump) {
    if( num_targ_vecs && sub_targ_vecs ) {
      *out << "\nInput/output vectors *after* transforamtion:\n";
      Teuchos::OSTab tab(out);
      for( int kc = 0; kc < num_cols; ++kc ) {
        for( int k = 0; k < num_targ_vecs; ++k ) {
          *out << "\nvecs["<<kc<<","<<k<<"] =\n";
          print(sub_targ_vecs[kc*num_targ_vecs+k],*out);
        }
      }
    }
    if(reduct_objs) {
      *out << "\nInput/output reduction objects *after* reduction:\n";
      Teuchos::OSTab tab(out);
      for( int kc = 0; kc < num_cols; ++kc ) {
        *out
          << "\nreduct_objs["<<kc<<"] =\n"
          << describe(*reduct_objs[kc],Teuchos::VERB_EXTREME);
      }
    }
    *out << "\nLeaving RTOpPack::SPMD_apply_op(...) ...\n";
  }
#endif // RTOPPACK_DEBUG
}


//
// Explicit Template Instaniation Macros
//


#define RTOPPACK_SPMD_APPLY_OP_INSTANT_SCALAR(SCALAR) \
  \
  template int serializedSize<SCALAR >( \
    int num_values, \
    int num_indexes, \
    int num_chars \
    ); \
  \
  template void serialize<SCALAR >( \
    const RTOpT<SCALAR > &op, \
    int num_values, \
    int num_indexes, \
    int num_chars, \
    const ReductTarget &reduct_obj, \
    char reduct_obj_ext[] \
    ); \
  \
  template void deserialize<SCALAR >( \
    const RTOpT<SCALAR > &op, \
    int num_values_in, \
    int num_indexes_in, \
    int num_chars_in, \
    const char reduct_obj_ext[], \
    ReductTarget *reduct_obj \
    ); \
  \
  template class ReductTargetSerializer<SCALAR >; \
  \
  template class ReductTargetReductionOp<SCALAR >; \
  \
  template void SPMD_all_reduce<SCALAR >( \
    const Teuchos::Comm<index_type> *comm, \
    const RTOpT<SCALAR > &op, \
    const int num_cols, \
    const ReductTarget*const i_reduct_objs[], \
    ReductTarget*const reduct_objs[] \
    ); \
  \
  template void SPMD_apply_op<SCALAR >( \
    const Teuchos::Comm<index_type> *comm, \
    const RTOpT<SCALAR > &op, \
    const int num_vecs, \
    const RTOpPack::ConstSubVectorView<SCALAR > sub_vecs[], \
    const int num_targ_vecs, \
    const RTOpPack::SubVectorView<SCALAR > targ_sub_vecs[], \
    ReductTarget *reduct_obj \
    ); \
  \
  template void SPMD_apply_op<SCALAR >( \
    const Teuchos::Comm<index_type> *comm, \
    const RTOpT<SCALAR > &op, \
    const int num_cols, \
    const int num_multi_vecs, \
    const RTOpPack::ConstSubMultiVectorView<SCALAR > sub_multi_vecs[], \
    const int num_targ_multi_vecs, \
    const RTOpPack::SubMultiVectorView<SCALAR > targ_sub_multi_vecs[], \
    RTOpPack::ReductTarget*const reduct_objs[] \
    ); \
  \
  template void SPMD_apply_op<SCALAR >( \
    const Teuchos::Comm<index_type> *comm, \
    const RTOpT<SCALAR > &op, \
    const int num_cols, \
    const int num_vecs, \
    const ConstSubVectorView<SCALAR > sub_vecs[], \
    const int num_targ_vecs, \
    const SubVectorView<SCALAR > sub_targ_vecs[], \
    ReductTarget*const reduct_objs[] \
    );


#endif // RTOPPACK_SPMD_APPLY_OP_DEF_HPP
