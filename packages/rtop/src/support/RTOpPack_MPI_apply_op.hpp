// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_MPI_APPLY_OP_HPP
#define RTOPPACK_MPI_APPLY_OP_HPP

#include "RTOpPack_MPI_apply_op_decl.hpp"
#include "Teuchos_RawMPITraits.hpp"
#include "Teuchos_Workspace.hpp"
#ifdef RTOp_USE_MPI
//#  include "Teuchos_MPIReductionOpBase.hpp"
#  include "Teuchos_MpiReductionOpSetter.hpp"
#endif

#ifdef RTOPPACK_MPI_APPLY_OP_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif // RTOPPACK_MPI_APPLY_OP_DUMP

namespace RTOpPack {

#ifdef RTOPPACK_MPI_APPLY_OP_DUMP
template<class Scalar>
void print( const ConstSubVectorView<Scalar> &v, Teuchos::FancyOStream &out_arg )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::rcp(&out_arg,false);
  Teuchos::OSTab tab(out);
  *out << "globalOffset="<<v.globalOffset()<<"\n";
  *out << "subDim="<<v.subDim()<<"\n";
  *out << "values:\n";
  tab.incrTab();
  for( int i = 0; i < v.subDim(); ++i )
    *out << " " << v(i) << ":" << (v.globalOffset()+i);
  *out << "\n";
}

#  include "Teuchos_VerboseObject.hpp"
#endif // RTOPPACK_MPI_APPLY_OP_DUMP

#ifdef RTOp_USE_MPI

template<class Scalar>
class MpiReductionOp : public Teuchos::MpiReductionOpBase {
public:

  MpiReductionOp(
    const Teuchos::RefCountPtr<const RTOpT<Scalar> >  &op
    )
    : op_(op)
    {
      *op; // Assert not null
    }

  void reduce(
    void              *invec
    ,void             *inoutvec
    ,int              *len
    ,RTOp_Datatype    *datatype
    ) const
    {
      typedef typename RTOpT<Scalar>::primitive_value_type  primitive_value_type;
      int num_reduct_type_values = 0, num_reduct_type_indexes = 0, num_reduct_type_chars = 0;
      op_->get_reduct_type_num_entries(
        &num_reduct_type_values, &num_reduct_type_indexes, &num_reduct_type_chars
        );
      const int reduct_obj_ext_size
        = RTOpPack::reduct_obj_ext_size<primitive_value_type>(
          num_reduct_type_values,num_reduct_type_indexes,num_reduct_type_chars
          );
      const char
        *in_reduct_obj_ext = reinterpret_cast<char*>(invec);
      char
        *inout_reduct_obj_ext = reinterpret_cast<char*>(inoutvec);
      Teuchos::RefCountPtr<ReductTarget>
        in_reduct_obj    = op_->reduct_obj_create(),
        inout_reduct_obj = op_->reduct_obj_create();
      for( int i = 0; i < (*len); ++i, in_reduct_obj_ext += reduct_obj_ext_size, inout_reduct_obj_ext += reduct_obj_ext_size  )
      {
        load_reduct_obj_ext_state( *op_, in_reduct_obj_ext,    &*in_reduct_obj    );
        load_reduct_obj_ext_state( *op_, inout_reduct_obj_ext, &*inout_reduct_obj );
        op_->reduce_reduct_objs( *in_reduct_obj, &*inout_reduct_obj );
        extract_reduct_obj_ext_state(
          *op_,*inout_reduct_obj,num_reduct_type_values,num_reduct_type_indexes,num_reduct_type_chars
          ,inout_reduct_obj_ext
          );
      }
    }
private:
  Teuchos::RefCountPtr<const RTOpT<Scalar> >  op_;
  // Not defined and not to be called!
  MpiReductionOp();
  MpiReductionOp(const MpiReductionOp<Scalar>&);
  MpiReductionOp<Scalar>& operator=(const MpiReductionOp<Scalar>&);
};


#endif // RTOp_USE_MPI

} // namespace RTOpPack

//
// Template implementation stuff
//

namespace {

// This class is designed to call MPI functions to free
// opaque objects.  This allows robust behavior in the
// face of exceptions.
template<class T>
class call_free_func {
public:
//  typedef extern "C" int (*free_func_ptr_t)(T*);
    typedef int (*free_func_ptr_t)(T*);
  call_free_func(T* opaque_obj, T null_value, free_func_ptr_t free_func_ptr)
    : opaque_obj_(opaque_obj), null_value_(null_value), free_func_ptr_(free_func_ptr)
    {}
  ~call_free_func()
    {
      if(*opaque_obj_ != null_value_)
        free_func_ptr_(opaque_obj_);
    }
private:
  T*                opaque_obj_;
  T                 null_value_;
  free_func_ptr_t   free_func_ptr_;
  // not defined and not to be called
  call_free_func();
  call_free_func(const call_free_func&);
  call_free_func& operator=(const call_free_func&);
}; 

} // end namespace

template<class primitive_value_type>
void RTOpPack::MPI_type_signature(
  const int num_values
  ,const int num_indexes
  ,const int num_chars
  ,int* num_entries
  ,int block_lengths[]
  ,MPI_Aint displacements[]
  ,MPI_Datatype datatypes[]
  )
{
  int k = 0, off = 0;
  // values
  block_lengths[k] = 3 + num_values; /* must carry size information */
  displacements[k] = 0;
  datatypes[k]     = Teuchos::RawMPITraits<primitive_value_type>::type();
  ++k;
  off += Teuchos::RawMPITraits<primitive_value_type>::adjustCount(3 + num_values) * sizeof(primitive_value_type);
  // indexes
  if( num_indexes ) {
    block_lengths[k] = num_indexes;
    displacements[k] = off;
    datatypes[k]     = Teuchos::RawMPITraits<index_type>::type();
    ++k;
    off += num_indexes * sizeof(index_type);
  }
  // chars
  if( num_chars ) {
    block_lengths[k] = num_chars;
    displacements[k] = off;
    datatypes[k]     = Teuchos::RawMPITraits<char_type>::type();
    ++k;
  }
  *num_entries = k;
}

template<class primitive_value_type>
int reduct_obj_ext_size(
  int   num_values
  ,int  num_indexes
  ,int  num_chars
  )
{
  TEST_FOR_EXCEPT(true);
  return 0;
}

template<class Scalar>
void RTOpPack::extract_reduct_obj_ext_state(
  const RTOpT<Scalar>    &op
  ,const ReductTarget    &reduct_obj
  ,int                   num_values
  ,int                   num_indexes
  ,int                   num_chars
  ,void                  *_reduct_obj_ext
  )
{
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  char *reduct_obj_ext = reinterpret_cast<char*>(_reduct_obj_ext);
  const int
    num_values_off  = 0,
    num_indexes_off = num_values_off  + sizeof(primitive_value_type),
    num_chars_off   = num_indexes_off + sizeof(primitive_value_type),
    values_off      = num_chars_off   + sizeof(primitive_value_type),
    indexes_off     = values_off      + num_values  * sizeof(primitive_value_type),
    chars_off       = indexes_off     + num_indexes * sizeof(index_type);
  *reinterpret_cast<primitive_value_type*>(reduct_obj_ext + num_values_off)  = static_cast<primitive_value_type>(num_values);
  *reinterpret_cast<primitive_value_type*>(reduct_obj_ext + num_indexes_off) = static_cast<primitive_value_type>(num_indexes);
  *reinterpret_cast<primitive_value_type*>(reduct_obj_ext + num_chars_off)   = static_cast<primitive_value_type>(num_chars);
  op.extract_reduct_obj_state(
    reduct_obj
    ,num_values,  num_values  ? reinterpret_cast<primitive_value_type*>(reduct_obj_ext + values_off) : NULL
     ,num_indexes, num_indexes ? reinterpret_cast<index_type*>(reduct_obj_ext + indexes_off)          : NULL
    ,num_chars,   num_chars   ? reinterpret_cast<char_type*>(reduct_obj_ext  + chars_off)            : NULL
    );
}

template<class Scalar>
void RTOpPack::load_reduct_obj_ext_state(
  const RTOpT<Scalar>    &op
  ,const void            *_reduct_obj_ext
  ,ReductTarget          *reduct_obj
  )
{
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  const char *reduct_obj_ext = reinterpret_cast<const char*>(_reduct_obj_ext);
  const int
    num_values_off  = 0,
    num_values      =                   static_cast<int>(*reinterpret_cast<const primitive_value_type*>(reduct_obj_ext + num_values_off)),
    num_indexes_off = num_values_off  + sizeof(primitive_value_type),
    num_indexes     =                   static_cast<int>(*reinterpret_cast<const primitive_value_type*>(reduct_obj_ext + num_indexes_off)),
    num_chars_off   = num_indexes_off + sizeof(primitive_value_type),
    num_chars       =                   static_cast<int>(*reinterpret_cast<const primitive_value_type*>(reduct_obj_ext + num_chars_off)),
    values_off      = num_chars_off + sizeof(primitive_value_type),
    indexes_off     = values_off  + num_values  * sizeof(primitive_value_type),
    chars_off       = indexes_off + num_indexes * sizeof(index_type);
  op.load_reduct_obj_state(
    num_values,   num_values  ? reinterpret_cast<const primitive_value_type*>(reduct_obj_ext + values_off) : NULL
     ,num_indexes, num_indexes ? reinterpret_cast<const index_type*>(reduct_obj_ext + indexes_off)          : NULL
    ,num_chars,   num_chars   ? reinterpret_cast<const char_type*>(reduct_obj_ext  + chars_off)            : NULL
    ,reduct_obj
    );
}

template<class Scalar>
void RTOpPack::MPI_apply_op(
  MPI_Comm                                      comm
  ,const RTOpT<Scalar>                          &op
  ,const int                                    root_rank
  ,const int                                    num_vecs
  ,const RTOpPack::ConstSubVectorView<Scalar>   sub_vecs[]
  ,const int                                    num_targ_vecs
  ,const RTOpPack::SubVectorView<Scalar>        targ_sub_vecs[]
  ,ReductTarget                                 *reduct_obj
  )
{
  ReductTarget* reduct_objs[] = { reduct_obj };
  MPI_apply_op(
    comm,op,root_rank,1,num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs
    ,reduct_obj ? reduct_objs : NULL
    );
}

/** \brief . */
template<class Scalar>
void RTOpPack::MPI_apply_op(
  MPI_Comm                                           comm
  ,const RTOpT<Scalar>                               &op
  ,const int                                         root_rank
  ,const int                                         num_cols
  ,const int                                         num_multi_vecs
  ,const RTOpPack::ConstSubMultiVectorView<Scalar>   sub_multi_vecs[]
  ,const int                                         num_targ_multi_vecs
  ,const RTOpPack::SubMultiVectorView<Scalar>        targ_sub_multi_vecs[]
  ,RTOpPack::ReductTarget*const                      reduct_objs[]
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
  MPI_apply_op(
    comm,op,root_rank,num_cols
    ,num_multi_vecs, num_multi_vecs && sub_multi_vecs ? &c_sub_vecs[0] : NULL
    ,num_targ_multi_vecs, num_targ_multi_vecs && targ_sub_multi_vecs ? &c_targ_sub_vecs[0] : NULL
    ,reduct_objs
    );
}

template<class Scalar>
void RTOpPack::MPI_all_reduce(
  MPI_Comm                            comm
  ,const RTOpT<Scalar>                &op
  ,const int                          root_rank
  ,const int                          num_cols
  ,const ReductTarget*const           i_reduct_objs[]
  ,ReductTarget*const                 reduct_objs[]
  )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  //
  // Get the description of the datatype of the target object.  We
  // need this in order to map it to an MPI user-defined data type so
  // that the reduction operations can be performed over all of the of
  // processors.
  //
  // Get the number of the entries in the array that describes target
  // type's datatypes
  int num_reduct_type_values = 0, num_reduct_type_indexes = 0,
    num_reduct_type_chars = 0,  num_reduct_type_entries = 0;
  op.get_reduct_type_num_entries(
    &num_reduct_type_values, &num_reduct_type_indexes, &num_reduct_type_chars
    );
  num_reduct_type_entries
    = (num_reduct_type_values ? 1 : 0)
    + (num_reduct_type_indexes ? 1 : 0)
    + (num_reduct_type_chars ? 1 : 0);
  //
  // Get the arrays that describe the reduction target type datatypes
  //
  if( num_reduct_type_values == 0 ) ++num_reduct_type_entries; // Add the block for the sizes
  int
    target_type_block_lengths[3];
  MPI_Aint
    target_type_displacements[3];
  RTOp_Datatype
    target_type_datatypes[3];
  MPI_type_signature<primitive_value_type>(
    num_reduct_type_values, num_reduct_type_indexes, num_reduct_type_chars
    ,&num_reduct_type_entries
    ,target_type_block_lengths, target_type_displacements, target_type_datatypes 
    );
  //
  // Translate this reduction target datatype description to an MPI datatype
  //
  MPI_Datatype
    mpi_reduct_ext_type = MPI_DATATYPE_NULL;
  TEST_FOR_EXCEPTION(
    0!=MPI_Type_struct(
      num_reduct_type_entries
      ,target_type_block_lengths, target_type_displacements
      ,target_type_datatypes, &mpi_reduct_ext_type
      )
    ,std::logic_error, "Error!"
    );
  call_free_func<MPI_Datatype>
    free_mpi_reduct_ext_type(&mpi_reduct_ext_type,MPI_DATATYPE_NULL,MPI_Type_free);
  TEST_FOR_EXCEPTION(
    0!=MPI_Type_commit( &mpi_reduct_ext_type )
    ,std::logic_error,"Error!"
    );
  //
  // Create external contiguous representations for the intermediate
  // reduction objects that MPI can use.  This is very low level but
  // at least the user does not have to do it.
  //
  const int reduct_obj_ext_size
    = RTOpPack::reduct_obj_ext_size<primitive_value_type>(
      num_reduct_type_values,num_reduct_type_indexes,num_reduct_type_chars
      );
  Workspace<char> // This is raw memory for external intermediate reduction object.
    i_reduct_objs_ext( wss, reduct_obj_ext_size*num_cols, false );
  for( int kc = 0; kc < num_cols; ++kc ) {
    extract_reduct_obj_ext_state(
      op,*i_reduct_objs[kc],num_reduct_type_values,num_reduct_type_indexes,num_reduct_type_chars
      ,&i_reduct_objs_ext[0]+kc*reduct_obj_ext_size
      );
  }
  //
  // Now perform the global reduction over all of the processors.
  //
  if( root_rank >= 0 ) {
    TEST_FOR_EXCEPTION( true, std::logic_error, "Error, we don't support root_rank >= 0 yet!" );
  }
  else {
    // Apply the reduction operation over all of the processors
    // and reduce to one target object on each processor
    //
    // Create output argument for MPI call (must be different from input)
    //
    Workspace<char> // This is raw memory for external intermediate reduction object.
      i_reduct_objs_tmp( wss, reduct_obj_ext_size*num_cols, false );
    Teuchos::RefCountPtr<ReductTarget> i_reduct_obj_tmp_uninit = op.reduct_obj_create(); // Just an unreduced object!
    for( int kc = 0; kc < num_cols; ++kc ) {
      extract_reduct_obj_ext_state( // Initialize to no initial reduction
        op,*i_reduct_obj_tmp_uninit,num_reduct_type_values,num_reduct_type_indexes,num_reduct_type_chars
        ,&i_reduct_objs_tmp[0]+kc*reduct_obj_ext_size
        );
    }
#ifdef RTOp_USE_MPI
    //
    // Setup the mpi-compatible global function object for this
    // reduction operator and create the MPI_Op object.
    //
    Teuchos::MpiReductionOpSetter
      mpi_op_setter( Teuchos::rcp(new MpiReductionOp<Scalar>(Teuchos::rcp(&op,false))) ); // Destructor will release everything!
    //
    // Call MPI
    //
    TEST_FOR_EXCEPT(
      0!=MPI_Allreduce(
        &i_reduct_objs_ext[0], &i_reduct_objs_tmp[0]
        ,num_cols, mpi_reduct_ext_type
        ,mpi_op_setter.mpi_op()
        ,comm
        )
      );
#else // RTOp_USE_MPI
    for( int k = 0; k < reduct_obj_ext_size; ++k )
      i_reduct_objs_tmp[k] = i_reduct_objs_ext[k];
#endif // RTOp_USE_MPI
    // Load the updated state of the reduction target object and reduce.
    if(1) {
      Teuchos::RefCountPtr<ReductTarget> tmp_reduct_obj = op.reduct_obj_create();
      for( int kc = 0; kc < num_cols; ++kc ) {
        load_reduct_obj_ext_state( op, &i_reduct_objs_tmp[0]+kc*reduct_obj_ext_size, &*tmp_reduct_obj );
        op.reduce_reduct_objs( *tmp_reduct_obj, reduct_objs[kc] );
      }
    }
  }
}

template<class Scalar>
void RTOpPack::MPI_apply_op(
  MPI_Comm                                  comm
  ,const RTOpT<Scalar>                      &op
  ,const int                                root_rank
  ,const int                                num_cols
  ,const int                                num_vecs
  ,const ConstSubVectorView<Scalar>         sub_vecs[]
  ,const int                                num_targ_vecs
  ,const SubVectorView<Scalar>              sub_targ_vecs[]
  ,ReductTarget*const                       reduct_objs[]
  )
{
#ifdef RTOPPACK_MPI_APPLY_OP_DUMP
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(show_mpi_apply_op_dump) {
    *out << "\nEntering RTOpPack::MPI_apply_op(...) ...\n";
    *out
      << "\ncomm = " << comm
      << "\nop = " << op.description()
      << "\nroot_rank = " << root_rank
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
#endif // RTOPPACK_MPI_APPLY_OP_DUMP
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
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
    if( comm == MPI_COMM_NULL || reduct_objs == NULL ) {
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
        ,"MPI_apply_op(...): Error, invalid arguments num_vecs = " << num_vecs
        << ", sub_vecs = " << sub_vecs << ", num_targ_vecs = " << num_targ_vecs
        << ", sub_targ_vecs = " << sub_targ_vecs
        );
      //
      // There is a non-null reduction target object and we are using
      // MPI so we need to reduce it across processors
      //
      // Allocate the intermediate target object and perform the
      // reduction for the vector elements on this processor.
      //
      Workspace<Teuchos::RefCountPtr<ReductTarget> >
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
      //
      // Reduce the local intermediate reduction objects into the global reduction objects
      //
      Workspace<const ReductTarget*>
        _i_reduct_objs( wss, num_cols );
      for( int kc = 0; kc < num_cols; ++kc )
        _i_reduct_objs[kc] = &*i_reduct_objs[kc];
      MPI_all_reduce(comm,op,root_rank,num_cols,&_i_reduct_objs[0],reduct_objs);
    }
  }
#ifdef RTOPPACK_MPI_APPLY_OP_DUMP
  if(show_mpi_apply_op_dump) {
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
    *out << "\nLeaving RTOpPack::MPI_apply_op(...) ...\n";
  }
#endif // RTOPPACK_MPI_APPLY_OP_DUMP
}

#endif // RTOPPACK_MPI_APPLY_OP_HPP
