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

#ifndef RTOPPACK_ROP_GET_SUB_VECTOR_HPP
#define RTOPPACK_ROP_GET_SUB_VECTOR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief . */
template<class Scalar> class ReductTargetSubVectorT;

/** \brief Reduction operator that extracts a sub-vector in the range of
 * global zero-based indexes [l,u].
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class ROpGetSubVector : public RTOpT<Scalar> {
public:

  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;

  /** \brief . */
  ROpGetSubVector( const index_type l = 0, const index_type u = 0 );

  /// Set the range of global indexes to extract elements for
  void set_range( const index_type l, const index_type u );

  /// Extract the subvector after all of the reductions are completed
  ReductTargetSubVectorT<Scalar>& operator()( ReductTarget& reduct_obj ) const;

  /// Extract the subvector after all of the reductions are completed
  const ReductTargetSubVectorT<Scalar>& operator()( const ReductTarget& reduct_obj ) const;

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const;
  /** \brief . */
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const;
  /** \brief . */
  void reduce_reduct_objs(
    const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
    ) const;
  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const;
  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,Teuchos_Index            index_data[]
    ,int                      num_chars
    ,::RTOpPack::char_type           char_data[]
    ) const;
  /** \brief . */
  void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const Teuchos_Index           index_data[]
    ,int                           num_chars
    ,const ::RTOpPack::char_type          char_data[]
    ,ReductTarget                  *reduct_obj
    ) const;
  /** \brief . */
  void get_op_type_num_entries(
    int*  num_values
    ,int* num_indexes
    ,int* num_chars
    ) const;
  /** \brief . */
  void extract_op_state(
    int                             num_values
    ,primitive_value_type           value_data[]
    ,int                            num_indexes
    ,Teuchos_Index                  index_data[]
    ,int                            num_chars
    ,::RTOpPack::char_type                 char_data[]
    ) const;
  /** \brief . */
  void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const Teuchos_Index          index_data[]
    ,int                          num_chars
    ,const ::RTOpPack::char_type         char_data[]
    );
  /** \brief . */
  bool coord_invariant() const;
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const;

  //@}

private:

  index_type  l_;
  index_type  u_;

}; // class ROpGetSubVector

/** \brief . */
template<class Scalar>
class ReductTargetSubVectorT : public ReductTarget {
public:
  /** \brief . */
  ReductTargetSubVectorT( const index_type l, const index_type u)
    {
      Scalar *values = NULL;
      try {
        const int subDim = u-l+1;
        sub_vec.initialize(
          l                    // global_offset
          ,subDim              // subDim
          ,new Scalar[subDim]  // values
          ,1                   // stride
          );
        reinit();
      }
      catch(...) {
        delete [] values;
      }
    }
  /** \brief . */
  ~ReductTargetSubVectorT()
    {
      free(&sub_vec);
    }
  /** \brief . */
  void reinit()
    {
      std::fill_n( const_cast<Scalar*>(sub_vec.values()), sub_vec.subDim(), Teuchos::ScalarTraits<Scalar>::zero() );
    }
  /** \brief . */
  void transfer( ConstSubVectorView<Scalar> *sub_vec_out )
    {
      sub_vec_out->initialize(sub_vec.globalOffset(),sub_vec.subDim(),sub_vec.values(),sub_vec.stride());
      sub_vec.set_uninitialized();
    }
  /** \brief . */
  static void free( ConstSubVectorView<Scalar> *sub_vec )
    {
      if(sub_vec->values())
        delete [] const_cast<Scalar*>(sub_vec->values());
      sub_vec->set_uninitialized();
    }
  /** \brief . */
  ConstSubVectorView<Scalar>  sub_vec;
private:
  // Not defined and not to be called!
  ReductTargetSubVectorT();
  ReductTargetSubVectorT(const ReductTargetSubVectorT&);
  ReductTargetSubVectorT<Scalar>& operator=(const ReductTargetSubVectorT&);
};

// ////////////////////////////////
// Template definitions

template<class Scalar>
ROpGetSubVector<Scalar>::ROpGetSubVector( const index_type l, const index_type u )
  :RTOpT<Scalar>("ROpGetSubVector"), l_(l), u_(u)
{}

template<class Scalar>
void ROpGetSubVector<Scalar>::set_range( const index_type l, const index_type u )
{
  l_ = l;
  u_ = u;
}

template<class Scalar>
ReductTargetSubVectorT<Scalar>&
ROpGetSubVector<Scalar>::operator()(ReductTarget& reduct_obj ) const
{
  using Teuchos::dyn_cast;
  return dyn_cast<ReductTargetSubVectorT<Scalar> >(reduct_obj);
}

template<class Scalar>
const ReductTargetSubVectorT<Scalar>&
ROpGetSubVector<Scalar>::operator()( const ReductTarget& reduct_obj ) const
{
  using Teuchos::dyn_cast;
  return dyn_cast<const ReductTargetSubVectorT<Scalar> >(reduct_obj);
}

// Overridden from RTOpT

template<class Scalar>
void ROpGetSubVector<Scalar>::get_reduct_type_num_entries(
  int*   num_values
  ,int*  num_indexes
  ,int*  num_chars
  ) const
{
  const int num_prim_objs_per_scalar = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
  *num_values  = (u_-l_+1)*num_prim_objs_per_scalar;
  *num_indexes = 0;
  *num_chars   = 0;
}

template<class Scalar>
Teuchos::RefCountPtr<ReductTarget>
ROpGetSubVector<Scalar>::reduct_obj_create() const
{
  return Teuchos::rcp(new ReductTargetSubVectorT<Scalar>(l_,u_));
}

template<class Scalar>
void ROpGetSubVector<Scalar>::reduce_reduct_objs(
  const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
  ) const
{
  using Teuchos::dyn_cast;
  const ConstSubVectorView<Scalar> &sub_vec_in = dyn_cast<const ReductTargetSubVectorT<Scalar> >(in_reduct_obj).sub_vec;
  ConstSubVectorView<Scalar> &sub_vec_inout = dyn_cast<ReductTargetSubVectorT<Scalar> >(*inout_reduct_obj).sub_vec;
  TEST_FOR_EXCEPT(
    sub_vec_in.subDim()!=sub_vec_inout.subDim()||sub_vec_in.globalOffset()!=sub_vec_inout.globalOffset()
    ||!sub_vec_in.values()||!sub_vec_inout.values()
    );
  Scalar *svio_values = const_cast<Scalar*>(sub_vec_inout.values());
  for( int k = 0; k < sub_vec_in.subDim(); ++k ) {
    svio_values[k] += sub_vec_in[k];
  }
}

template<class Scalar>
void ROpGetSubVector<Scalar>::reduct_obj_reinit( ReductTarget* reduct_obj ) const
{
  using Teuchos::dyn_cast;
  dyn_cast<ReductTargetSubVectorT<Scalar> >(*reduct_obj).reinit();
}

template<class Scalar>
void ROpGetSubVector<Scalar>::extract_reduct_obj_state(
  const ReductTarget        &reduct_obj
  ,int                      num_values
  ,primitive_value_type     value_data[]
  ,int                      num_indexes
  ,Teuchos_Index          index_data[]
  ,int                      num_chars
  ,::RTOpPack::char_type           char_data[]
  ) const
{
  using Teuchos::dyn_cast;
  typedef Teuchos::PrimitiveTypeTraits<Scalar> PTT;
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  const ConstSubVectorView<Scalar> &sub_vec = dyn_cast<const ReductTargetSubVectorT<Scalar> >(reduct_obj).sub_vec;
  int value_data_off = 0;
  for( int k = 0; k < sub_vec.subDim(); ++k, value_data_off += num_prim_objs_per_scalar )
    PTT::extractPrimitiveObjs( sub_vec[k], num_prim_objs_per_scalar, value_data+value_data_off );
}

template<class Scalar>
void ROpGetSubVector<Scalar>::load_reduct_obj_state(
  int                            num_values
  ,const primitive_value_type    value_data[]
  ,int                           num_indexes
  ,const Teuchos_Index         index_data[]
  ,int                           num_chars
  ,const ::RTOpPack::char_type          char_data[]
  ,ReductTarget                  *reduct_obj
  ) const
{
  using Teuchos::dyn_cast;
  typedef Teuchos::PrimitiveTypeTraits<Scalar> PTT;
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  ConstSubVectorView<Scalar> &sub_vec = dyn_cast<ReductTargetSubVectorT<Scalar> >(*reduct_obj).sub_vec;
  Scalar *sv_values = const_cast<Scalar*>(sub_vec.values()); 
  int value_data_off = 0;
  for( int k = 0; k < sub_vec.subDim(); ++k, value_data_off += num_prim_objs_per_scalar )
    PTT::loadPrimitiveObjs( num_prim_objs_per_scalar, value_data+value_data_off, &sv_values[k] );
}

template<class Scalar>
void ROpGetSubVector<Scalar>::get_op_type_num_entries(
  int*  num_values
  ,int* num_indexes
  ,int* num_chars
  ) const
{
  TEST_FOR_EXCEPT( !num_values || !num_indexes || !num_chars ); 
  *num_values  = 0;
  *num_indexes = 2; // l, u
  *num_chars   = 0;
}

template<class Scalar>
void ROpGetSubVector<Scalar>::extract_op_state(
  int                             num_values
  ,primitive_value_type           value_data[]
  ,int                            num_indexes
  ,Teuchos_Index                index_data[]
  ,int                            num_chars
  ,::RTOpPack::char_type                 char_data[]
  ) const
{
  TEST_FOR_EXCEPT( num_values!=0 || num_indexes!=2 || num_chars!=0 ); 
  TEST_FOR_EXCEPT( !index_data ); 
  index_data[0] = l_;
  index_data[1] = u_;
}

template<class Scalar>
void ROpGetSubVector<Scalar>::load_op_state(
  int                           num_values
  ,const primitive_value_type   value_data[]
  ,int                          num_indexes
  ,const Teuchos_Index        index_data[]
  ,int                          num_chars
  ,const ::RTOpPack::char_type         char_data[]
  )
{
  TEST_FOR_EXCEPT( num_values!=0 || num_indexes!=2 || num_chars!=0 ); 
  TEST_FOR_EXCEPT( !index_data ); 
  l_ = index_data[0];
  u_ = index_data[1];
}

template<class Scalar>
bool ROpGetSubVector<Scalar>::coord_invariant() const
{
  return false;
}

template<class Scalar>
void ROpGetSubVector<Scalar>::apply_op(
  const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
  ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
  ,ReductTarget *reduct_obj
  ) const
{
  using Teuchos::dyn_cast;

  RTOP_APPLY_OP_1_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);

  if( u_ < globalOffset || globalOffset + subDim - 1 < l_ )
    return; // None of the sub-vector elements that we are looking for is in this vector chunk!
  
  index_type
    i_l = ( l_ <= globalOffset           ? 0          : l_ - globalOffset  ),
    i_u = ( u_ >= globalOffset+subDim-1  ? subDim-1   : u_ - globalOffset  );

  ConstSubVectorView<Scalar> &sub_vec_targ = dyn_cast<ReductTargetSubVectorT<Scalar> >(*reduct_obj).sub_vec;
  Scalar *svt_values = const_cast<Scalar*>(sub_vec_targ.values());

  for( index_type i = i_l; i <= i_u; ++i )
    svt_values[i+(globalOffset-l_)] = v0_val[i*v0_s];

}

} // namespace RTOpPack

#endif // RTOPPACK_ROP_GET_SUB_VECTOR_HPP
