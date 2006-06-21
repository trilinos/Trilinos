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

#ifndef RTOPPACK_TOP_SET_SUB_VECTOR_HPP
#define RTOPPACK_TOP_SET_SUB_VECTOR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Advanced transformation operator that assigns elements from a
 * sparse explicit vector.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class TOpSetSubVector : public RTOpT<Scalar> {
public:

  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;

  /** \brief . */
  TOpSetSubVector( const SparseSubVectorT<Scalar> &sub_vec );

  /** \brief . */
  void set_sub_vec( const SparseSubVectorT<Scalar> &sub_vec );

  /** @name Overridden from RTOpT */
  //@{

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
    ,index_type                     index_data[]
    ,int                            num_chars
    ,char_type                      char_data[]
    ) const;
  /** \brief . */
  void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const index_type             index_data[]
    ,int                          num_chars
    ,const char_type              char_data[]
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

  // ////////////////////////////
  // Private types

  enum { num_sub_vec_members = 6 };

  // ////////////////////////////
  // Private data members

  SparseSubVectorT<Scalar>  sub_vec_;  // We do not own memory it its arrays!
  bool                      ownsMem_;

}; // class TOpSetSubVector

// ////////////////////////////////
// Template definitions

template<class Scalar>
TOpSetSubVector<Scalar>::TOpSetSubVector( const SparseSubVectorT<Scalar> &sub_vec )
  :RTOpT<Scalar>("TOpSetSubVector"), sub_vec_(sub_vec),ownsMem_(false)
{}

template<class Scalar>
void TOpSetSubVector<Scalar>::set_sub_vec( const SparseSubVectorT<Scalar> &sub_vec )
{
  if( ownsMem_ ) {
    if(sub_vec_.values())  delete [] const_cast<Scalar*>(sub_vec_.values());
    if(sub_vec_.indices()) delete [] const_cast<index_type*>(sub_vec_.indices());
  }
  sub_vec_ = sub_vec;
  ownsMem_ = false;
}

// Overridden from RTOpT

template<class Scalar>
void TOpSetSubVector<Scalar>::get_op_type_num_entries(
  int*  num_values
  ,int* num_indexes
  ,int* num_chars
  ) const
{
  typedef Teuchos::PrimitiveTypeTraits<Scalar> PTT;
  TEST_FOR_EXCEPT( !num_values || !num_indexes || !num_chars ); 
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  *num_values = num_prim_objs_per_scalar*sub_vec_.subNz();  // values[]
  *num_indexes =
    num_sub_vec_members // globalOffset,subDim,subNz,localOffset,isSorted,ownsMem
    + (sub_vec_.indices() ? sub_vec_.subNz() : 0 ); // indices[]
  *num_chars   = 0;
}

template<class Scalar>
void TOpSetSubVector<Scalar>::extract_op_state(
  int                             num_values
  ,primitive_value_type           value_data[]
  ,int                            num_indexes
  ,index_type                     index_data[]
  ,int                            num_chars
  ,char_type                      char_data[]
  ) const
{
  typedef Teuchos::PrimitiveTypeTraits<Scalar> PTT;
  register index_type k, j;
#ifdef RTOp_DEBUG
  TEST_FOR_EXCEPT(!(num_values==sub_vec_.subNz()));
  TEST_FOR_EXCEPT(!(num_indexes==num_sub_vec_members+(sub_vec_.indices()?sub_vec_.subNz():0)));
  TEST_FOR_EXCEPT(!num_chars==0);
#endif
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  index_type vd_off = 0;
  for( k = 0; k < sub_vec_.subNz(); ++k, vd_off += num_prim_objs_per_scalar ) {
    PTT::extractPrimitiveObjs(
      *(sub_vec_.values() + k*sub_vec_.valuesStride())
      ,num_prim_objs_per_scalar, value_data+vd_off
      );
  }
  index_data[k=0] = sub_vec_.globalOffset();
  index_data[++k] = sub_vec_.subDim();
  index_data[++k] = sub_vec_.subNz();
  index_data[++k] = sub_vec_.localOffset();
  index_data[++k] = sub_vec_.isSorted();
  index_data[++k] = ownsMem_ ? 1 : 0;
  if( sub_vec_.indices() ) {
    for( j = 0; j < sub_vec_.subNz(); ++j )
      index_data[++k] = *(sub_vec_.indices() + j*sub_vec_.indicesStride());
  }
}

template<class Scalar>
void TOpSetSubVector<Scalar>::load_op_state(
  int                           num_values
  ,const primitive_value_type   value_data[]
  ,int                          num_indexes
  ,const index_type             index_data[]
  ,int                          num_chars
  ,const char_type              char_data[]
  )
{
  typedef Teuchos::PrimitiveTypeTraits<Scalar> PTT;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( num_chars!=0 );
  TEST_FOR_EXCEPT( !value_data );
  TEST_FOR_EXCEPT( !index_data );
#endif
  const int num_prim_objs_per_scalar = PTT::numPrimitiveObjs();
  index_type k, j;
  index_type Nz            = num_values / num_prim_objs_per_scalar;
  Scalar     *scalars      = const_cast<Scalar*>(sub_vec_.values());
  ptrdiff_t  scalarsStride = sub_vec_.valuesStride();
  index_type *indices      = const_cast<index_type*>(sub_vec_.indices());
  ptrdiff_t  indicesStride = sub_vec_.indicesStride();
  // Reallocate storage if we have to
  if( Nz != sub_vec_.subNz() || (indices==NULL && num_indexes > num_sub_vec_members ) ) {
    // The current sub_vec_ does not have storage setup to hold the incoming subvector.
    // Delete current storage if owned.
    if( ownsMem_ ) {
      if(scalars) delete [] scalars;
      if(indices) delete [] indices;
    }
    // We need to reallocate storage for scalars[] and perhaps indices[]
    scalars = new Scalar[Nz];
    if( num_indexes > num_sub_vec_members )
      indices = new index_type[Nz];
    ownsMem_ = true;
    scalarsStride = 1;
    indicesStride = 1;
  }
  else {
    // The storage in sub_vec_ is already correct to hold the incoming subvector
  }
  // Set the internal sub_vec
  index_type v_off = 0;
  for( k = 0; k < Nz; ++k, v_off += num_prim_objs_per_scalar )
    PTT::loadPrimitiveObjs( num_prim_objs_per_scalar, value_data+v_off, scalars+k*scalarsStride );
  const index_type globalOffset  = index_data[k=0];
  const index_type subDim        = index_data[++k];
  const index_type subNz         = index_data[++k];
  const ptrdiff_t  localOffset   = index_data[++k];
  const int        isSorted      = index_data[++k];
  ownsMem_                       = ( index_data[++k]==1 ? true : false );
  if( num_indexes > num_sub_vec_members ) {
    for( j = 0; j < num_values; ++j )
      *(indices+j*indicesStride) = index_data[++k];
  }
  TEST_FOR_EXCEPT( subNz != Nz );
  sub_vec_.initialize(globalOffset,subDim,subNz,scalars,scalarsStride,indices,indicesStride,localOffset,isSorted);
}

template<class Scalar>
bool TOpSetSubVector<Scalar>::coord_invariant() const
{
  return false;
}

template<class Scalar>
void TOpSetSubVector<Scalar>::apply_op(
  const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
  ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
  ,ReductTarget *reduct_obj
  ) const
{
  // Get the local vector chunk data
  RTOP_APPLY_OP_0_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
  const ptrdiff_t  z_global_offset = targ_sub_vecs[0].globalOffset();
  const index_type z_sub_dim       = subDim; // From macro above
  Scalar           *z_val          = z0_val; // From macro above
  const ptrdiff_t  z_val_s         = z0_s;   // From macro above
  // Get the sparse subvector data
  const index_type v_global_offset  = sub_vec_.globalOffset();
  const index_type v_sub_dim        = sub_vec_.subDim();
  const index_type v_sub_nz         = sub_vec_.subNz();
  const Scalar     *v_val           = sub_vec_.values();
  const ptrdiff_t  v_val_s          = sub_vec_.valuesStride();
  const index_type *v_ind           = sub_vec_.indices();
  const ptrdiff_t  v_ind_s          = sub_vec_.indicesStride();
  const ptrdiff_t  v_l_off          = sub_vec_.localOffset();
  //const bool       v_sorted         = sub_vec_.isSorted();
  
  //
  // Set part of the sub-vector v for this chunk z.
  //

  if( v_global_offset + v_sub_dim < z_global_offset + 1
    || z_global_offset + z_sub_dim < v_global_offset + 1 )
    return; // The sub-vector that we are setting does not overlap with this vector chunk!

  if( v_sub_nz == 0 )
    return; // The sparse sub-vector we are reading from is empty?

  // Get the number of elements that overlap
  index_type num_overlap;
  if( v_global_offset <= z_global_offset ) {
    if( v_global_offset + v_sub_dim >= z_global_offset + z_sub_dim )
      num_overlap = z_sub_dim;
    else
      num_overlap = (v_global_offset + v_sub_dim) - z_global_offset;
  }
  else {
    if( z_global_offset + z_sub_dim >= v_global_offset + v_sub_dim )
      num_overlap = v_sub_dim;
    else
      num_overlap = (z_global_offset + z_sub_dim) - v_global_offset;
  }
  
  // Set the part of the sub-vector that overlaps
  if( v_ind != NULL ) {
    // Sparse elements
    // Set the overlapping elements to zero first.
    if( v_global_offset >= z_global_offset )
      z_val += (v_global_offset - z_global_offset) * z_val_s;
    for( index_type k = 0; k < num_overlap; ++k, z_val += z_val_s )
      *z_val = 0.0;
    // Now set the sparse entries
    z_val = targ_sub_vecs[0].values();
    for( index_type k = 0; k < v_sub_nz; ++k, v_val += v_val_s, v_ind += v_ind_s ) {
      const index_type i = v_global_offset + v_l_off + (*v_ind);
      if( z_global_offset < i && i <= z_global_offset + z_sub_dim )
        z_val[ z_val_s * (i - z_global_offset - 1) ] = *v_val;
    }
    // ToDo: Implement a faster version for v sorted and eliminate the
    // if statement in the above loop.
  }
  else {
    // Dense elements
    if( v_global_offset <= z_global_offset )
      v_val += (z_global_offset - v_global_offset) * v_val_s;
    else
      z_val += (v_global_offset - z_global_offset) * z_val_s;
    for( index_type k = 0; k < num_overlap; ++k, v_val += v_val_s, z_val += z_val_s )
      *z_val = *v_val;
  }
}

} // namespace RTOpPack

#endif // RTOPPACK_TOP_SET_SUB_VECTOR_HPP
