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

#ifndef RTOPPACK_TOP_SET_SUB_VECTOR_DEF_HPP
#define RTOPPACK_TOP_SET_SUB_VECTOR_DEF_HPP


namespace RTOpPack {


template<class Scalar>
TOpSetSubVector<Scalar>::TOpSetSubVector()
  :RTOpT<Scalar>("TOpSetSubVector")
{}


template<class Scalar>
TOpSetSubVector<Scalar>::TOpSetSubVector( const SparseSubVectorT<Scalar> &sub_vec )
  :RTOpT<Scalar>("TOpSetSubVector")
{
  set_sub_vec(sub_vec);
}


template<class Scalar>
void TOpSetSubVector<Scalar>::set_sub_vec( const SparseSubVectorT<Scalar> &sub_vec )
{
  sub_vec_ = sub_vec;
}


// Overridden from RTOpT


template<class Scalar>
bool TOpSetSubVector<Scalar>::coord_invariant_impl() const
{
  return false;
}


template<class Scalar>
void TOpSetSubVector<Scalar>::apply_op_impl(
  const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
  const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
  const Ptr<ReductTarget> &reduct_obj
  ) const
{

  typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
  typedef typename Teuchos::ArrayRCP<const Teuchos_Index>::iterator const_indices_iter_t;
  
  validate_apply_op( *this, 0, 1, false, sub_vecs, targ_sub_vecs, reduct_obj.getConst() );

  // Get the dense subvector data that we will set
  const SubVectorView<Scalar> &z = targ_sub_vecs[0];
  const index_type z_global_offset = z.globalOffset();
  const index_type z_sub_dim = z.subDim();
  iter_t z_val = z.values().begin();
  const ptrdiff_t z_val_s = z.stride();

  // Get the sparse subvector data
  const index_type v_global_offset = sub_vec_.globalOffset();
  const index_type v_sub_dim = sub_vec_.subDim();
  const index_type v_sub_nz = sub_vec_.subNz();
  const_iter_t v_val = sub_vec_.values().begin();
  const ptrdiff_t v_val_s = sub_vec_.valuesStride();
  const bool has_v_ind = !is_null(sub_vec_.indices());
  const_indices_iter_t v_ind = sub_vec_.indices().begin();
  const ptrdiff_t v_ind_s = sub_vec_.indicesStride();
  const ptrdiff_t v_l_off = sub_vec_.localOffset();
  //const bool v_sorted = sub_vec_.isSorted();
 
  //
  // Set part of the sub-vector v for this chunk z.
  //

  if( v_global_offset + v_sub_dim < z_global_offset + 1
    || z_global_offset + z_sub_dim < v_global_offset + 1 )
  {
     // The sub-vector that we are setting does not overlap with this vector
     // chunk!
    return;
  }

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
  if (has_v_ind) {
    // Sparse elements
    // Set the overlapping elements to zero first.
    if( v_global_offset >= z_global_offset )
      z_val += (v_global_offset - z_global_offset) * z_val_s;
    for( index_type k = 0; k < num_overlap; ++k, z_val += z_val_s )
      *z_val = 0.0;
    // Now set the sparse entries
    z_val = targ_sub_vecs[0].values().begin();
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


#endif // RTOPPACK_TOP_SET_SUB_VECTOR_DEF_HPP
