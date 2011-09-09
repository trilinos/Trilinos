/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

template<typename Scalar , class DeviceType>
struct CRSMatrixGatherFill;

template<typename Scalar>
struct CRSMatrixGatherFill<Scalar ,KOKKOS_MACRO_DEVICE>{
  
  typedef KOKKOS_MACRO_DEVICE     device_type;
  typedef device_type::size_type  index_type;

  typedef Kokkos::MultiVectorView<Scalar , device_type>     scalar_vector_d;
  typedef Kokkos::MultiVectorView<index_type , device_type> index_vector_d;

  typedef Kokkos::MDArrayView<Scalar,device_type>      scalar_array_d;
  typedef Kokkos::MDArrayView<index_type,device_type>  index_array_d;    
  
  scalar_vector_d  A ;
  scalar_vector_d  b ;
  index_vector_d   A_col_offset;
  index_vector_d   A_col_index;

  index_array_d  node_elemIDs;
  index_array_d  elem_nodeIDs;
  index_array_d  elems_per_node;
  
  scalar_array_d  element_stiffness;
  scalar_array_d  element_load;

  CRSMatrixGatherFill(
    const scalar_vector_d   & arg_A,
    const scalar_vector_d   & arg_b,
    const index_vector_d    & arg_A_col_offset,
    const index_vector_d    & arg_A_col_index,
    const index_array_d     & arg_elems_per_node,
    const index_array_d     & arg_node_elemIDs,
    const index_array_d     & arg_elem_nodeIDs,
    const scalar_array_d    & arg_element_stiffness,
    const scalar_array_d    & arg_element_load)
  : A(arg_A), 
    b(arg_b),
    A_col_offset(arg_A_col_offset),
    A_col_index(arg_A_col_index),
    node_elemIDs(arg_node_elemIDs),
    elem_nodeIDs(arg_elem_nodeIDs),
    elems_per_node(arg_elems_per_node),
    element_stiffness(arg_element_stiffness),
    element_load(arg_element_load)
  {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()(index_type irow) const {

    const index_type base_index = A_col_offset(irow);
    const index_type last_index = A_col_offset(irow + 1);

    const index_type node_elem_begin = elems_per_node( irow );
    const index_type node_elem_end   = elems_per_node( irow + 1 );

    //  for each element that a node belongs to

    for(index_type i = node_elem_begin ; i < node_elem_end ; i++){

      //  elems_per_node is a cumulative structure, so 
      //  elems_per_node(irow) should be the index where
      //  a particular row's elem_IDs begin

      const index_type nelem          = node_elemIDs( i, 0);
      const index_type elem_row_index = node_elemIDs( i, 1);

      b(irow) += element_load(nelem, elem_row_index);

      //  for each node in a particular related element  
      for(index_type j = 0; j < 8; j++){

        //  gather the contents of the element stiffness
        //  matrix that belong in irow

        index_type column_search = base_index;
  
        const index_type node_id = elem_nodeIDs(nelem, j);

        for ( index_type len = last_index - base_index ; 0 < len ; ) {
  
          const index_type half = len >> 1;
          const index_type middle = column_search + half ;

          if ( A_col_index(middle) < node_id ){
            column_search = middle + 1 ;
            len -= half + 1 ;
          }
          else {
            len = half ;
          }
        }

        A(column_search) += element_stiffness(nelem, elem_row_index, j);
      }
    }
  }

}; //CRSMatrixGatherFill


