/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

//----------------------------------------------------------------------------

namespace HybridFEM {

//----------------------------------------------------------------------------

template< typename ScalarType ,
          unsigned ElemNode ,
          typename CoordScalarType ,
          class elem_matrices_type ,
          class elem_vectors_type >
struct GatherFill< 
  KokkosArray::CrsMatrix< ScalarType , KOKKOSARRAY_MACRO_DEVICE > ,
  FEMesh< CoordScalarType , ElemNode , KOKKOSARRAY_MACRO_DEVICE > ,
  elem_matrices_type , elem_vectors_type >
{
  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  static const size_type ElemNodeCount = ElemNode ;

  typedef KokkosArray::CrsMatrix< ScalarType , device_type >    matrix_type ;
  typedef typename matrix_type::coefficients_type   coefficients_type ;
  typedef KokkosArray::View< ScalarType[] , device_type >  vector_type ;
  typedef KokkosArray::View< size_type[][ElemNodeCount][ElemNodeCount] , device_type >       elem_graph_type ;

  typedef FEMesh< CoordScalarType , ElemNodeCount , device_type > mesh_type ;
  typedef typename mesh_type::node_elem_ids_type node_elem_ids_type ;

private:

  node_elem_ids_type  node_elem_ids ;
  elem_graph_type     elem_graph ;
  elem_matrices_type  elem_matrices ;
  elem_vectors_type   elem_vectors ;
  coefficients_type   system_coeff ;
  vector_type         system_rhs ;

public:

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( size_type irow ) const
  {
    const size_type node_elem_begin = node_elem_ids.row_map[irow];
    const size_type node_elem_end   = node_elem_ids.row_map[irow+1];

    //  for each element that a node belongs to 

    for ( size_type i = node_elem_begin ; i < node_elem_end ; i++ ) {

      const size_type elem_id   = node_elem_ids.entries( i, 0);
      const size_type row_index = node_elem_ids.entries( i, 1);

      system_rhs(irow) += elem_vectors(elem_id, row_index);

      //  for each node in a particular related element  
      //  gather the contents of the element stiffness
      //  matrix that belong in irow

      for ( size_type j = 0 ; j < ElemNodeCount ; ++j ){
        const size_type A_index = elem_graph( elem_id , row_index , j );

        system_coeff( A_index ) += elem_matrices( elem_id, row_index, j );
      }
    }
  }


  static void apply( const matrix_type & matrix ,
                     const vector_type & rhs ,
                     const mesh_type   & mesh ,
                     const elem_graph_type    & elem_graph ,
                     const elem_matrices_type & elem_matrices ,
                     const elem_vectors_type  & elem_vectors )
  {
    const size_t row_count = matrix.graph.row_map.dimension(0) - 1 ;
    GatherFill op ;
    op.node_elem_ids = mesh.node_elem_ids ;
    op.elem_graph    = elem_graph ;
    op.elem_matrices = elem_matrices ;
    op.elem_vectors  = elem_vectors ;
    op.system_coeff  = matrix.coefficients ;
    op.system_rhs    = rhs ;

    parallel_for( row_count , op );
  }
};

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


