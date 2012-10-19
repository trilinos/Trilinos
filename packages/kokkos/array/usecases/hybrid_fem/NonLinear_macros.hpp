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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

namespace HybridFEM {
namespace NonLinear {

//----------------------------------------------------------------------------

template< typename ScalarCoordType , unsigned ElemNode , typename ScalarType >
struct DirichletSolution<
  FEMesh< ScalarCoordType , ElemNode , KOKKOSARRAY_MACRO_DEVICE > , ScalarType >
{
  typedef KOKKOSARRAY_MACRO_DEVICE  device_type;

  static const unsigned ElementNodeCount = ElemNode ;

  typedef KokkosArray::View< ScalarType[] , device_type >  vector_type ;

  typedef FEMesh< ScalarCoordType , ElementNodeCount , device_type > mesh_type ;

  typename mesh_type::node_coords_type node_coords ;

  vector_type     solution ;
  ScalarCoordType bc_lower_z ;
  ScalarCoordType bc_upper_z ;
  ScalarType      bc_lower_value ;
  ScalarType      bc_upper_value ;

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const unsigned inode ) const
  {
   
  // Apply dirichlet boundary condition on the Solution vector.
  // Define boundary node values to be either bc_lower_value or
  // bc_upper_value, depending on which boundary face they lie on.
  // Non-boundary terms will be left at their previous value.

    const ScalarCoordType z = node_coords(inode,2);
    const bool bc_lower = z <= bc_lower_z ;
    const bool bc_upper = bc_upper_z <= z ;

    if ( bc_lower || bc_upper ) {
      const ScalarType bc_value = bc_lower ? bc_lower_value
                                           : bc_upper_value ;

      solution(inode) = bc_value ; //  set the solution vector
    }
  }

  static void apply( const vector_type    & solution ,
                     const mesh_type      & mesh ,
                     const ScalarCoordType  bc_lower_z ,
                     const ScalarCoordType  bc_upper_z ,
                     const ScalarType       bc_lower_value ,
                     const ScalarType       bc_upper_value )
  {
    DirichletSolution op ;
    op.node_coords    = mesh.node_coords ;
    op.solution       = solution ;
    op.bc_lower_z     = bc_lower_z ;
    op.bc_upper_z     = bc_upper_z ;
    op.bc_lower_value = bc_lower_value ;
    op.bc_upper_value = bc_upper_value ;
    parallel_for( solution.dimension_0() , op );
  }
};

template< typename ScalarCoordType , unsigned ElemNode , typename ScalarType >
struct DirichletResidual<
  FEMesh< ScalarCoordType , ElemNode , KOKKOSARRAY_MACRO_DEVICE > , ScalarType >
{
  typedef KOKKOSARRAY_MACRO_DEVICE     device_type;
  typedef device_type::size_type  size_type ;

  static const unsigned ElementNodeCount = ElemNode ;

  typedef KokkosArray::CrsMatrix< ScalarType , device_type >    matrix_type ;
  typedef KokkosArray::View< ScalarType[] , device_type >  vector_type ;

  typedef FEMesh< ScalarCoordType , ElementNodeCount , device_type > mesh_type ;

  typename mesh_type::node_coords_type node_coords ;
  matrix_type     matrix ;
  vector_type     rhs ;
  ScalarCoordType bc_lower_z ;
  ScalarCoordType bc_upper_z ;

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const unsigned inode ) const
  {
    //  Apply a dirichlet boundary condition to 'irow'
    //  to maintain the symmetry of the original 
    //  global stiffness matrix, zero out the columns
    //  that correspond to boundary conditions, and
    //  adjust the load vector accordingly

    const size_type iBeg = matrix.graph.row_map[inode];
    const size_type iEnd = matrix.graph.row_map[inode+1];

    const ScalarCoordType z = node_coords(inode,2);
    const bool bc_lower = z <= bc_lower_z ;
    const bool bc_upper = bc_upper_z <= z ;

    if ( bc_lower || bc_upper ) {
      rhs(inode) = 0 ; //  set the residual vector

      //  zero each value on the row, and leave a one
      //  on the diagonal

      for( size_type i = iBeg ; i < iEnd ; i++) {
        matrix.coefficients(i) =
          (int) inode == matrix.graph.entries(i) ? 1 : 0 ;
      }
    }
    else {

      //  Find any columns that are boundary conditions.
      //  Clear them and adjust the load vector

      for( size_type i = iBeg ; i < iEnd ; i++ ) {
        const size_type cnode = matrix.graph.entries(i) ;

        const ScalarCoordType zc = node_coords(cnode,2);
        const bool c_bc_lower = zc <= bc_lower_z ;
        const bool c_bc_upper = bc_upper_z <= zc ;

        if ( c_bc_lower || c_bc_upper ) {

	   matrix.coefficients(i) = 0 ;
        }
      }
    }
  }


  static void apply( const matrix_type & linsys_matrix ,
                     const vector_type & linsys_rhs ,
                     const mesh_type   & mesh ,
                     const ScalarCoordType  bc_lower_z ,
                     const ScalarCoordType  bc_upper_z)
  {
    const size_t row_count = linsys_matrix.graph.row_map.dimension(0) - 1 ;

    DirichletResidual op ;
    op.node_coords    = mesh.node_coords ;
    op.matrix         = linsys_matrix ;
    op.rhs            = linsys_rhs ;
    op.bc_lower_z     = bc_lower_z ;
    op.bc_upper_z     = bc_upper_z ;
    parallel_for( row_count , op );
  }
};

//----------------------------------------------------------------------------

} /* namespace NonLinear */
} /* namespace HybridFEM */

