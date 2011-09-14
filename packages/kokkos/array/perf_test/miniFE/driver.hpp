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

#include <cstddef>
#include <SparseOutput.hpp>
#include <BoxMeshFixture.hpp>
#include <CRSGraph.hpp>

#define  PRINT_SAMPLE_OF_SOLUTION  0

namespace Test {

template<class DeviceType >
void run_kernel(int, int, int, double*);

template<>
void run_kernel<KOKKOS_MACRO_DEVICE>(int x, int y, int z, double* times) 
{
  typedef KOKKOS_MACRO_DEVICE    device_type;
  typedef device_type::size_type index_type ;
  typedef double                 Scalar ;

  typedef Kokkos::MDArrayView<Scalar,     device_type>  scalar_array_d;
  typedef Kokkos::MDArrayView<index_type, device_type>  index_array_d;    

  typedef Kokkos::MultiVectorView<Scalar,     device_type>  scalar_vector_d;
  typedef Kokkos::MultiVectorView<index_type, device_type>  index_vector_d;

  typedef scalar_array_d::HostView  scalar_array_h ;
  typedef index_array_d ::HostView  index_array_h ;

  typedef scalar_vector_d::HostView  scalar_vector_h ;
  typedef index_vector_d ::HostView  index_vector_h ;

  // Problem coefficients

  const Scalar elem_coeff_K = 2 ;
  const Scalar elem_load_Q  = 1 ;

  //  Host Data Structures

  index_vector_h  A_row_h , A_col_h ;

  //  Device Data Structures

  scalar_array_d  elem_stiffness, elem_load ;

  index_vector_d A_row_d , A_col_d, dirichlet_flag_d ;
  scalar_vector_d A, b, X , dirichlet_value_d ;

  //------------------------------
  // Generate mesh and corresponding sparse matrix graph

  Kokkos::Impl::Timer wall_clock ;

  const BoxMeshFixture< Scalar , device_type > mesh( x , y , z );

  mesh.init_dirichlet_z( dirichlet_flag_d , dirichlet_value_d );

  init_crsgraph( mesh.h_mesh.node_elem_offset ,
                 mesh.h_mesh.node_elem_ids ,
                 mesh.h_mesh.elem_node_ids ,
                 A_row_h , A_col_h );

  // Copy sparse matrix graph to device

  A_row_d = Kokkos::create_labeled_multivector< index_array_d >("A_row_d",A_row_h.length());
  A_col_d = Kokkos::create_labeled_multivector< index_array_d >("A_col_d",A_col_h.length());

  Kokkos::deep_copy(A_row_d, A_row_h);
  Kokkos::deep_copy(A_col_d, A_col_h);

  device_type::wait_functor_completion();

  times[0] = wall_clock.seconds(); // Mesh and graph allocation and population.

  //------------------------------
  // Allocate device memory for linear system and element contributions.

  A = Kokkos::create_labeled_multivector< scalar_vector_d > ("A",A_col_h.length());  
  b = Kokkos::create_labeled_multivector< scalar_vector_d > ("b",mesh.elem_count, 8);  
  X = Kokkos::create_labeled_multivector< scalar_vector_d > ("X",mesh.node_count);

  elem_stiffness =  Kokkos::create_mdarray< scalar_array_d > (mesh.elem_count, 8, 8);
  elem_load      =  Kokkos::create_mdarray< scalar_array_d > (mesh.elem_count, 8);

  wall_clock.reset();

  Kokkos::parallel_for( mesh.elem_count,
    assembleFE<Scalar, device_type>( mesh.d_mesh.elem_node_ids ,
                                     mesh.d_mesh.node_coords ,
                                     elem_stiffness, elem_load ,
                                     elem_coeff_K , elem_load_Q ) );

  Kokkos::parallel_for( mesh.node_count,
    CRSMatrixGatherFill<Scalar, device_type>( A, b, A_row_d, A_col_d,
                                              mesh.d_mesh.node_elem_offset ,
                                              mesh.d_mesh.node_elem_ids,
                                              mesh.d_mesh.elem_node_ids,
                                              elem_stiffness,
                                              elem_load) );

  Kokkos::parallel_for(mesh.node_count ,
    Dirichlet<Scalar , device_type>(A, A_row_d ,A_col_d, b,
                                    dirichlet_flag_d , dirichlet_value_d ) );

  device_type::wait_functor_completion();

  times[1] = wall_clock.seconds(); // Matrix computation and assembly

//  printSparse< scalar_vector_d , index_vector_d>("A.txt",A,A_row_d,A_col_d);

  //------------------------------
  // Solve linear sytem

  times[2] = CG_Solve<Scalar, device_type>::run(A , A_row_d, A_col_d , b , X ,times );

#if  PRINT_SAMPLE_OF_SOLUTION

  scalar_vector_h X_h Kokkos::mirror_create( X );

  Kokkos::mirror_update( X_h , X );

  for ( int i = 0 ; i < (int) mesh.node_count_z ; ++i ) {
    const int ix = mesh.node_count_x - 1 ;
    const int iy = mesh.node_count_y - 1 ;
    const Scalar val00 = X_h( mesh.node_id( 0 , 0 , i ) );
    const Scalar valX0 = X_h( mesh.node_id( ix , 0 , i ) );
    const Scalar valXY = X_h( mesh.node_id( ix , iy , i ) );
    const Scalar val0Y = X_h( mesh.node_id( 0 , iy , i ) );
    std::cout << "corners_00_X0_XY_Y0("<<i<<") = "
              << val00 << " " << valX0 << " "
              << valXY << " " << val0Y << std::endl ;
  }

#endif

//  printGLUT<Scalar , scalar_vector_d , scalar_array_h , index_array_h>
//      ("X.txt", X , elem_coords_h , elem_nodeIDs_h,x,y,z);
}

} // namespace Test

