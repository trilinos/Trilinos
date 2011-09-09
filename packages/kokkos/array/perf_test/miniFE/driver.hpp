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
  typedef index_array_d   ::HostView  index_array_h ;

  typedef scalar_vector_d::HostView  scalar_vector_h ;
  typedef index_vector_d   ::HostView  index_vector_h ;

  //  Host Data Structures

  index_vector_h  A_row_h , A_col_h;

  //  Device Data Structures

  scalar_array_d  elem_stiffness, elem_load;

  index_vector_d A_row_d , A_col_d;
  scalar_vector_d A, b, X;
  
  timeval start,stop,result;
  double time = 0.0;
  double total_time = 0.0;

  gettimeofday(&start, NULL);

  const BoxMeshFixture< Scalar , device_type > mesh( x , y , z );

  init_crsgraph( mesh.h_mesh.node_elem_offset ,
                 mesh.h_mesh.node_elem_ids ,
                 mesh.h_mesh.elem_node_ids ,
                 A_row_h , A_col_h );

  gettimeofday(&stop, NULL);
  timersub(&stop, &start, &result);
  time = (result.tv_sec + result.tv_usec/1000000.0);

  total_time += time;  
  times[0] = time;

//  copy host data to device

  A_row_d = Kokkos::create_labeled_multivector< index_array_d >("A_row_d",A_row_h.length());
  A_col_d = Kokkos::create_labeled_multivector< index_array_d >("A_col_d",A_col_h.length());

  elem_stiffness =  Kokkos::create_mdarray< scalar_array_d > (mesh.elem_count, 8, 8);
  elem_load      =  Kokkos::create_mdarray< scalar_array_d > (mesh.elem_count, 8);

  A = Kokkos::create_labeled_multivector< scalar_vector_d > ("A",A_col_h.length());  
  b = Kokkos::create_labeled_multivector< scalar_vector_d > ("b",mesh.elem_count, 8);  
  X = Kokkos::create_labeled_multivector< scalar_vector_d > ("X",mesh.node_count);

  gettimeofday(&start, NULL);

  Kokkos::deep_copy(A_row_d, A_row_h);
  Kokkos::deep_copy(A_col_d, A_col_h);

  gettimeofday(&stop, NULL);
  timersub(&stop, &start, &result);
  time = (result.tv_sec + result.tv_usec/1000000.0);

  total_time += time;  


  Kokkos::parallel_for( mesh.elem_count,
    assembleFE<Scalar, device_type>( mesh.d_mesh.elem_node_ids ,
                                     mesh.d_mesh.node_coords ,
                                     elem_stiffness, elem_load ), time);

  total_time += time;  
  times[1] = time;


  Kokkos::parallel_for(mesh.node_count,
    CRSMatrixGatherFill<Scalar, device_type>( A, b, A_row_d, A_col_d,
                                              mesh.d_mesh.node_elem_offset ,
                                              mesh.d_mesh.node_elem_ids,
                                              mesh.d_mesh.elem_node_ids,
                                              elem_stiffness,
                                              elem_load), time);

  total_time += time;  
  times[2] = time;

  Kokkos::parallel_for(mesh.node_count , Dirichlet<Scalar , device_type>(A, A_row_d ,A_col_d, b,x+1,y+1,z+1,1.0) , time);
  total_time += time;
  times[3] = time;

//  printSparse< scalar_vector_d , index_vector_d>("A.txt",A,A_row_d,A_col_d);

  time = CG_Solve<Scalar, device_type>::run(A , A_row_d, A_col_d , b , X ,times );
  total_time += time;
  
  times[6] = total_time;

//  printGLUT<Scalar , scalar_vector_d , scalar_array_h , index_array_h>
//      ("X.txt", X , elem_coords_h , elem_nodeIDs_h,x,y,z);
}

} //Test
