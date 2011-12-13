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
#include "GMRES_Solve.hpp"

#define  PRINT_SAMPLE_OF_SOLUTION  0

namespace Test {

template< typename Scalar , class Device >
struct MiniImplTherm ;

struct PerformanceData {
  double mesh_time ;
  double elem_time ;
  size_t elem_flop ;
  double fill_time ;
  size_t fill_flop ;
  size_t solve_flops ;
  double gmres_solve_time ;

  PerformanceData()
  : mesh_time(0)
  , elem_time(0)
  , elem_flop(0)
  , fill_time(0)
  , fill_flop(0)
  , solve_flops(0)
  , gmres_solve_time(0)
  {}

};

template< typename Scalar >
struct MiniImplTherm< Scalar , KOKKOS_MACRO_DEVICE > {

static void run(int x, int y, int z, PerformanceData & perf )
{
  typedef KOKKOS_MACRO_DEVICE    device_type;
  typedef device_type::size_type index_type ;

  typedef Kokkos::MDArrayView<Scalar,     device_type>  scalar_array_d;
  typedef Kokkos::MDArrayView<index_type, device_type>  index_array_d;    

  typedef Kokkos::MultiVectorView<Scalar,     device_type>  scalar_vector_d;
  typedef Kokkos::MultiVectorView<index_type, device_type>  index_vector_d;

  typedef typename scalar_array_d::HostView  scalar_array_h ;
  typedef typename index_array_d ::HostView  index_array_h ;

  typedef typename scalar_vector_d::HostView  scalar_vector_h ;
  typedef typename index_vector_d ::HostView  index_vector_h ;

  // Problem coefficients

  const Scalar elem_coeff_K = 2 ;
  const Scalar elem_load_Q  = 1 ;

  //  Host Data Structures

  index_array_h   elem_graph_col_h ;
  index_vector_h  A_row_h , A_col_h ;

  //  Device Data Structures

  scalar_array_d  elem_stiffness, elem_load ;

  index_array_d  elem_graph_col_d ;
  index_vector_d A_row_d , A_col_d, dirichlet_flag_d ;
  scalar_vector_d A, b, X , dirichlet_value_d ;

  //------------------------------
  // Generate mesh and corresponding sparse matrix graph

  Kokkos::Impl::Timer wall_clock ;

  const BoxMeshFixture< double , device_type > mesh( x , y , z );

  mesh.init_dirichlet_z( dirichlet_flag_d , dirichlet_value_d );

  elem_graph_col_d =
    Kokkos::create_mdarray< index_array_d >(
      mesh.h_mesh.elem_node_ids.dimension(0) ,
      mesh.h_mesh.elem_node_ids.dimension(1) ,
      mesh.h_mesh.elem_node_ids.dimension(1) );

  elem_graph_col_h = mirror_create( elem_graph_col_d );

  init_crsgraph( mesh.h_mesh.node_elem_offset ,
                 mesh.h_mesh.node_elem_ids ,
                 mesh.h_mesh.elem_node_ids ,
                 elem_graph_col_h ,
                 A_row_h , A_col_h );

  mirror_update( elem_graph_col_d , elem_graph_col_h );

  // Copy sparse matrix graph to device

  A_row_d = Kokkos::create_labeled_multivector< index_array_d >("A_row_d",A_row_h.length());
  A_col_d = Kokkos::create_labeled_multivector< index_array_d >("A_col_d",A_col_h.length());

  Kokkos::deep_copy(A_row_d, A_row_h);
  Kokkos::deep_copy(A_col_d, A_col_h);

  device_type::wait_functor_completion();

  perf.mesh_time = wall_clock.seconds(); // Mesh and graph allocation and population.

  //------------------------------
  // Allocate device memory for linear system and element contributions.

  A = Kokkos::create_labeled_multivector< scalar_vector_d > ("A",A_col_h.length());  
  b = Kokkos::create_labeled_multivector< scalar_vector_d > ("b",mesh.elem_count, 8);  
  X = Kokkos::create_labeled_multivector< scalar_vector_d > ("X",mesh.node_count);

  elem_stiffness =  Kokkos::create_mdarray< scalar_array_d > (mesh.elem_count, 8, 8);
  elem_load      =  Kokkos::create_mdarray< scalar_array_d > (mesh.elem_count, 8);

  wall_clock.reset();

  typedef ElementComp< Scalar , double , device_type > ElementFunctor ;
  typedef CRSMatrixGatherFill<Scalar, device_type> GatherFillFunctor ;

  Kokkos::parallel_for( mesh.elem_count,
    ElementFunctor( mesh.d_mesh.elem_node_ids ,
                    mesh.d_mesh.node_coords ,
                    elem_stiffness, elem_load ,
                    elem_coeff_K , elem_load_Q ) );

  device_type::wait_functor_completion();

  // Element computation time and flops
  perf.elem_time = wall_clock.seconds();
  perf.elem_flop = (size_t) ElementFunctor::FLOPS_operator *
                   (size_t) mesh.elem_count ;

  wall_clock.reset();

  Kokkos::parallel_for( mesh.node_count,
    GatherFillFunctor( A, b, A_row_d, A_col_d,
                       mesh.d_mesh.node_elem_offset ,
                       mesh.d_mesh.node_elem_ids,
                       elem_graph_col_d,
                       elem_stiffness,
                       elem_load ) );

  device_type::wait_functor_completion();

  // Matrix gather-fill time and flops
  perf.fill_time = wall_clock.seconds();
  perf.fill_flop = (size_t) elem_stiffness.size() + (size_t) elem_load.size();

  //------------------------------

  Kokkos::parallel_for(mesh.node_count ,
    Dirichlet<Scalar , device_type>(A, A_row_d ,A_col_d, b,
                                    dirichlet_flag_d , dirichlet_value_d ) );


//  printSparse< scalar_vector_d , index_vector_d>("A.txt",A,A_row_d,A_col_d);

  //------------------------------
  // Solve linear sytem

  const size_t num_iters = 50 ;
  GMRES_Solve<Scalar, device_type>::run(A , A_row_d, A_col_d , b , X, num_iters , perf.solve_flops, perf.gmres_solve_time );

#if  PRINT_SAMPLE_OF_SOLUTION

  scalar_vector_h X_h = Kokkos::mirror_create( X );

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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

static void driver( const char * label , int beg , int end , int runs )
{

  PerformanceData ** perf = new PerformanceData*[end-beg];
  PerformanceData * mins = new PerformanceData[end-beg];
  PerformanceData * maxes = new PerformanceData[end-beg];
  PerformanceData * medians = new PerformanceData[end-beg];
  int *num_elements = new int[end-beg];
   
  for(int i = beg ; i < end; ++i )
  {
    const int ix = (int) cbrt( (double) ( 1 << i ) );
    const int iy = ix + 1 ;
    const int iz = iy + 1 ;
    const int n  = ix * iy * iz ;
    num_elements[i-beg] = n;
   
    perf[i-beg] = new PerformanceData[runs];

    for(int j = 0; j < runs; j++){

      run(ix,iy,iz,perf[i-beg][j]);

    }
  }

  double *mesh_times = new double[runs];
  double *elem_times = new double[runs];
  double *fill_times = new double[runs];
  double *gmres_solve_times = new double[runs];
  for(int i = 0 ; i < end-beg; ++i){
 
    mins[i].elem_flop = perf[i][0].elem_flop;
    maxes[i].elem_flop = perf[i][0].elem_flop;
    medians[i].elem_flop = perf[i][0].elem_flop;
    
    mins[i].fill_flop = perf[i][0].fill_flop;
    maxes[i].fill_flop = perf[i][0].fill_flop;
    medians[i].fill_flop = perf[i][0].fill_flop;
    
    mins[i].solve_flops = perf[i][0].solve_flops;
    maxes[i].solve_flops = perf[i][0].solve_flops;
    medians[i].solve_flops = perf[i][0].solve_flops;

    for(int j=0;j<runs; ++j){
      mesh_times[j] = perf[i][j].mesh_time;
      elem_times[j] = perf[i][j].elem_time;
      fill_times[j] = perf[i][j].fill_time;
      gmres_solve_times[j] = perf[i][j].gmres_solve_time;
    }
    std::sort(mesh_times, mesh_times+runs);
    std::sort(elem_times, elem_times+runs);
    std::sort(fill_times, fill_times+runs);
    std::sort(gmres_solve_times, gmres_solve_times+runs);

    mins[i].mesh_time = mesh_times[0];
    maxes[i].mesh_time = mesh_times[runs-1];
    medians[i].mesh_time = mesh_times[runs/2];

    mins[i].elem_time = elem_times[0];
    maxes[i].elem_time = elem_times[runs-1];
    medians[i].elem_time = elem_times[runs/2];

    mins[i].fill_time = fill_times[0];
    maxes[i].fill_time = fill_times[runs-1];
    medians[i].fill_time = fill_times[runs/2];

    mins[i].gmres_solve_time = gmres_solve_times[0];
    maxes[i].gmres_solve_time = gmres_solve_times[runs-1];
    medians[i].gmres_solve_time = gmres_solve_times[runs/2];
  }
  delete [] mesh_times;
  delete [] elem_times;
  delete [] fill_times;
  delete [] gmres_solve_times;



  std::cout << std::endl ;
  std::cout << "Minimum" << std::endl;
  std::cout << "\"MiniImplTherm with Kokkos " << label << "\"" << std::endl;
  std::cout << "\"Size\" ,     \"Setup\" ,    \"Element\" ,  \"Element\" , \"Fill\" ,   \"Fill\" ,      \"Solve\" ,    \"Solve\" ,     \"Solve\" "  << std::endl
            << "\"elements\" , \"millisec\" , \"millisec\" , \"flops\" , \"millisec\" , \"flops\" ,     \"flops\" ,     \"time\" , \"Mflops/sec\" " << std::endl ; 

  for(int i=0; i<end-beg; ++i){

   std::cout << std::setw(8) << num_elements[i] << " , "
             << std::setw(10) << mins[i].mesh_time * 1000 << " , "
             << std::setw(10) << mins[i].elem_time * 1000 << " , "
             << std::setw(10) << mins[i].elem_flop << " , "
             << std::setw(10) << mins[i].fill_time * 1000 << " , "
             << std::setw(10) << mins[i].fill_flop << " , "
             << std::setw(10) << mins[i].solve_flops << " , "
             << std::setw(10) << mins[i].gmres_solve_time << " , "
             << std::setw(10) << 
              (double)mins[i].solve_flops / ( mins[i].gmres_solve_time * 1.0e6)
             << std::endl ;
  }

  std::cout << std::endl ;
  std::cout << "Maximum" << std::endl;
  std::cout << "\"MiniImplTherm with Kokkos " << label << "\"" << std::endl;
  std::cout << "\"Size\" ,     \"Setup\" ,    \"Element\" ,  \"Element\" , \"Fill\" ,   \"Fill\" ,      \"Solve\" ,    \"Solve\" ,     \"Solve\" "  << std::endl
            << "\"elements\" , \"millisec\" , \"millisec\" , \"flops\" , \"millisec\" , \"flops\" ,     \"flops\" ,     \"time\" , \"Mflops/sec\" " << std::endl ; 

  for(int i=0; i<end-beg; ++i){

   std::cout << std::setw(8) << num_elements[i] << " , "
             << std::setw(10) << maxes[i].mesh_time * 1000 << " , "
             << std::setw(10) << maxes[i].elem_time * 1000 << " , "
             << std::setw(10) << maxes[i].elem_flop << " , "
             << std::setw(10) << maxes[i].fill_time * 1000 << " , "
             << std::setw(10) << maxes[i].fill_flop << " , "
             << std::setw(10) << maxes[i].solve_flops << " , "
             << std::setw(10) << maxes[i].gmres_solve_time << " , "
             << std::setw(10) << 
              (double)maxes[i].solve_flops / ( maxes[i].gmres_solve_time * 1.0e6)
             << std::endl ;
  }

  std::cout << std::endl ;
  std::cout << "Median" << std::endl;
  std::cout << "\"MiniImplTherm with Kokkos " << label << "\"" << std::endl;
  std::cout << "\"Size\" ,     \"Setup\" ,    \"Element\" ,  \"Element\" , \"Fill\" ,   \"Fill\" ,      \"Solve\" ,    \"Solve\" ,     \"Solve\" "  << std::endl
            << "\"elements\" , \"millisec\" , \"millisec\" , \"flops\" , \"millisec\" , \"flops\" ,     \"flops\" ,     \"time\" , \"Mflops/sec\" " << std::endl ; 

  for(int i=0; i<end-beg; ++i){

   std::cout << std::setw(8) << num_elements[i] << " , "
             << std::setw(10) << medians[i].mesh_time * 1000 << " , "
             << std::setw(10) << medians[i].elem_time * 1000 << " , "
             << std::setw(10) << medians[i].elem_flop << " , "
             << std::setw(10) << medians[i].fill_time * 1000 << " , "
             << std::setw(10) << medians[i].fill_flop << " , "
             << std::setw(10) << medians[i].solve_flops << " , "
             << std::setw(10) << medians[i].gmres_solve_time << " , "
             << std::setw(10) << 
              (double)medians[i].solve_flops / ( medians[i].gmres_solve_time * 1.0e6)
             << std::endl ;
  }

  for(int i = beg ; i < end; ++i )
  {
    delete[] perf[i];
  }
  delete[] perf;
  delete[] num_elements;
  delete[] mins;
  delete[] maxes;
  delete[] medians;
}

};

} // namespace Test

