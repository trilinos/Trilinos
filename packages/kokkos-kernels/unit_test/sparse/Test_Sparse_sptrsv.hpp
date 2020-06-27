/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <Kokkos_Concepts.hpp>
#include <string>
#include <stdexcept>

#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_SparseUtils.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include "KokkosSparse_sptrsv.hpp"
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
#include "KokkosSparse_sptrsv_supernode.hpp"
#endif

#include<gtest/gtest.h>


using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Impl;
using namespace KokkosKernels::Experimental;

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

namespace Test {

#if 0
template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_sptrsv_mtx() {

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsmat_t;
  typedef typename crsmat_t::StaticCrsGraphType graph_t;

  //typedef Kokkos::View< size_type*, device >  RowMapType;
  //typedef Kokkos::View< lno_t*, device >      EntriesType;
  typedef Kokkos::View< scalar_t*, device >     ValuesType;

  // Lower tri
  std::cout << "LowerTriTest Begin" << std::endl;
  {

//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/L-offshore-amd.mtx";
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/L-Transport-amd.mtx";
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/L-Fault_639amd.mtx";
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/L-thermal2-amd.mtx";
    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/L-dielFilterV2real-amd.mtx";
    std::cout << "Matrix file: " << mtx_filename << std::endl;
    crsmat_t triMtx = KokkosKernels::Impl::read_kokkos_crst_matrix<crsmat_t>(mtx_filename.c_str()); //in_matrix
    graph_t  lgraph  = triMtx.graph; // in_graph

    auto row_map = lgraph.row_map;
    auto entries = lgraph.entries;
    auto values  = triMtx.values;

    const size_type nrows = lgraph.numRows();
//    const size_type nnz   = triMtx.nnz();

    scalar_t ZERO = scalar_t(0);
    scalar_t ONE = scalar_t(1);

    typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_type, lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

    std::cout << "UnitTest nrows = " << nrows << std::endl;

    KernelHandle kh;
    bool is_lower_tri = true;
    std::cout << "Create handle" << std::endl;
    kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
    
    std::cout << "Prepare linear system" << std::endl;
    // Create known_lhs, generate rhs, then solve for lhs to compare to known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

//    typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
//    crsMat_t triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);

    std::cout << "SPMV" << std::endl;
    KokkosSparse::spmv( "N", ONE, triMtx, known_lhs, ZERO, rhs);

    std::cout << "TriSolve Symbolic" << std::endl;
    Kokkos::Timer timer;
    sptrsv_symbolic( &kh, row_map, entries );
    std::cout << "LTRI Symbolic Time: " << timer.seconds() << std::endl;

    std::cout << "TriSolve Solve" << std::endl;
    kh.get_sptrsv_handle()->print_algorithm();
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    std::cout << "LTRI Solve TEAMPOLICY! Time: " << timer.seconds() << std::endl;

    scalar_t sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
    if ( sum != lhs.extent(0) ) {
      std::cout << "Lower Tri Solve FAILURE" << std::endl;
    }
    else {
      std::cout << "Lower Tri Solve SUCCESS!" << std::endl;
      //std::cout << "Num-levels = " << kh->get_sptrsv_handle()->get_num_levels() << std::endl;
    }
    EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

    Kokkos::deep_copy(lhs, 0);
    kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHD_RP);
    kh.get_sptrsv_handle()->print_algorithm();
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    std::cout << "LTRI Solve SEQLVLSCHD_RP Time: " << timer.seconds() << std::endl;

    sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
    if ( sum != lhs.extent(0) ) {
      std::cout << "Lower Tri Solve FAILURE" << std::endl;
    }
    else {
      std::cout << "Lower Tri Solve SUCCESS!" << std::endl;
      //std::cout << "Num-levels = " << kh->get_sptrsv_handle()->get_num_levels() << std::endl;
    }
    EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

    Kokkos::deep_copy(lhs, 0);
    kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHED_TP2);
    kh.get_sptrsv_handle()->print_algorithm();
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    std::cout << "LTRI Solve SEQLVLSCHED_TP2 Time: " << timer.seconds() << std::endl;

    sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
    if ( sum != lhs.extent(0) ) {
      std::cout << "Lower Tri Solve FAILURE" << std::endl;
    }
    else {
      std::cout << "Lower Tri Solve SUCCESS!" << std::endl;
      //std::cout << "Num-levels = " << kh->get_sptrsv_handle()->get_num_levels() << std::endl;
    }
    EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );


    kh.destroy_sptrsv_handle();
  }
  // Upper tri
  std::cout << "UpperTriTest Begin" << std::endl;
  {
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/U-offshore-amd.mtx";
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/U-Transport-amd.mtx";
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/U-Fault_639amd.mtx";
//    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/U-thermal2-amd.mtx";
    std::string mtx_filename = "/ascldap/users/ndellin/TestCodes-GitlabEx/KokkosEcoCodes/KokkosKernels-DevTests/Matrices/U-dielFilterV2real-amd.mtx";
    std::cout << "Matrix file: " << mtx_filename << std::endl;
    crsmat_t triMtx = KokkosKernels::Impl::read_kokkos_crst_matrix<crsmat_t>(mtx_filename.c_str()); //in_matrix
    graph_t  lgraph  = triMtx.graph; // in_graph

    auto row_map = lgraph.row_map;
    auto entries = lgraph.entries;
    auto values  = triMtx.values;

    const size_type nrows = lgraph.numRows();
//    const size_type nnz   = triMtx.nnz();

    scalar_t ZERO = scalar_t(0);
    scalar_t ONE = scalar_t(1);

    typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_type, lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

    std::cout << "UnitTest nrows = " << nrows << std::endl;

    KernelHandle kh;
    bool is_lower_tri = false;
    std::cout << "Create handle" << std::endl;
    kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
    
    std::cout << "Prepare linear system" << std::endl;
    // Create known_lhs, generate rhs, then solve for lhs to compare to known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

//    typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
//    crsMat_t triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);
    std::cout << "SPMV" << std::endl;
    KokkosSparse::spmv( "N", ONE, triMtx, known_lhs, ZERO, rhs);

    std::cout << "TriSolve Symbolic" << std::endl;
    Kokkos::Timer timer;
    sptrsv_symbolic( &kh, row_map, entries );
    std::cout << "UTRI Symbolic Time: " << timer.seconds() << std::endl;

    std::cout << "TriSolve Solve" << std::endl;
    kh.get_sptrsv_handle()->print_algorithm();
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    std::cout << "UTRI Solve SEQLVLSCHD_TP1 Time: " << timer.seconds() << std::endl;

    scalar_t sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
    if ( sum != lhs.extent(0) ) {
      std::cout << "Upper Tri Solve FAILURE" << std::endl;
    }
    else {
      std::cout << "Upper Tri Solve SUCCESS!" << std::endl;
      //std::cout << "Num-levels = " << kh->get_sptrsv_handle()->get_num_levels() << std::endl;
    }
    EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

    Kokkos::deep_copy(lhs, 0);
    kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHD_RP);
    kh.get_sptrsv_handle()->print_algorithm();
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    std::cout << "UTRI Solve SEQLVLSCHD_RP Time: " << timer.seconds() << std::endl;

    sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
    if ( sum != lhs.extent(0) ) {
      std::cout << "Upper Tri Solve FAILURE" << std::endl;
    }
    else {
      std::cout << "Upper Tri Solve SUCCESS!" << std::endl;
      //std::cout << "Num-levels = " << kh->get_sptrsv_handle()->get_num_levels() << std::endl;
    }
    EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

    Kokkos::deep_copy(lhs, 0);
    kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHED_TP2);
    kh.get_sptrsv_handle()->print_algorithm();
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    std::cout << "UTRI Solve SEQLVLSCHED_TP2 Time: " << timer.seconds() << std::endl;

    sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
    if ( sum != lhs.extent(0) ) {
      std::cout << "Upper Tri Solve FAILURE" << std::endl;
    }
    else {
      std::cout << "Upper Tri Solve SUCCESS!" << std::endl;
      //std::cout << "Num-levels = " << kh->get_sptrsv_handle()->get_num_levels() << std::endl;
    }
    EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );


    kh.destroy_sptrsv_handle();
  }

}
#endif


template < class ViewType, typename ValueType, typename OrdinalType >
struct ReductionCheck {

  using lno_t = OrdinalType;
  using value_type = ValueType;

  ViewType lhs;

  ReductionCheck( const ViewType & lhs_ ) : lhs(lhs_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( lno_t i, value_type &tsum ) const {
    tsum += lhs(i);
  }

};


template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_sptrsv() {

  typedef Kokkos::View< size_type*, device >  RowMapType;
  typedef Kokkos::View< lno_t*, device >      EntriesType;
  typedef Kokkos::View< scalar_t*, device >   ValuesType;

  scalar_t ZERO = scalar_t(0);
  scalar_t ONE = scalar_t(1);

  const size_type nrows = 5;
  const size_type nnz   = 10;

  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle <size_type, lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space, typename device::memory_space>;

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
  using host_crsmat_t = typename KernelHandle::SPTRSVHandleType::host_crsmat_t;
  using host_graph_t = typename host_crsmat_t::StaticCrsGraphType;

  using row_map_view_t  = typename host_graph_t::row_map_type::non_const_type;
  using cols_view_t     = typename host_graph_t::entries_type::non_const_type;
  using values_view_t   = typename host_crsmat_t::values_type::non_const_type;

  using in_row_map_view_t = typename RowMapType::HostMirror;
  using in_cols_view_t    = typename EntriesType::HostMirror;
  using in_values_view_t  = typename ValuesType::HostMirror;

  // L & U handle for supernodal SpTrsv
  KernelHandle khL;
  KernelHandle khU;

  // right-hand-side and solution
  ValuesType B ("rhs", nrows);
  ValuesType X ("sol", nrows);

  // host CRS for L & U
  host_crsmat_t L, U;
#endif

  // Upper tri
  {
    RowMapType  row_map("row_map", nrows+1);
    EntriesType entries("entries", nnz);
    ValuesType  values("values", nnz);

    auto hrow_map = Kokkos::create_mirror_view(row_map);
    auto hentries = Kokkos::create_mirror_view(entries);
    auto hvalues  = Kokkos::create_mirror_view(values);

    hrow_map(0) = 0;
    hrow_map(1) = 2;
    hrow_map(2) = 4;
    hrow_map(3) = 7;
    hrow_map(4) = 9;
    hrow_map(5) = 10;

    hentries(0) = 0;
    hentries(1) = 2;
    hentries(2) = 1;
    hentries(3) = 4;
    hentries(4) = 2;
    hentries(5) = 3;
    hentries(6) = 4;
    hentries(7) = 3;
    hentries(8) = 4;
    hentries(9) = 4;

    for ( size_type i = 0; i < nnz; ++i ) {
      hvalues(i) = ONE;
    }

    Kokkos::deep_copy(row_map, hrow_map);
    Kokkos::deep_copy(entries, hentries);
    Kokkos::deep_copy(values,  hvalues);

    // Create known_lhs, generate rhs, then solve for lhs to compare to known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
    crsMat_t triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);
    KokkosSparse::spmv( "N", ONE, triMtx, known_lhs, ZERO, rhs);

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    Kokkos::deep_copy (B, ONE);
    KokkosSparse::spmv ("N", ONE, triMtx, B, ZERO, X);
#endif

    {
      KernelHandle kh;
      bool is_lower_tri = false;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);

      sptrsv_symbolic( &kh, row_map, entries );
      Kokkos::fence();

      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      Kokkos::deep_copy(lhs, 0);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHD_RP);
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      //FIXME Issues with various integral type combos - algorithm currently unavailable and commented out until fixed
      /*
      Kokkos::deep_copy(lhs, 0);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHED_TP2);
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );
      */


      kh.destroy_sptrsv_handle();
    }

    {
      Kokkos::deep_copy(lhs, 0);
      KernelHandle kh;
      bool is_lower_tri = false;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN, nrows, is_lower_tri);
      auto chain_threshold = 1;
      kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);

      sptrsv_symbolic( &kh, row_map, entries );
      Kokkos::fence();

      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      kh.destroy_sptrsv_handle();
    }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if (std::is_same<size_type,int>::value && std::is_same<lno_t,int>::value && std::is_same<typename device::execution_space, Kokkos::Cuda>::value)
    {
      Kokkos::deep_copy(lhs, 0);
      KernelHandle kh;
      bool is_lower_tri = false;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows, is_lower_tri);

      sptrsv_symbolic(&kh, row_map, entries, values);
      Kokkos::fence();

      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      kh.destroy_sptrsv_handle();
    }
#endif

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    {
      // copy and save U for Supernodal Sptrsv
      size_type nnzL = hrow_map (nrows);
      row_map_view_t Urow_map ("rowmap_view", nrows+1);
      cols_view_t    Uentries ("colmap_view", nnzL);
      values_view_t  Uvalues  ("values_view", nnzL);
      Kokkos::deep_copy (Urow_map, hrow_map);
      Kokkos::deep_copy (Uentries, hentries);
      Kokkos::deep_copy (Uvalues,  hvalues);

      host_graph_t static_graph(Uentries, Urow_map);
      U = host_crsmat_t ("CrsMatrixU", nrows, Uvalues, static_graph);

      // create handle for Supernodal Sptrsv
      bool is_lower_tri = false;
      khU.create_sptrsv_handle (SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, is_lower_tri);
    }
#endif
  }

  // Lower tri
  {
    RowMapType  row_map("row_map", nrows+1);
    EntriesType entries("entries", nnz);
    ValuesType  values("values", nnz);

    auto hrow_map = Kokkos::create_mirror_view(row_map);
    auto hentries = Kokkos::create_mirror_view(entries);
    auto hvalues  = Kokkos::create_mirror_view(values);

    hrow_map(0) = 0;
    hrow_map(1) = 1;
    hrow_map(2) = 2;
    hrow_map(3) = 4;
    hrow_map(4) = 6;
    hrow_map(5) = 10;

    hentries(0) = 0;
    hentries(1) = 1;
    hentries(2) = 0;
    hentries(3) = 2;
    hentries(4) = 2;
    hentries(5) = 3;
    hentries(6) = 1;
    hentries(7) = 2;
    hentries(8) = 3;
    hentries(9) = 4;

    for ( size_type i = 0; i < nnz; ++i ) {
      hvalues(i) = ONE;
    }

    Kokkos::deep_copy(row_map, hrow_map);
    Kokkos::deep_copy(entries, hentries);
    Kokkos::deep_copy(values,  hvalues);

    // Create known_lhs, generate rhs, then solve for lhs to compare to known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
    crsMat_t triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);
    KokkosSparse::spmv( "N", ONE, triMtx, known_lhs, ZERO, rhs);

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    KokkosSparse::spmv( "N", ONE, triMtx, X, ZERO, B);
#endif

    {
      KernelHandle kh;
      bool is_lower_tri = true;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);

      sptrsv_symbolic( &kh, row_map, entries );
      Kokkos::fence();

      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      Kokkos::deep_copy(lhs, 0);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHD_RP);
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      //FIXME Issues with various integral type combos - algorithm currently unavailable and commented out until fixed
      /*
      Kokkos::deep_copy(lhs, 0);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHED_TP2);
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );
      */

      kh.destroy_sptrsv_handle();
    }

    {
      Kokkos::deep_copy(lhs, 0);
      KernelHandle kh;
      bool is_lower_tri = true;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN, nrows, is_lower_tri);
      auto chain_threshold = 1;
      kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);

      sptrsv_symbolic( &kh, row_map, entries );
      Kokkos::fence();

      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      kh.destroy_sptrsv_handle();
    }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if (std::is_same<size_type,int>::value && std::is_same<lno_t,int>::value && std::is_same<typename device::execution_space, Kokkos::Cuda>::value)
    {
      Kokkos::deep_copy(lhs, 0);
      KernelHandle kh;
      bool is_lower_tri = true;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows, is_lower_tri);

      sptrsv_symbolic(&kh, row_map, entries, values);
      Kokkos::fence();

      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );

      kh.destroy_sptrsv_handle();
    }
#endif

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    {
      // convert L into csc
      using host_exec_space = typename KernelHandle::SPTRSVHandleType::supercols_host_execution_space;

      size_type nnzL = hrow_map (nrows);
      row_map_view_t Lrow_map ("rowmap_view", nrows+1);
      cols_view_t    Lentries ("colmap_view", nnzL);
      values_view_t  Lvalues  ("values_view", nnzL);
      transpose_matrix <in_row_map_view_t, in_cols_view_t, in_values_view_t,
                           row_map_view_t,    cols_view_t,    values_view_t,
                        row_map_view_t, host_exec_space>
        (nrows, nrows, hrow_map, hentries, hvalues,
                       Lrow_map, Lentries, Lvalues);
      Kokkos::fence();

      // sptrsv assumes the diagonal "block" is stored before off-diagonal blocks
      for (size_type j = 0; j < nrows; j++) {
        for (size_type k = Lrow_map (j); k < Lrow_map (j+1); k++) {
          if (Lentries (k) == (lno_t)j) {
            scalar_t val = Lvalues (k);

            Lvalues (k) = Lvalues (Lrow_map(j));
            Lentries (k) = Lentries (Lrow_map(j));

            Lvalues (Lrow_map(j)) = val;
            Lentries (Lrow_map(j)) = j;
            break;
          }
        }
      }

      // store L in crsmat
      host_graph_t static_graph(Lentries, Lrow_map);
      L = host_crsmat_t ("CrsMatrixL", nrows, Lvalues, static_graph);

      bool is_lower_tri = true;
      khL.create_sptrsv_handle (SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, is_lower_tri);
    }

    {
      // unit-test for correctness
      // > set up supernodes (block size = one)
      int *etree = NULL;
      Kokkos::View<int*, Kokkos::HostSpace> supercols ("supercols", 1+nrows);
      for (size_type i = 0; i <= nrows; i++) {
        supercols (i) = i;
      }

      // invert diagonal blocks
      bool invert_diag = true;
      khL.set_sptrsv_invert_diagonal (invert_diag);
      khU.set_sptrsv_invert_diagonal (invert_diag);

      // > symbolic (on host)
      sptrsv_supernodal_symbolic (nrows, supercols.data (), etree, L.graph, &khL, U.graph, &khU);
      // > numeric (on host)
      sptrsv_compute (&khL, L);
      sptrsv_compute (&khU, U);
      Kokkos::fence();

      // > solve 
      Kokkos::deep_copy (X, 0);
      sptrsv_solve (&khL, &khU, X, B);
      Kokkos::fence();

      // > check
      scalar_t sum = 0.0;
      Kokkos::parallel_reduce (Kokkos::RangePolicy<typename device::execution_space>(0, X.extent(0)), ReductionCheck<ValuesType, scalar_t, lno_t>(X), sum);
      if ( sum != lhs.extent(0) ) {
        std::cout << "Supernode Tri Solve FAILURE" << std::endl;
        khL.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(X.extent(0)) );

      khL.destroy_sptrsv_handle();
      khU.destroy_sptrsv_handle();
    }
#endif
  }
}

} // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_sptrsv() {
  Test::run_test_sptrsv<scalar_t, lno_t, size_type, device>();
//  Test::run_test_sptrsv_mtx<scalar_t, lno_t, size_type, device>();
}


#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, sparse ## _ ## sptrsv ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_sptrsv<SCALAR,ORDINAL,OFFSET,DEVICE>(); \
}


#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif

#if 0

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

#endif
