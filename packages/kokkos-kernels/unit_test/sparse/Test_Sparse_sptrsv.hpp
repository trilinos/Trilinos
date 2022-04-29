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

#include <gtest/gtest.h>

using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Impl;
using namespace KokkosKernels::Experimental;

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #endif
// #ifndef kokkos_complex_float
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

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

namespace {
template <class ViewType, typename ValueType, typename OrdinalType>
struct ReductionCheck {
  using lno_t      = OrdinalType;
  using value_type = ValueType;

  ViewType lhs;

  ReductionCheck(const ViewType &lhs_) : lhs(lhs_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(lno_t i, value_type &tsum) const { tsum += lhs(i); }
};
}  // namespace

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void run_test_sptrsv() {
  typedef Kokkos::View<size_type *, device> RowMapType;
  typedef Kokkos::View<lno_t *, device> EntriesType;
  typedef Kokkos::View<scalar_t *, device> ValuesType;

  scalar_t ZERO = scalar_t(0);
  scalar_t ONE  = scalar_t(1);

  const size_type nrows = 5;
  const size_type nnz   = 10;

  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, typename device::execution_space,
      typename device::memory_space, typename device::memory_space>;

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
  using host_crsmat_t = typename KernelHandle::SPTRSVHandleType::host_crsmat_t;
  using host_graph_t  = typename host_crsmat_t::StaticCrsGraphType;

  using row_map_view_t = typename host_graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename host_graph_t::entries_type::non_const_type;
  using values_view_t  = typename host_crsmat_t::values_type::non_const_type;

  // L & U handle for supernodal SpTrsv
  KernelHandle khL;
  KernelHandle khU;

  // right-hand-side and solution
  ValuesType B("rhs", nrows);
  ValuesType X("sol", nrows);

  // host CRS for L & U
  host_crsmat_t L, U, Ut;
#endif

  // Upper tri
  {
    RowMapType row_map("row_map", nrows + 1);
    EntriesType entries("entries", nnz);
    ValuesType values("values", nnz);

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

    for (size_type i = 0; i < nnz; ++i) {
      hvalues(i) = ONE;
    }

    Kokkos::deep_copy(row_map, hrow_map);
    Kokkos::deep_copy(entries, hentries);
    Kokkos::deep_copy(values, hvalues);

    // Create known_lhs, generate rhs, then solve for lhs to compare to
    // known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
    crsMat_t triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);
    KokkosSparse::spmv("N", ONE, triMtx, known_lhs, ZERO, rhs);

    {
      KernelHandle kh;
      bool is_lower_tri = false;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows,
                              is_lower_tri);

      sptrsv_symbolic(&kh, row_map, entries);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      Kokkos::deep_copy(lhs, ZERO);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHD_RP);
      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      // FIXME Issues with various integral type combos - algorithm currently
      // unavailable and commented out until fixed
      /*
      Kokkos::deep_copy(lhs, ZERO);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHED_TP2);
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename
      device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType,
      scalar_t, lno_t>(lhs), sum); if ( sum != lhs.extent(0) ) { std::cout <<
      "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );
      */

      kh.destroy_sptrsv_handle();
    }

    {
      Kokkos::deep_copy(lhs, ZERO);
      KernelHandle kh;
      bool is_lower_tri = false;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN, nrows,
                              is_lower_tri);
      auto chain_threshold = 1;
      kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);

      sptrsv_symbolic(&kh, row_map, entries);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      kh.destroy_sptrsv_handle();
    }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if (std::is_same<size_type, int>::value &&
        std::is_same<lno_t, int>::value &&
        std::is_same<typename device::execution_space, Kokkos::Cuda>::value) {
      Kokkos::deep_copy(lhs, ZERO);
      KernelHandle kh;
      bool is_lower_tri = false;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows,
                              is_lower_tri);

      sptrsv_symbolic(&kh, row_map, entries, values);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Upper Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      kh.destroy_sptrsv_handle();
    }
#endif

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    const scalar_t FIVE    = scalar_t(5);
    const size_type nnz_sp = 14;
    {
      // U in csr
      row_map_view_t hUrowptr("hUrowptr", nrows + 1);
      cols_view_t hUcolind("hUcolind", nnz_sp);
      values_view_t hUvalues("hUvalues", nnz_sp);

      // rowptr
      hUrowptr(0) = 0;
      hUrowptr(1) = 4;
      hUrowptr(2) = 8;
      hUrowptr(3) = 11;
      hUrowptr(4) = 13;
      hUrowptr(5) = 14;

      // colind
      // first row (first supernode)
      hUcolind(0) = 0;
      hUcolind(1) = 1;
      hUcolind(2) = 2;
      hUcolind(3) = 4;
      // second row (first supernode)
      hUcolind(4) = 0;
      hUcolind(5) = 1;
      hUcolind(6) = 2;
      hUcolind(7) = 4;
      // third row (second supernode)
      hUcolind(8)  = 2;
      hUcolind(9)  = 3;
      hUcolind(10) = 4;
      // fourth row (third supernode)
      hUcolind(11) = 3;
      hUcolind(12) = 4;
      // fifth row (fourth supernode)
      hUcolind(13) = 4;

      // values
      // first row (first supernode)
      hUvalues(0) = FIVE;
      hUvalues(1) = ONE;
      hUvalues(2) = ONE;
      hUvalues(3) = ZERO;
      // second row (first supernode)
      hUvalues(4) = ZERO;
      hUvalues(5) = FIVE;
      hUvalues(6) = ZERO;
      hUvalues(7) = ONE;
      // third row (second supernode)
      hUvalues(8)  = FIVE;
      hUvalues(9)  = ONE;
      hUvalues(10) = ONE;
      // fourth row (third supernode)
      hUvalues(11) = FIVE;
      hUvalues(12) = ONE;
      // fifth row (fourth supernode)
      hUvalues(13) = FIVE;

      // save U for Supernodal Sptrsv
      host_graph_t static_graph(hUcolind, hUrowptr);
      U = host_crsmat_t("CrsMatrixU", nrows, hUvalues, static_graph);

      // create handle for Supernodal Sptrsv
      bool is_lower_tri = false;
      khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows,
                               is_lower_tri);

      // X = U*ONES to generate B = A*ONES (on device)
      {
        RowMapType Urowptr("Urowptr", nrows + 1);
        EntriesType Ucolind("Ucolind", nnz_sp);
        ValuesType Uvalues("Uvalues", nnz_sp);

        Kokkos::deep_copy(Urowptr, hUrowptr);
        Kokkos::deep_copy(Ucolind, hUcolind);
        Kokkos::deep_copy(Uvalues, hUvalues);

        crsMat_t mtxU("mtxU", nrows, nrows, nnz_sp, Uvalues, Urowptr, Ucolind);
        Kokkos::deep_copy(B, ONE);
        KokkosSparse::spmv("N", ONE, mtxU, B, ZERO, X);
      }
    }

    {
      // U in csc (for inverting off-diag)
      row_map_view_t hUcolptr("hUcolptr", nrows + 1);
      cols_view_t hUrowind("hUrowind", nnz_sp);
      values_view_t hUvalues("hUvalues", nnz_sp);

      // colptr
      hUcolptr(0) = 0;
      hUcolptr(1) = 2;
      hUcolptr(2) = 4;
      hUcolptr(3) = 7;
      hUcolptr(4) = 9;
      hUcolptr(5) = 14;

      // colind
      // first column (first supernode)
      hUrowind(0) = 0;
      hUrowind(1) = 1;
      // second column (first supernode)
      hUrowind(2) = 0;
      hUrowind(3) = 1;
      // third column (second supernode)
      hUrowind(4) = 2;
      hUrowind(5) = 0;
      hUrowind(6) = 1;
      // fourth column (third supernode)
      hUrowind(7) = 3;
      hUrowind(8) = 2;
      // fifth column (fourth supernode)
      hUrowind(9)  = 4;
      hUrowind(10) = 0;
      hUrowind(11) = 1;
      hUrowind(12) = 2;
      hUrowind(13) = 3;

      // values
      // first column (first supernode)
      hUvalues(0) = FIVE;
      hUvalues(1) = ZERO;
      // second column (first supernode)
      hUvalues(2) = ONE;
      hUvalues(3) = FIVE;
      // third column (second supernode)
      hUvalues(4) = FIVE;
      hUvalues(5) = ONE;
      hUvalues(6) = ZERO;
      // fourth column (third supernode)
      hUvalues(7) = FIVE;
      hUvalues(8) = ONE;
      // fifth column (fourth supernode)
      hUvalues(9)  = FIVE;
      hUvalues(10) = ZERO;
      hUvalues(11) = ONE;
      hUvalues(12) = ONE;
      hUvalues(13) = ONE;

      // store Ut in crsmat
      host_graph_t static_graph(hUrowind, hUcolptr);
      Ut = host_crsmat_t("CrsMatrixUt", nrows, hUvalues, static_graph);
    }
#endif
  }

  // Lower tri
  {
    RowMapType row_map("row_map", nrows + 1);
    EntriesType entries("entries", nnz);
    ValuesType values("values", nnz);

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

    for (size_type i = 0; i < nnz; ++i) {
      hvalues(i) = ONE;
    }

    Kokkos::deep_copy(row_map, hrow_map);
    Kokkos::deep_copy(entries, hentries);
    Kokkos::deep_copy(values, hvalues);

    // Create known_lhs, generate rhs, then solve for lhs to compare to
    // known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
    crsMat_t triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);
    KokkosSparse::spmv("N", ONE, triMtx, known_lhs, ZERO, rhs);

    {
      KernelHandle kh;
      bool is_lower_tri = true;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows,
                              is_lower_tri);

      sptrsv_symbolic(&kh, row_map, entries);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      Kokkos::deep_copy(lhs, ZERO);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHD_RP);
      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      // FIXME Issues with various integral type combos - algorithm currently
      // unavailable and commented out until fixed
      /*
      Kokkos::deep_copy(lhs, ZERO);
      kh.get_sptrsv_handle()->set_algorithm(SPTRSVAlgorithm::SEQLVLSCHED_TP2);
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
      Kokkos::fence();

      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<typename
      device::execution_space>(0, lhs.extent(0)), ReductionCheck<ValuesType,
      scalar_t, lno_t>(lhs), sum); if ( sum != lhs.extent(0) ) { std::cout <<
      "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE( sum == scalar_t(lhs.extent(0)) );
      */

      kh.destroy_sptrsv_handle();
    }

    {
      Kokkos::deep_copy(lhs, ZERO);
      KernelHandle kh;
      bool is_lower_tri = true;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN, nrows,
                              is_lower_tri);
      auto chain_threshold = 1;
      kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);

      sptrsv_symbolic(&kh, row_map, entries);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      kh.destroy_sptrsv_handle();
    }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if (std::is_same<size_type, int>::value &&
        std::is_same<lno_t, int>::value &&
        std::is_same<typename device::execution_space, Kokkos::Cuda>::value) {
      Kokkos::deep_copy(lhs, ZERO);
      KernelHandle kh;
      bool is_lower_tri = true;
      kh.create_sptrsv_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows,
                              is_lower_tri);

      sptrsv_symbolic(&kh, row_map, entries, values);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0,
                                                                lhs.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(lhs), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Lower Tri Solve FAILURE" << std::endl;
        kh.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(lhs.extent(0)));

      kh.destroy_sptrsv_handle();
    }
#endif

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    {
      // L in csc
      const scalar_t TWO     = scalar_t(2);
      const scalar_t FIVE    = scalar_t(5);
      const size_type nnz_sp = 14;

      row_map_view_t hLcolptr("hUcolptr", nrows + 1);
      cols_view_t hLrowind("hUrowind", nnz_sp);
      values_view_t hLvalues("hUvalues", nnz_sp);

      // colptr
      hLcolptr(0) = 0;
      hLcolptr(1) = 4;
      hLcolptr(2) = 8;
      hLcolptr(3) = 11;
      hLcolptr(4) = 13;
      hLcolptr(5) = 14;

      // rowind
      // first column (first supernode)
      hLrowind(0) = 0;
      hLrowind(1) = 1;
      hLrowind(2) = 2;
      hLrowind(3) = 4;
      // second column (first supernode)
      hLrowind(4) = 0;
      hLrowind(5) = 1;
      hLrowind(6) = 2;
      hLrowind(7) = 4;
      // third column (second supernode)
      hLrowind(8)  = 2;
      hLrowind(9)  = 3;
      hLrowind(10) = 4;
      // fourth column (third supernode)
      hLrowind(11) = 3;
      hLrowind(12) = 4;
      // fifth column (fourth supernode)
      hLrowind(13) = 4;

      // values
      // first column (first supernode)
      hLvalues(0) = FIVE;
      hLvalues(1) = TWO;
      hLvalues(2) = ONE;
      hLvalues(3) = ZERO;
      // second column (first supernode)
      hLvalues(4) = ZERO;
      hLvalues(5) = FIVE;
      hLvalues(6) = ZERO;
      hLvalues(7) = ONE;
      // third column (second supernode)
      hLvalues(8)  = FIVE;
      hLvalues(9)  = ONE;
      hLvalues(10) = ONE;
      // fourth column (third supernode)
      hLvalues(11) = FIVE;
      hLvalues(12) = ONE;
      // fifth column (fourth supernode)
      hLvalues(13) = FIVE;

      // store Lt in crsmat
      host_graph_t static_graph(hLrowind, hLcolptr);
      L = host_crsmat_t("CrsMatrixL", nrows, hLvalues, static_graph);

      bool is_lower_tri = true;
      khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows,
                               is_lower_tri);

      // generate B = A*ONES = L*(U*ONES), where X = U*ONES (on device)
      {
        RowMapType Lcolptr("Lcolptr", nrows + 1);
        EntriesType Lrowind("Lrowind", nnz_sp);
        ValuesType Lvalues("Lvalues", nnz_sp);

        Kokkos::deep_copy(Lcolptr, hLcolptr);
        Kokkos::deep_copy(Lrowind, hLrowind);
        Kokkos::deep_copy(Lvalues, hLvalues);

        crsMat_t mtxL("mtxL", nrows, nrows, nnz_sp, Lvalues, Lcolptr, Lrowind);
        KokkosSparse::spmv("T", ONE, mtxL, X, ZERO, B);
      }
    }

    {
      // unit-test for supernode SpTrsv (default)
      // > set up supernodes (block size = one)
      size_type nsupers = 4;
      Kokkos::View<int *, Kokkos::HostSpace> supercols("supercols",
                                                       1 + nsupers);
      supercols(0) = 0;
      supercols(1) = 2;     // two columns
      supercols(2) = 3;     // one column
      supercols(3) = 4;     // one column
      supercols(4) = 5;     // one column
      int *etree   = NULL;  // we generate graph internally

      // invert diagonal blocks
      bool invert_diag = true;
      khL.set_sptrsv_invert_diagonal(invert_diag);
      khU.set_sptrsv_invert_diagonal(invert_diag);

      // > symbolic (on host)
      sptrsv_supernodal_symbolic(nsupers, supercols.data(), etree, L.graph,
                                 &khL, U.graph, &khU);
      // > numeric (on host)
      sptrsv_compute(&khL, L);
      sptrsv_compute(&khU, U);
      Kokkos::fence();

      // > solve
      ValuesType b("b", nrows);
      Kokkos::deep_copy(b, B);
      Kokkos::deep_copy(X, ZERO);
      sptrsv_solve(&khL, &khU, X, b);
      Kokkos::fence();

      // > check
      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0, X.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(X), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Supernode Tri Solve FAILURE : " << sum << " vs."
                  << lhs.extent(0) << std::endl;
        khL.get_sptrsv_handle()->print_algorithm();
      } else {
        std::cout << "Supernode Tri Solve SUCCESS" << std::endl;
        khL.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(X.extent(0)));

      khL.destroy_sptrsv_handle();
      khU.destroy_sptrsv_handle();
    }

    {
      // unit-test for supernode SpTrsv (running TRMM on device for compute)
      // > set up supernodes
      size_type nsupers = 4;
      Kokkos::View<int *, Kokkos::HostSpace> supercols("supercols",
                                                       1 + nsupers);
      supercols(0) = 0;
      supercols(1) = 2;     // two columns
      supercols(2) = 3;     // one column
      supercols(3) = 4;     // one column
      supercols(4) = 5;     // one column
      int *etree   = NULL;  // we generate tree internally

      // > create handles
      KernelHandle khLd;
      KernelHandle khUd;
      khLd.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, true);
      khUd.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, false);

      // > invert diagonal blocks
      bool invert_diag = true;
      khLd.set_sptrsv_invert_diagonal(invert_diag);
      khUd.set_sptrsv_invert_diagonal(invert_diag);

      // > invert off-diagonal blocks
      bool invert_offdiag = true;
      khUd.set_sptrsv_column_major(true);
      khLd.set_sptrsv_invert_offdiagonal(invert_offdiag);
      khUd.set_sptrsv_invert_offdiagonal(invert_offdiag);

      // > forcing sptrsv compute to perform TRMM on device
      khLd.set_sptrsv_diag_supernode_sizes(1, 1);
      khUd.set_sptrsv_diag_supernode_sizes(1, 1);

      // > symbolic (on host)
      sptrsv_supernodal_symbolic(nsupers, supercols.data(), etree, L.graph,
                                 &khLd, Ut.graph, &khUd);
      // > numeric (on host)
      sptrsv_compute(&khLd, L);
      sptrsv_compute(&khUd, Ut);
      Kokkos::fence();

      // > solve
      ValuesType b("b", nrows);
      Kokkos::deep_copy(b, B);
      Kokkos::deep_copy(X, ZERO);
      sptrsv_solve(&khLd, &khUd, X, b);
      Kokkos::fence();

      // > check
      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<typename device::execution_space>(0, X.extent(0)),
          ReductionCheck<ValuesType, scalar_t, lno_t>(X), sum);
      if (sum != lhs.extent(0)) {
        std::cout << "Supernode Tri Solve FAILURE : " << sum << " vs."
                  << lhs.extent(0) << std::endl;
        khLd.get_sptrsv_handle()->print_algorithm();
      } else {
        std::cout << "Supernode Tri Solve SUCCESS" << std::endl;
        khLd.get_sptrsv_handle()->print_algorithm();
      }
      EXPECT_TRUE(sum == scalar_t(X.extent(0)));

      khLd.destroy_sptrsv_handle();
      khUd.destroy_sptrsv_handle();
    }
#endif
  }
}

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_sptrsv() {
  Test::run_test_sptrsv<scalar_t, lno_t, size_type, device>();
  //  Test::run_test_sptrsv_mtx<scalar_t, lno_t, size_type, device>();
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                      \
  TEST_F(TestCategory,                                                     \
         sparse##_##sptrsv##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_sptrsv<SCALAR, ORDINAL, OFFSET, DEVICE>();                        \
  }

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&         \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

#undef EXECUTE_TEST
