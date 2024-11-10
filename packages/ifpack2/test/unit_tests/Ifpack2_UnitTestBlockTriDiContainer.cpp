// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestBlockTriDiContainer.cpp

\brief Ifpack2 Unit and performance test for the BlockTriDiContainter template.
*/


#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_Core.hpp>

#ifdef HAVE_TEUCHOS_COMPLEX
# include <complex>
#endif
#include <memory>

#include "Ifpack2_UnitTestBlockCrsUtil.hpp"
#include "Ifpack2_UnitTestBlockTriDiContainerUtil.hpp"

struct Input {
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  bool quiet, teuchos_test, contiguous, verbose;
  int isplit, jsplit;
  int ni, nj, nk;
  int bs; // block size
  int nrhs; // #vectors in multivector
  bool nonuniform_lines; // Test lines having different sizes.
  int repeat, iterations, check_tol_every;
  double tol;
  bool use_seq_method, use_overlap_comm, jacobi;

  Input (const Teuchos::RCP<const Teuchos::Comm<int> >& icomm) { init(icomm); }

  Input (int argc, char** argv, const Teuchos::RCP<const Teuchos::Comm<int> >& icomm) {
    init(icomm);
    Teuchos::CommandLineProcessor clp(false);
    clp.setDocString("Unit and performance tester for BlockTriDiContainer. The performance test\n"
                     "constructs a 3D 7-point stencil in a block domain of arbitrary size. The\n"
                     "block is distributed by splitting it in each axis. The index space is\n"
                     "denoted (I,J,K) in the documentation. The lines are in the K dimension.");
    clp.setOption("ni", &ni, "Number of cells in the I dimension.");
    clp.setOption("nj", &nj, "Number of cells in the J dimension.");
    clp.setOption("nk", &nk, "Number of cells in the K dimension.");
    clp.setOption("isplit", &isplit, "Cut the I dimension into this many pieces for MPI distribution. "
                  "isplit*jsplit must equal the number of ranks.");
    clp.setOption("jsplit", &jsplit, "Cut the J dimension into this many pieces for MPI distribution. "
                  "isplit*jsplit must equal the number of ranks.");
    clp.setOption("bs", &bs, "Block size. This might be the number of degrees of freedom in a system of PDEs, "
                  "for example.");
    clp.setOption("nrhs", &nrhs, "Number of right hand sides to solve for.");
    clp.setOption("repeat", &repeat, "Number of times to run the performance test for timing. "
                  "To run the performance test, set it to a number > 0.");
    clp.setOption("iterations", &iterations, "Max number of fixed-point iterations.");
    clp.setOption("tol", &tol, "Tolerance for norm-based termination. Set to <= 0 to turn off norm-based termination.");
    clp.setOption("check-tol-every", &check_tol_every, "Check norm every CE iterations.");
    clp.setOption("verbose", "quiet", &verbose, "Verbose output.");
    clp.setOption("test", "notest", &teuchos_test, "Run unit tests.");
    clp.setOption("jacobi", "tridiag", &jacobi, "Run with little-block Jacobi instead of block tridagional.");
    clp.setOption("contiguous", "noncontiguous", &contiguous, "Use contiguous GIDs.");
    clp.setOption("seq-method", "non-seq-method", &use_seq_method, "Developer option.");
    clp.setOption("overlap-comm", "non-overlap-comm", &use_overlap_comm, "Developer option.");
    auto out = clp.parse(argc, argv, icomm->getRank() == 0 ? &std::cerr : nullptr);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(out > Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED,
                                "Parse error.");
    check();
  }

  void check () {
    TEUCHOS_ASSERT( ! comm.is_null());
    TEUCHOS_ASSERT(ni >= 1 && nj >= 1);
    TEUCHOS_ASSERT(isplit >= 1 && jsplit >= 1);
    TEUCHOS_ASSERT(bs >= 1);
    TEUCHOS_ASSERT(nrhs >= 1);
    TEUCHOS_ASSERT(repeat >= 0);
    TEUCHOS_ASSERT(iterations >= 1);
    TEUCHOS_ASSERT(check_tol_every >= 1);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(nk <= 1, "k dimension is <= 1; must be >= 2.");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
      isplit*jsplit != comm->getSize(),
      "isplit*jsplit must be == to #ranks; isplit = " << isplit << " jsplit = " << jsplit <<
      " but #ranks = " << comm->getSize() << ".");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(std::max(isplit, jsplit) > comm->getSize(),
                                "i,jsplit must be <= #ranks; isplit = " << isplit <<
                                " jsplit = " << jsplit);
    if (comm->getRank() == 0) print(std::cout);
  }

  void print (std::ostream& os) const {
    os << "<I> nranks " << comm->getSize()
       << " ni " << ni << " nj " << nj << " nk " << nk
       << " bs " << bs
       << " nrhs " << nrhs
       << " isplit " << isplit << " jsplit " << jsplit;
#ifdef KOKKOS_ENABLE_OPENMP
    os << " nthreads " << omp_get_max_threads();
#endif
    if (nonuniform_lines) os << " nonuniform-lines";
    if ( ! contiguous) os << " noncontiguous";
    if (use_seq_method) os << " seq";
    if (use_overlap_comm) os << " overlap";
    if (jacobi) os << " jacobi";
    os << "\n";
  }

private:
  void init (const Teuchos::RCP<const Teuchos::Comm<int> >& icomm) {
    comm = icomm;
    quiet = false;
    teuchos_test = true;
    ni = nj = nk = 10;
    isplit = comm->getSize();
    jsplit = 1;
    bs = 5;
    nrhs = 1;
    nonuniform_lines = false;
    repeat = 0;
    iterations = 20;
    check_tol_every = 1;
    tol = 0;
    contiguous = true;
    verbose = false;
    use_seq_method = false;
    use_overlap_comm = false;
    jacobi = false;
  }
};

#ifdef HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
// Configure with Ifpack2_ENABLE_BlockTriDiContainer_Timers:BOOL=ON to get timer
// output for each piece of the computation.
template <typename Scalar,
          typename LO = Tpetra::Map<>::local_ordinal_type,
          typename GO = Tpetra::Map<>::global_ordinal_type>
static void run_performance_test (const Input& in) {
  typedef LO Int;
  typedef tif_utest::BlockCrsMatrixMaker<Scalar, LO, GO> bcmm;
  typedef tif_utest::BlockTriDiContainerTester<Scalar, LO, GO> btdct;

  // Performance test. Use BlockTriDiContainer directly.
  typename bcmm::StructuredBlock sb(in.ni, in.nj, in.nk, in.contiguous);
  typename bcmm::StructuredBlockPart sbp = bcmm::make_StructuredBlockPart(sb, in.isplit, in.jsplit,
                                                                          in.comm->getRank());
  auto g = bcmm::make_crs_graph(in.comm, sb, sbp);
  auto A = bcmm::make_bcrs_matrix(sb, g, in.bs);
  const auto T = btdct::make_BTDC(sb, sbp, A, in.use_overlap_comm, in.nonuniform_lines, in.jacobi,
                                  in.use_seq_method);
  const auto X = bcmm::make_multivector(sb, A, in.bs, in.nrhs);
  const auto B = Teuchos::rcp(new typename bcmm::Tpetra_MultiVector(A->getRangeMap(), in.nrhs));
  A->apply(*X, *B);
  const auto X_solve = Teuchos::rcp(new typename bcmm::Tpetra_MultiVector(A->getDomainMap(), in.nrhs));
  Teuchos::TimeMonitor::zeroOutTimers();
  T->initialize();
  auto input = T->createDefaultApplyParameters();
  input.zeroStartingSolution = true;
  input.maxNumSweeps = in.iterations;
  input.tolerance = in.tol;
  input.checkToleranceEvery = in.check_tol_every;
  int sweep = 0;
  {
    TEUCHOS_FUNC_TIME_MONITOR("Performance test");
    for (Int repeat = 0; repeat < in.repeat; ++repeat) {
      T->compute();
      sweep = T->applyInverseJacobi(*B, *X_solve, input);
    }
  }
  const auto rd = bcmm::reldif(*X, *X_solve);
  if (in.verbose && in.comm->getRank() == 0)
    std::cout << "Performance test rd = " << rd << "\n";
  if (in.comm->getRank() == 0) {
    const auto n0 = T->getNorms0();
    const auto nf = T->getNormsFinal();
    bool ok = true;
    if (n0 > 0) {
      std::stringstream ss;
      ss << "Norm reduction:";
      const auto r = nf / n0;
      ss << " " << r;
      if (r > in.tol && sweep != in.iterations) ok = false;
      ss << "\n";
      if (in.verbose)
        std::cout << ss.str();
      if ( ! ok)
        std::cout << "FAIL: norm-based termination\n";
    }
  }
  Teuchos::TimeMonitor::summarize();
}

#define TEUCHOS_TEST(result, msg) do {                          \
    const bool TEUCHOS_TEST_result = (result);                  \
    out << "TEST: " << msg << "\n";                             \
    out << TEUCHOS_PASS_FAIL(TEUCHOS_TEST_result) << "\n";      \
    if ( ! TEUCHOS_TEST_result) success = false;                \
  } while (0)

template <typename Scalar, typename LO, typename GO>
static LO run_teuchos_tests (const Input& in, Teuchos::FancyOStream& out, bool& success) {
  typedef LO Int;
  typedef tif_utest::BlockCrsMatrixMaker<Scalar, LO, GO> bcmm;
  typedef tif_utest::BlockTriDiContainerTester<Scalar, LO, GO> btdct;

  success = true;
  Int nerr = 0;
  Int ne = 0;
  const bool different_maps = false;
  for (const Int bs : {5, 1}) {
    for (const bool contiguous : {true, false}) {
#ifndef HAVE_TPETRA_BCRS_DO_POINT_IMPORT
      // BlockMultiVector's doImport seems not to support noncontiguous GIDs. But
      // if Trilinos is configured with Tpetra_BCRS_Point_Import:BOOL=ON, then
      // this is OK.
      if ( ! contiguous) continue;
#endif
      {
        ne = bcmm::test_StructuredBlockPart(20, 20, contiguous);
        TEUCHOS_TEST(ne == 0, "test_StructuredBlockPart");
        nerr += ne;
      }
      typename bcmm::StructuredBlock sb(2*in.isplit, 2*in.jsplit, 3, contiguous);
      typename bcmm::StructuredBlockPart sbp = bcmm::make_StructuredBlockPart(sb, in.isplit, in.jsplit,
                                                                              in.comm->getRank());
      {
        ne = bcmm::test_bcrs_matrix(in.comm, sb, sbp, bs, false);
        TEUCHOS_TEST(ne == 0, "test_bcrs_matrix");
        nerr += ne;
        ne = bcmm::test_bcrs_matrix(in.comm, sb, sbp, bs, true /* tridiags_only */);
        TEUCHOS_TEST(ne == 0, "test_bcrs_matrix tridiags only");
        nerr += ne;
      }
      for (const bool jacobi : {false, true})
        for (const bool seq_method : {false, true}) {
          //The sequential method only works when 'device' memory space is host-accessible (aka UVM semantics)
          if(seq_method && !Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace,
              typename bcmm::DeviceType::memory_space>::accessible) {
            continue;
          }
          for (const bool overlap_comm : {false, true}) { // temporary disabling overlap comm version
            if (seq_method && overlap_comm) continue;
            for (const bool nonuniform_lines : {false, true}) {
              for (const bool pointwise : {false, true}) {
                for (const bool explicitConversion : {false, true}) {
                  if (jacobi && nonuniform_lines) continue;
                  if (!pointwise && explicitConversion) continue;
                  for (const int nvec : {1, 3}) {
                    std::stringstream ss;
                    ss << "test_BR_BTDC:"
                      << " bs " << bs
                      << (contiguous ? " contig" : " noncontig")
                      << (jacobi ? " jacobi" : " tridiag")
                      << (seq_method ? " seq_method" : "")
                      << (overlap_comm ? " overlap_comm" : "")
                      << (pointwise ? " point_wise" : "")
                      << (explicitConversion ? " explicit_block_conversion" : "")
                      << (nonuniform_lines ? " nonuniform_lines" : " uniform_lines")
                      << " nvec " << nvec;
                    const std::string details = ss.str();
                    bool threw = false;
                    try {
                      ne = btdct::test_BR_BTDC(in.comm, sb, sbp, bs, nvec, nonuniform_lines,
                                              different_maps, jacobi, overlap_comm, seq_method, pointwise,
                                              explicitConversion, details);
                      nerr += ne;
                    } catch (const std::exception& e) {
                      threw = true;
                    }
                    if (threw)
                      printf("Exception threw from rank %d, %s\n", in.comm->getRank(), details.c_str());

                    TEUCHOS_TEST(ne == 0 && ! threw, details);
                  }
                }
              }
            }
          }
        }
    }
  }
  return nerr;
}
#endif // HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2BlockTriDi, Unit, Scalar, LO, GO) {
#ifdef HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
  Input in(Tpetra::getDefaultComm());
  run_teuchos_tests<Scalar, LO, GO>(in, out, success);
#endif // HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
}

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(Ifpack2BlockTriDi, Unit, SC, LO, GO )
#include "Ifpack2_ETIHelperMacros.h"
IFPACK2_ETI_MANGLING_TYPEDEFS()
IFPACK2_INSTANTIATE_SLG(UNIT_TEST_GROUP_SC_LO_GO)

static Teuchos::RCP<const Teuchos::Comm<int> > make_comm () {
  return Tpetra::getDefaultComm();
}

int main (int argc, char** argv) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(IFPACK2_BLOCKTRIDICONTAINER_ENABLE_PROFILE)
  cudaProfilerStop();
#endif

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  int ret = 0;
  do {
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    std::shared_ptr<Input> in;
    try {
      comm = make_comm();
      in = std::make_shared<Input>(argc, argv, comm);
    } catch (const std::exception& e) {
      if (comm.is_null() || comm->getRank() == 0) {
        std::cerr << e.what() << "\n";
      }
      ret = -1;
      break;
    }
    try {
      if (comm->getRank() == 0) {
        const bool detail = false;
        Kokkos::print_configuration(std::cout, detail);
      }
      if (in->teuchos_test) {
        Teuchos::UnitTestRepository::runUnitTestsFromMain(1, argv);
      }
#ifdef HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
      if (in->repeat) {
        run_performance_test<double>(*in);
      }
#endif
    } catch (const std::exception& e) {
      std::cerr << e.what() << "\n";
      ret = -1;
      break;
    }
  } while (0);

  return ret;
}
