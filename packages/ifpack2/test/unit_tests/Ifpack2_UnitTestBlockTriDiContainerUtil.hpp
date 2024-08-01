// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestBlockTriDiContainerUtil.hpp

\brief Ifpack2 Unit and performance test utilities for BlockTriDiContainter.
*/

#ifndef IFPACK2_UNITTEST_BLOCKTRIDICONTAINER_UTIL_HPP
#define IFPACK2_UNITTEST_BLOCKTRIDICONTAINER_UTIL_HPP

#include <Ifpack2_BlockRelaxation.hpp>
#ifdef HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
#include <Ifpack2_BlockTriDiContainer_decl.hpp>
#include "Tpetra_BlockCrsMatrix_Helpers_decl.hpp"

namespace tif_utest {

template <typename Scalar, typename LO, typename GO>
struct BlockTriDiContainerTester {
  typedef LO Int;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
#if defined(HAVE_IFPACK2_BLOCKTRIDICONTAINER_SMALL_SCALAR) 
  typedef typename std::conditional<std::is_same<Magnitude,double>::value,float,Magnitude>::type SmallMagnitude;
#else 
  typedef Magnitude SmallMagnitude;  
#endif
  typedef tif_utest::BlockCrsMatrixMaker<Scalar, LO, GO> bcmm;
  typedef typename bcmm::Tpetra_RowMatrix Tpetra_RowMatrix;
  typedef typename bcmm::Tpetra_BlockCrsMatrix Tpetra_BlockCrsMatrix;
  typedef typename bcmm::Tpetra_MultiVector Tpetra_MultiVector;
  typedef typename bcmm::StructuredBlock StructuredBlock;
  typedef typename bcmm::StructuredBlockPart StructuredBlockPart;
  typedef typename bcmm::StencilShape StencilShape;

  static void
  zero_block (Tpetra_BlockCrsMatrix& A, const LO& row_lid, const LO& row_lid_to_match) {
    const auto& g = A.getCrsGraph();
    const auto row_map = g.getRowMap();
    const auto col_map = g.getColMap();
    const auto gid = row_map->getGlobalElement(row_lid_to_match);
    const auto col_lid = col_map->getLocalElement(gid);
    auto block = A.getLocalBlockHostNonConst(row_lid, col_lid);
    const Int bs = block.extent(1);
    for (Int bi = 0; bi < bs; ++bi)
      for (Int bj = 0; bj < bs; ++bj)
        block(bi,bj) = 0;
  }

  // Template on array type because the PL requires ArrayRCP<LO> while the
  // container ctor requires Array<LO>.
  template <typename Array>
  static void
  make_parts (const StructuredBlock& sb, const StructuredBlockPart& sbp,
              Tpetra_BlockCrsMatrix& A, const bool nonuniform_lines, const bool jacobi,
              Teuchos::Array<Array>& parts) {
    const auto& g = A.getCrsGraph();
    if ( ! jacobi) {
      // For the parts.
      if ( ! nonuniform_lines) {
        const Int n_lcl_parts = (sbp.ie - sbp.is)*(sbp.je - sbp.js);
        parts.resize(n_lcl_parts);
        const auto rmap = g.getRowMap();
        for (Int i = sbp.is, ip = 0; i < sbp.ie; ++i)
          for (Int j = sbp.js; j < sbp.je; ++j, ++ip) {
            auto& part = parts[ip];
            part.resize(sbp.ke - sbp.ks);
            for (Int k = sbp.ks; k < sbp.ke; ++k) {
              const auto gid = sb.ijk2id(i,j,k);
              part[k - sbp.ks] = rmap->getLocalElement(gid);
            }
          }
      } else {
        // Break some lines into 2. This requires modifying A in order to make the
        // 2 line solves equivalent to 1.
        const auto rmap = g.getRowMap();
        for (Int i = sbp.is, ip = 0, cnt = 0; i < sbp.ie; ++i) {
          for (Int j = sbp.js; j < sbp.je; ++j, ++ip, ++cnt) {
            if (cnt % 3 == 0) {
              Int ks = sbp.ks, ke = ks + (sbp.ke - sbp.ks)/2;
              LO lids[2];
              for (Int pi = 0; pi < 2; ++pi) {
                parts.push_back(Array());
                auto& part = parts.back();
                part.resize(ke - ks);
                for (Int k = ks; k < ke; ++k) {
                  const auto gid = sb.ijk2id(i,j,k);
                  part[k - ks] = rmap->getLocalElement(gid);
                }
                lids[pi] = pi == 0 ? part[part.size() - 1] : part[0];
                ks = ke;
                ke = sbp.ke;
              }
              // Zero the corresponding offdiags so the solve is exact in the
              // test.
              zero_block(A, lids[0], lids[1]);
              zero_block(A, lids[1], lids[0]);
            } else {
              parts.push_back(Array());
              auto& part = parts.back();
              part.resize(sbp.ke - sbp.ks);
              for (Int k = sbp.ks; k < sbp.ke; ++k) {
                const auto gid = sb.ijk2id(i,j,k);
                part[k - sbp.ks] = rmap->getLocalElement(gid);
              }
            }
          }
        }
      }
    }
  }

  // Make BlockRelaxation smoother with BlockTriDiContainer
  // with a pointwise matrix.
  // N.B. Modifies A if nonuniform_lines is true.
  static Teuchos::RCP<Ifpack2::BlockRelaxation<Tpetra_RowMatrix> >
  make_BR_BTDC_PW (const StructuredBlock& sb, const StructuredBlockPart& sbp,
                const Teuchos::RCP<Tpetra_BlockCrsMatrix>& A,
                const bool nonuniform_lines = false,
                const bool zero_starting_soln = true,
                const int num_sweeps = 1,
                const bool jacobi = false,
                const bool explicitConversion = false) {
    Teuchos::Array<Teuchos::ArrayRCP<LO> > parts;
    // make_parts modifies entries of A so the call to convertToCrsMatrix
    // needs to happen after make_parts
    make_parts(sb, sbp, *A, nonuniform_lines, jacobi, parts);
    auto A_pw = Tpetra::convertToCrsMatrix(*A);
    const auto T = Teuchos::rcp(new Ifpack2::BlockRelaxation<Tpetra_RowMatrix>(A_pw));
    {
      Teuchos::ParameterList p;
      p.set("relaxation: container", "BlockTriDi");
      p.set("relaxation: type", "MT Split Jacobi");
      p.set("relaxation: sweeps", 1);
      p.set("partitioner: type", "user");
      p.set("relaxation: zero starting solution", zero_starting_soln);
      p.set<int>("relaxation: sweeps", num_sweeps);
      p.set<LO>("partitioner: local parts", parts.size());
      p.set("partitioner: parts", parts);
      p.set("partitioner: subparts per part", 1);
      p.set("partitioner: block size", A->getBlockSize());
      p.set("partitioner: explicit convert to BlockCrs", explicitConversion);
      T->setParameters(p);
    }
    return T;
  }

  // Make BlockRelaxation smoother with BlockTriDiContainer.
  // N.B. Modifies A if nonuniform_lines is true.
  static Teuchos::RCP<Ifpack2::BlockRelaxation<Tpetra_RowMatrix> >
  make_BR_BTDC (const StructuredBlock& sb, const StructuredBlockPart& sbp,
                const Teuchos::RCP<Tpetra_BlockCrsMatrix>& A,
                const bool nonuniform_lines = false,
                const bool zero_starting_soln = true,
                const int num_sweeps = 1,
                const bool jacobi = false) {
    const auto T = Teuchos::rcp(new Ifpack2::BlockRelaxation<Tpetra_RowMatrix>(A));
    {
      Teuchos::ParameterList p;
      p.set("relaxation: container", "BlockTriDi");
      p.set("relaxation: type", "MT Split Jacobi");
      p.set("relaxation: sweeps", 1);
      p.set("partitioner: type", "user");
      p.set("relaxation: zero starting solution", zero_starting_soln);
      p.set<int>("relaxation: sweeps", num_sweeps);
      Teuchos::Array<Teuchos::ArrayRCP<LO> > parts;
      make_parts(sb, sbp, *A, nonuniform_lines, jacobi, parts);
      p.set<LO>("partitioner: local parts", parts.size());
      p.set("partitioner: parts", parts);
      p.set("partitioner: subparts per part", 1);
      p.set("partitioner: block size", -1);
      p.set("partitioner: explicit convert to BlockCrs", false);
      T->setParameters(p);
    }
    return T;
  }

  // Make a bare BlockTriDiContainer
  // with a pointwise matrix.
  // N.B. Modifies A if nonuniform_lines is true.
  static Teuchos::RCP<Ifpack2::BlockTriDiContainer<Tpetra_RowMatrix> >
  make_BTDC_PW (const StructuredBlock& sb, const StructuredBlockPart& sbp,
             const Teuchos::RCP<Tpetra_BlockCrsMatrix>& A,
             const bool overlap_comm = false, const bool nonuniform_lines = false,
             const bool jacobi = false, const bool seq_method = false,
             const bool explicitConversion = false) {
    Teuchos::Array<Teuchos::Array<LO> > parts;
    // make_parts modifies entries of A so the call to convertToCrsMatrix
    // needs to happen after make_parts
    make_parts(sb, sbp, *A, nonuniform_lines, jacobi, parts);
    auto A_pw = Tpetra::convertToCrsMatrix(*A);

    return Teuchos::rcp(new Ifpack2::BlockTriDiContainer<Tpetra_RowMatrix>(
                          A_pw, parts, 1, overlap_comm, seq_method, A->getBlockSize(), explicitConversion));
  }

  // Make a bare BlockTriDiContainer.
  // N.B. Modifies A if nonuniform_lines is true.
  static Teuchos::RCP<Ifpack2::BlockTriDiContainer<Tpetra_RowMatrix> >
  make_BTDC (const StructuredBlock& sb, const StructuredBlockPart& sbp,
             const Teuchos::RCP<Tpetra_BlockCrsMatrix>& A,
             const bool overlap_comm = false, const bool nonuniform_lines = false,
             const bool jacobi = false, const bool seq_method = false) {
    Teuchos::Array<Teuchos::Array<LO> > parts;
    make_parts(sb, sbp, *A, nonuniform_lines, jacobi, parts);
    return Teuchos::rcp(new Ifpack2::BlockTriDiContainer<Tpetra_RowMatrix>(
                          A, parts, 1, overlap_comm, seq_method));
  }

  static Int
  test_BR_BTDC (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                const StructuredBlock& sb, const StructuredBlockPart& sbp,
                const Int bs, const Int nvec, const bool nonuniform_lines,
                const bool different_maps, const bool jacobi, const bool overlap_comm,
                const bool seq_method, const bool pointwise, const bool explicitConversion, 
                const std::string& details) {
#define TEST_BR_BTDC_FAIL(msg) do {             \
      ++nerr;                                   \
      if (comm->getRank() == 0) {               \
        std::stringstream _ss_;                 \
        _ss_ << msg << "\n";                    \
        std::cerr << _ss_.str();                \
      }                                         \
    } while (0)
#define TEST_BR_BTDC_SUCCESS(msg) do {          \
      if (comm->getRank() == 0) {               \
        std::stringstream _ss_;                 \
        _ss_ << msg << "\n";                    \
        std::cerr << _ss_.str();                \
      }                                         \
    } while (0)
    struct Parameters {
      const char* name;
      bool tridiag_only, tridiag_is_identity;
    };
    static const Parameters parms[] = {    
      {"x - R y", false, true }, // Test the matvec. D = I, so inv(D) (b - R x) = x - R y.
      {"solve"  , true , false}, // Test the line solve. R = 0, so inv(D) (b - R x) = inv(D) b.
      {"general", false, false}, // Test the general case.
      {"I"      , true , true }, // Internal test of tridiag_is_identity.
    };
    int nerr = 0;
    for (size_t pi = 0; pi < sizeof(parms)/sizeof(Parameters); ++pi) {
      const auto& p = parms[pi];
      auto g = bcmm::make_crs_graph(comm, sb, sbp, p.tridiag_only, different_maps);
      auto A = bcmm::make_bcrs_matrix(sb, g, bs, p.tridiag_is_identity, jacobi && p.tridiag_only);
      const bool solve_case = ! (p.tridiag_is_identity || p.tridiag_only);
      const bool zero_starting = solve_case;
      const int num_sweeps = solve_case ? (jacobi ? 40 : 20) : 1;
      const bool use_br = ! (overlap_comm || seq_method);
      const Magnitude tol = 1e-3;
      const auto T_br = use_br ?
        ( pointwise ?
          make_BR_BTDC_PW(sb, sbp, A, nonuniform_lines, zero_starting, num_sweeps, jacobi, explicitConversion): 
          make_BR_BTDC(sb, sbp, A, nonuniform_lines, zero_starting, num_sweeps, jacobi) 
        ):
        Teuchos::null;
      const auto T_bare = use_br ? Teuchos::null :
        ( pointwise ?
          make_BTDC_PW(sb, sbp, A, overlap_comm, nonuniform_lines, jacobi, seq_method, explicitConversion): 
          make_BTDC(sb, sbp, A, overlap_comm, nonuniform_lines, jacobi, seq_method) 
        );
      if ( ! T_br.is_null()) {
        T_br->initialize();
        T_br->compute();
      } else {
        T_bare->initialize();
        T_bare->compute(T_bare->createDefaultComputeParameters());
      }
      auto apply = [=] (const Tpetra_MultiVector& B, Tpetra_MultiVector& X,
                        const bool norm_based) -> int {
        if ( ! T_br.is_null()) {
          T_br->apply(B, X);
          return 0;
        } else {
          auto input = T_bare->createDefaultApplyParameters();
          input.zeroStartingSolution = zero_starting;
          input.maxNumSweeps = num_sweeps;
          input.tolerance = norm_based ? tol : 0;
          return T_bare->applyInverseJacobi(B, X, input);
        }
      };
      const auto X = bcmm::make_multivector(sb, A, bs, nvec);
      const auto B = Teuchos::rcp(new Tpetra_MultiVector(A->getRangeMap(), nvec));
      A->apply(*X, *B);
      const auto X_solve = Teuchos::rcp(new Tpetra_MultiVector(A->getDomainMap(), nvec));
      Magnitude rd = 0;
      if (p.tridiag_only && p.tridiag_is_identity) {
        // Test that we formed I.
        rd = bcmm::reldif(*X, *B);
        if (rd != 0) TEST_BR_BTDC_FAIL("FAIL: test_BR_BTDC (A = I) " << details << " rd " << rd);
        else TEST_BR_BTDC_SUCCESS("SUCCESS: test_BR_BTDC (A = I) " << details << " rd " << rd);
      } else if (p.tridiag_only) {
        apply(*B, *X_solve, false);
        rd = bcmm::reldif(*X, *X_solve);
        // D can be small scalar (float when scalar is double)
        // without norm termination, the error should be bounded by small scalar limit
        if (rd > 1e2*std::numeric_limits<SmallMagnitude>::epsilon())
          TEST_BR_BTDC_FAIL("FAIL: test_BR_BTDC (A = D) " << details << " rd " << rd);
        else
          TEST_BR_BTDC_SUCCESS("SUCCESS: test_BR_BTDC (A = D) " << details << " rd " << rd);
      } else if (p.tridiag_is_identity) {
        // A = I + R. In this case, T->apply(X, Y) computes
        //   Y := inv(D) (X - R Y) = I (X + R X) = A X.
        const auto& Y = X_solve;
        Y->update(-1, *X, 0);
        apply(*X, *Y, false);
        rd = bcmm::reldif(*B, *Y);
        // D can be small scalar (float when scalar is double)
        // without norm termination, the error should be bounded by small scalar limit
        if (rd > 1e2*std::numeric_limits<SmallMagnitude>::epsilon())
          TEST_BR_BTDC_FAIL("FAIL: test_BR_BTDC (A = I + R) " << details << " rd " << rd);
        else 
          TEST_BR_BTDC_SUCCESS("SUCCESS: test_BR_BTDC (A = I + R) " << details << " rd " << rd);
      } else {
        // Test that we can solve a problem.
        apply(*B, *X_solve, false);
        rd = bcmm::reldif(*X, *X_solve);
        if (rd > 1e-4)
          TEST_BR_BTDC_FAIL("FAIL: test_BR_BTDC (A = D + R) " << details << " rd " << rd);
        else
          TEST_BR_BTDC_SUCCESS("SUCCESS: test_BR_BTDC (A = D + R) " << details << " rd " << rd);
        { // Advanced options only the bare object supports.
          auto T_bare_advanced = T_bare;
          if (T_bare_advanced.is_null()) {
            T_bare_advanced = make_BTDC(sb, sbp, A, overlap_comm, nonuniform_lines, jacobi,
                                        seq_method);
            T_bare_advanced->initialize();
            T_bare_advanced->compute(T_bare_advanced->createDefaultComputeParameters());
          }
          if ( ! seq_method) { // Test norm-based termination.
            auto input = T_bare_advanced->createDefaultApplyParameters();
            input.zeroStartingSolution = zero_starting;
            input.maxNumSweeps = num_sweeps;
            input.tolerance = tol;
            const int nits = T_bare_advanced->applyInverseJacobi(*B, *X, input);
            if (nits < num_sweeps) {
              const auto n0 = T_bare_advanced->getNorms0();
              const auto nf = T_bare_advanced->getNormsFinal();

              const Magnitude r = nf/n0;
              const bool ok = r <= tol;
              if ( ! ok)
                TEST_BR_BTDC_FAIL("FAIL: test_BR_BTDC (A = D + R, norm) " << details << " r " << r);
              else
                TEST_BR_BTDC_SUCCESS("SUCCESS: test_BR_BTDC (A = D + R, norm) " << details << " r " << r);
            }
          }
          { // Test damping factor and ! zeroStartingSolution.
            const Magnitude df = 0.75;
            // Start by mimicking damping factor manually. First sweep.
            const auto X1 = Teuchos::rcp(new Tpetra_MultiVector(A->getRangeMap(), nvec));
            auto input = T_bare_advanced->createDefaultApplyParameters();
            input.zeroStartingSolution = true;
            input.maxNumSweeps = 1;
            T_bare_advanced->applyInverseJacobi(*B, *X1, input);
            X1->scale(df);
            // Second sweep. This sweep also tests ! zeroStartingSolution.
            auto X_true = Tpetra::createCopy(*X1);
            input = T_bare_advanced->createDefaultApplyParameters();
            T_bare_advanced->applyInverseJacobi(*B, *X1, input);
            X_true.update(df, *X1, 1 - df);
            // Norm-based calcs have a different code path; need to test both.
            for (const auto tol_test : {0.0, 1e-20}) {
              if (seq_method) continue;
              // Damping factor computation. 2 sweeps, starting with 0.
              input = T_bare_advanced->createDefaultApplyParameters();
              input.zeroStartingSolution = true;
              input.maxNumSweeps = 2;
              input.tolerance = tol_test;
              input.dampingFactor = df;
              T_bare_advanced->applyInverseJacobi(*B, *X, input);
              // Check if X agrees with the manually computed X_true.
              rd = bcmm::reldif(*X, X_true);
              const auto eps = 1e1*std::numeric_limits<Magnitude>::epsilon();
              if (rd > eps)
                TEST_BR_BTDC_FAIL("FAIL: test_BR_BTDC (A = D + R, damping factor) " << details << " rd " << rd << " eps " << eps);
              else
                TEST_BR_BTDC_SUCCESS("SUCCESS: test_BR_BTDC (A = D + R, damping factor) " << details << " rd " << rd << " eps " << eps);
            }
          }
        }
      }
    }
#undef TEST_BR_BTDC_FAIL
    return nerr;
  }
}; // struct BlockTriDiContainerTester

} // namespace tif_utest

#endif // HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
#endif // IFPACK2_UNITTEST_BLOCKTRIDICONTAINER_UTIL_HPP
