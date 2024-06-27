// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestBlockCrsUtil.hpp

\brief Ifpack2 Unit and performance test utilities for preconditioners that use
       Tpetra::BlockCrsMatrix.
*/

#ifndef IFPACK2_UNITTEST_BLOCKCRS_UTIL_HPP
#define IFPACK2_UNITTEST_BLOCKCRS_UTIL_HPP

#include <Tpetra_Map.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockMultiVector.hpp>

namespace tif_utest {

template <typename Scalar, typename LO, typename GO>
struct BlockCrsMatrixMaker {
  typedef LO Int;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef Tpetra::Map<LO, GO> Tpetra_Map;
  typedef typename Tpetra_Map::node_type::device_type DeviceType;
  typedef Tpetra::Import<LO, GO> Tpetra_Import;
  typedef Tpetra::Export<LO, GO> Tpetra_Export;
  typedef Tpetra::MultiVector<Scalar, LO, GO> Tpetra_MultiVector;
  typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO>
          Tpetra_MultiVector_Magnitude;
  typedef Tpetra::BlockMultiVector<Scalar, LO, GO> Tpetra_BlockMultiVector;
  typedef Tpetra::CrsGraph<LO, GO> Tpetra_CrsGraph;
  typedef Tpetra::RowMatrix<Scalar, LO, GO> Tpetra_RowMatrix;
  typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO> Tpetra_BlockCrsMatrix;

  // Representation of a structured block mesh. The fastest index is k.
  struct StructuredBlock {
    const Int ni, nj, nk;
    StructuredBlock (const Int ni_, const Int nj_, const Int nk_, const bool contiguous = true)
      : ni(ni_), nj(nj_), nk(nk_), njnk_(nj_*nk_), contiguous_(contiguous) {}
    KOKKOS_INLINE_FUNCTION GO size () const { return ni*njnk_; }
    KOKKOS_INLINE_FUNCTION bool is_contiguous () const { return contiguous_; }
    KOKKOS_INLINE_FUNCTION GO ijk2id (const Int i, const Int j, const Int k) const {
      return (keep_same_gids_ || contiguous_ ?
              (i*nj + j)*nk + k :
              (k*nj + j)*ni + i);
    }
    KOKKOS_INLINE_FUNCTION void id2ijk (const GO id, Int& i, Int& j, Int& k) const {
      if (keep_same_gids_ || contiguous_) {
        i = id / njnk_;
        k = id % njnk_;
        j = k / nk;
        k = k % nk;
      } else {
        const auto njni = nj*ni;
        k = id / njni;
        i = id % njni;
        j = i / ni;
        i = i % ni;
      }
    }
  private:
    const GO njnk_;
    const bool contiguous_;
    // Originally, I tested with different GIDs in the noncontiguous case. But
    // that means the math problem ends up being different. I prefer the
    // contiguous and noncontiguous tests to test the same linear equation.
    enum : bool { keep_same_gids_ = true };
  };

  // Part of the block owned by this process. Split the block in i and j directions.
  struct StructuredBlockPart {
    Int is, ie, js, je, ks, ke;
    StructuredBlockPart (const Int iis, const Int iie, const Int ijs, const Int ije,
                         const Int iks, const Int ike)
      : is(iis), ie(iie), js(ijs), je(ije), ks(iks), ke(ike)
    {}
  };

  struct StencilShape { enum Enum { cross }; };

  static void partition_n_uniformly (const Int n, const Int nparts, std::vector<Int>& p) {
    p.resize(nparts + 1);
    const Int base = n / nparts;
    Int rem = n - base*nparts;
    Int extra = rem > 0 ? 1 : 0;
    p[0] = 0;
    for (Int i = 1; i <= nparts; ++i) {
      p[i] = p[i-1] + base + extra;
      if (rem > 0) {
        --rem;
        if (rem == 0) extra = 0;
      }
    }
  }

  static StructuredBlockPart make_StructuredBlockPart (
    const StructuredBlock& sb, const Int isplit, const Int jsplit, const Int rank)
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(isplit > sb.ni || jsplit > sb.nj,
                                "splits must be <= dimension; isplit jsplit ni nj " <<
                                isplit << " " << jsplit << " " << sb.ni << " " << sb.nj);
    const Int ios = rank / jsplit;
    const Int jos = rank % jsplit;
    std::vector<Int> p;
    p.reserve(std::max(isplit, jsplit) + 1);
    partition_n_uniformly(sb.ni, isplit, p);
    const Int is = p[ios];
    const Int ie = p[ios+1];
    TEUCHOS_ASSERT(is >= 0 && is < ie && ie <= sb.ni);
    partition_n_uniformly(sb.nj, jsplit, p);
    const Int js = p[jos];
    const Int je = p[jos+1];
    TEUCHOS_ASSERT(js >= 0 && js < je && je <= sb.nj);
    return StructuredBlockPart(is, ie, js, je, 0, sb.nk);
  }

  // Correctness check the testing utility struct StructuredBlockPart.
  static Int test_StructuredBlockPart (const Int ni_max, const Int nj_max, const bool contiguous) {
    Int nerr = 0;
    const Int nk = 10;
    std::vector<Int> hit;
    hit.reserve(ni_max*nj_max);
    for (Int ni = 1; ni <= ni_max; ++ni) {
      for (Int nj = 1; nj <= nj_max; ++nj) {
        hit.resize(ni*nj*nk, 0);
        StructuredBlock sb(ni, nj, nk, contiguous);
        for (Int nrank = 0, max_nrank = ni*nj; nrank < max_nrank; ++nrank) {
          for (Int isplit = 1, i_max = std::min(nrank, ni); isplit <= i_max; ++isplit) {
            const Int jsplit = nrank / isplit;
            if (jsplit >= nj || isplit*jsplit != nrank) continue;
            std::fill(hit.begin(), hit.end(), 0);
            for (Int rank = 0; rank < nrank; ++rank) {
              auto sbp = make_StructuredBlockPart(sb, isplit, jsplit, rank);
              for (Int i = sbp.is; i < sbp.ie; ++i)
                for (Int j = sbp.js; j < sbp.je; ++j)
                  for (Int k = sbp.ks; k < sbp.ke; ++k)
                    ++hit[sb.ijk2id(i, j, k)];
            }
            bool all1 = true;
            for (const auto e : hit)
              if (e != 1) {
                all1 = false;
                break;
              }
            if ( ! all1) ++nerr;
          }
        }
      }
    }
    return nerr;
  }

  // blockrow(b,i,j) is the set of blocks in this block row. (I,J,K) is the cell
  // in the mesh. n_blocks is the number of blocks in this row. diag_idx is the
  // index into 0 : n_blocks-1 corresponding to the diagonal block. max_blocks is
  // the max number of blocks in a row or column.
  //   Fill the values with reasonably interesting, deterministically generated
  // (not random), parallel consistent, values.
  template <typename T> static void
  make_entry (const Int& B1, const Int& B2, const Int& i, const Int& j, const GO& c,
              T& e) {
    e = ((c+i-j) % 14)/12.0 + ((B1-i) % 13)/14.0 - ((B2+j) % 17)/12.0 + 1e-7;
  }
#ifdef HAVE_TEUCHOS_COMPLEX
  template <typename T> static void
  make_entry (const Int& B1, const Int& B2, const Int& i, const Int& j, const GO& c,
              std::complex<T>& e) {
    T re, im;
    make_entry(B1, B2, i, j, c  , re);
    make_entry(B2, B1, j, i, c+1, im);
    e = std::complex<T>(re, im);
  }
#endif

  template <typename Array3D> static void fill_block_row (
    const Int I, const Int J, const Int K, const GO* gids, const Int n_blocks,
    const Int diag_idx, const Int max_blocks, Array3D blockrow)
  {
    const Int bs = blockrow.extent_int(2);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(blockrow.extent_int(1) != bs, "Blocks must be square.");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(blockrow.extent_int(0) < n_blocks,
                                "blockrow is not consistent with n_blocks: blockrow.extent(0) = "
                                << blockrow.extent(0) << " but n_blocks = " << n_blocks);
    const Int B1 = I + 3*J + 5*K;
    const Int B2 = 5*I + J + 3*K;
    for (Int b = 0; b < n_blocks; ++b)
      for (Int i = 0; i < bs; ++i)
        for (Int j = 0; j < bs; ++j) {
          const auto c = gids[b];
          make_entry(B1, B2, i, j, c, blockrow(b,i,j));
        }
  }

  template <typename T> static KOKKOS_INLINE_FUNCTION T abs (const T& v) { return std::abs(v); }
  template <typename T> static KOKKOS_INLINE_FUNCTION T abs (const Kokkos::complex<T>& v) { return Kokkos::abs(v); }
  template <typename T> static KOKKOS_INLINE_FUNCTION T signof (const T& v) { return v >= 0 ? 1 : -1; }
  template <typename T> static KOKKOS_INLINE_FUNCTION Kokkos::complex<T> signof (const Kokkos::complex<T>& v)
  { return v == 0 ? 1 : v / Kokkos::abs(v); }

  // Given a matrix A with a 1D (row) partition, find new diagonal values so that
  // it is row and column diagonally dominant. This requires global communication
  // of local column abs sums.
  static void make_row_and_col_diag_dominant (Tpetra_BlockCrsMatrix& a) {
    const auto& g = a.getCrsGraph();
    const auto& rowptr = g.getLocalGraphHost().row_map;
    const auto& colidx = g.getLocalGraphHost().entries;
    const auto& values = a.getValuesHostNonConst();

    const auto row_map = g.getRowMap();
    const auto col_map = g.getColMap();
    const LO bs = a.getBlockSize(), bs2 = bs*bs;
    const LO nrows = rowptr.size() - 1;

    // Accumulate into these arrays.
    std::vector<Magnitude> rowsum(bs*a.getLocalNumRows(), 0);
    const auto cpm = Teuchos::rcp(
      new Tpetra_Map(Tpetra_BlockMultiVector::makePointMap(*col_map, bs)));
    Tpetra_MultiVector_Magnitude colsum_mv(cpm, 1);
    colsum_mv.putScalar(0);
    {
      auto colsum = colsum_mv.getLocalViewHost(Tpetra::Access::ReadWrite);

      // Get off-diag 1-norms.
      Kokkos::fence(); // uvm access
      for (LO r = 0; r < nrows; ++r) {
        const auto rgid = row_map->getGlobalElement(r);
        for (size_t j = rowptr(r); j < rowptr(r+1); ++j) {
          const LO c = colidx(j);
          const auto cgid = col_map->getGlobalElement(c);
          const bool diag_block = cgid == rgid;
          auto* const block = &values(j*bs2);
          for (Int bi = 0; bi < bs; ++bi)
            for (Int bj = 0; bj < bs; ++bj) {
              if (diag_block && bj == bi) continue;
              const auto e = abs(block[bi*bs + bj]);
              rowsum[bs*r + bi] += e;
              colsum(bs*c + bj, 0) += e;
            }
        }
      }
    }

    { // Is this the best way to do overlap -ADD-> overlap?
      // overlap -ADD-> nonoverlap
      Tpetra_MultiVector_Magnitude d(a.getDomainMap(), 1); {
        Tpetra_Export exporter(cpm, a.getDomainMap());
        d.doExport(colsum_mv, exporter, Tpetra::ADD);
      }
      // nonoverlap -REPLACE-> overlap
      Tpetra_Import importer(a.getDomainMap(), cpm);
      colsum_mv.doImport(d, importer, Tpetra::REPLACE);
    }

    {
      auto colsum = colsum_mv.getLocalViewHost(Tpetra::Access::ReadOnly);
      // Modify diag entries.
      for (LO r = 0; r < nrows; ++r) {
        const auto rgid = row_map->getGlobalElement(r);
        for (size_t j = rowptr(r); j < rowptr(r+1); ++j) {
          const LO c = colidx(j);
          const auto cgid = col_map->getGlobalElement(c);
          const bool diag_block = cgid == rgid;
          if ( ! diag_block) continue;
          auto* const block = &values(j*bs2);
          for (Int bi = 0; bi < bs; ++bi) {
            auto& e = block[bi*bs + bi];
            e = Magnitude(1.01)*std::max(rowsum[bs*r + bi], colsum(bs*c + bi, 0))*signof(e);
          }
        }
      }
    }
  }

  static Teuchos::RCP<Tpetra_CrsGraph>
  make_crs_graph (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                  const StructuredBlock& sb, const StructuredBlockPart& sbp,
                  const bool tridiags_only = false, const bool different_maps = false) {
    const GO idx_base = 0;
    const GO num_gbl_cell = sb.ni*sb.nj*sb.nk;

    Teuchos::RCP<Tpetra_Map> row_map;
    std::vector<GO> my_row_gids;
    {
      const GO num_my_cell = (sbp.ie - sbp.is)*(sbp.je - sbp.js)*(sbp.ke - sbp.ks);
      my_row_gids.resize(num_my_cell, 0);
      {
        GO cnt = 0;
        for (Int i = sbp.is; i < sbp.ie; ++i) {
          if (sb.is_contiguous()) {
            for (Int j = sbp.js; j < sbp.je; ++j)
              for (Int k = sbp.ks; k < sbp.ke; ++k, ++cnt)
                my_row_gids[cnt] = sb.ijk2id(i,j,k);
          } else {
            // The noncontiguity test modifies multiple things: NC1. Make tridiags
            // pull noncontiguously from the LID space.
            for (Int k = sbp.ks; k < sbp.ke; ++k)
              for (Int j = sbp.js; j < sbp.je; ++j, ++cnt)
                my_row_gids[cnt] = sb.ijk2id(i,j,k);
          }
        }
        row_map = Teuchos::rcp(new Tpetra_Map(num_gbl_cell, my_row_gids.data(),
                                              my_row_gids.size(), idx_base, comm));
      }
    }

    Teuchos::RCP<Tpetra_Map> col_map;
    {
      const auto pad_extent = [] (Int& s, Int& e, const Int ncell) {
        Int n = e - s;
        if (s > 0) { ++n; --s; }
        if (e < ncell) { ++n; ++e; }
        return n;
      };
      StructuredBlockPart sbppad(sbp);
      pad_extent(sbppad.is, sbppad.ie, sb.ni);
      pad_extent(sbppad.js, sbppad.je, sb.nj);
      pad_extent(sbppad.ks, sbppad.ke, sb.nk);
      std::vector<GO> my_col_gids;
      {
        my_col_gids.insert(my_col_gids.begin(), my_row_gids.begin(), my_row_gids.end());
        if ( ! sb.is_contiguous()) {
          // NC2. Do a different ordering of (i,j,k) than in the row case. Thus,
          // lclrow != lclcol.
          std::vector<GO> v(my_col_gids);
          size_t j = 0;
          for (size_t i = 0; i < v.size(); i += 2) my_col_gids[j++] = v[i];
          for (size_t i = 1; i < v.size(); i += 2) my_col_gids[j++] = v[i];
        }
        const Int lims[6][3][2] =
          {{{sbppad.is, sbp.is}, {sbp.js, sbp.je}, {sbp.ks, sbp.ke}},
           {{sbp.ie, sbppad.ie}, {sbp.js, sbp.je}, {sbp.ks, sbp.ke}},
           {{sbp.is, sbp.ie}, {sbppad.js, sbp.js}, {sbp.ks, sbp.ke}},
           {{sbp.is, sbp.ie}, {sbp.je, sbppad.je}, {sbp.ks, sbp.ke}},
           {{sbp.is, sbp.ie}, {sbp.js, sbp.je}, {sbppad.ks, sbp.ks}},
           {{sbp.is, sbp.ie}, {sbp.js, sbp.je}, {sbp.ke, sbppad.ke}}};
        for (Int li = 0; li < 6; ++li)
          for (Int i = lims[li][0][0]; i < lims[li][0][1]; ++i)
            for (Int j = lims[li][1][0]; j < lims[li][1][1]; ++j)
              for (Int k = lims[li][2][0]; k < lims[li][2][1]; ++k)
                my_col_gids.push_back(sb.ijk2id(i,j,k));
        col_map = Teuchos::rcp(new Tpetra_Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                              my_col_gids.data(), my_col_gids.size(), idx_base, comm));
      }
    }

    typename Tpetra_CrsGraph::local_graph_device_type g;
    {
      typedef typename Tpetra_CrsGraph::local_graph_device_type::row_map_type row_map_type;
      typedef typename Tpetra_CrsGraph::local_graph_device_type::entries_type entries_type;
      const GO nr = my_row_gids.size();
      typename row_map_type::non_const_type::HostMirror rowptr("rowptr", nr + 1);
      typename entries_type::HostMirror colidx;
      GO nnz = 0;
      for (Int pass = 0; pass < 2; ++pass) {
        if (pass == 0) {
          // Counting pass.
          Kokkos::deep_copy(rowptr, 0);
        } else {
          // Fill pass.
          colidx = typename entries_type::HostMirror("colidx", nnz);
          // Cumsum.
          for (Int i = 2; i <= nr; ++i)
            rowptr(i) += rowptr(i-1);
        }
        // For each owned cell:
        for (Int i = sbp.is; i < sbp.ie; ++i)
          for (Int j = sbp.js; j < sbp.je; ++j)
            for (Int k = sbp.ks; k < sbp.ke; ++k) {
              const GO my_gid = sb.ijk2id(i,j,k);
              const LO row_lid = row_map->getLocalElement(my_gid);
              // For each nbr and myself:
              for (Int ni = std::max(i-1, 0); ni < std::min(i+2, sb.ni); ++ni)
                for (Int nj = std::max(j-1, 0); nj < std::min(j+2, sb.nj); ++nj) {
                  if (tridiags_only && (ni != i || nj != j)) continue;
                  for (Int nk = std::max(k-1, 0); nk < std::min(k+2, sb.nk); ++nk) {
                    if (Int(ni == i) + Int(nj == j) + Int(nk == k) < 2) continue;
                    if (pass == 0) {
                      if (row_lid+2 <= nr)
                        ++rowptr(row_lid+2);
                      ++nnz;
                    } else {
                      const GO nbr_gid = sb.ijk2id(ni,nj,nk);
                      const LO col_lid = col_map->getLocalElement(nbr_gid);
                      colidx(rowptr(row_lid+1)++) = col_lid;
                    }
                  }
                }
            }
      }
      TEUCHOS_ASSERT(static_cast<GO>(rowptr(nr)) == nnz);

      // Sort columns in each row.
      for (Int r = 0; r < nr; ++r)
        std::sort(colidx.data() + rowptr(r), colidx.data() + rowptr(r+1));

      {
        typename row_map_type::non_const_type row_map_tmp("row_map", rowptr.size());
        Kokkos::deep_copy(row_map_tmp, rowptr);
        entries_type entries("entries", colidx.size());
        Kokkos::deep_copy(entries, colidx);
        g = typename Tpetra_CrsGraph::local_graph_device_type(entries, row_map_tmp);
      }

      if ( ! tridiags_only) {
        const GO n_full_bpr = std::max(0, (sb.nk - 2)*(sb.nj - 2)*(sb.ni - 2));
        const GO n_4_bpr = (sb.nk >= 2 ? 2 : 1)*(sb.nj >= 2 ? 2 : 1)*(sb.ni >= 2 ? 2 : 1);
        Int i_full_bpr = 0, i_4_bpr = 0;
        for (Int r = 0; r < nr; ++r) {
          const Int bpr = rowptr(r+1) - rowptr(r);
          if (bpr == 7) ++i_full_bpr;
          else if (bpr == 4) ++i_4_bpr;
        }
        GO lcl[] = {i_full_bpr, i_4_bpr}, gbl[2];
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 2, lcl, gbl);
        TEUCHOS_ASSERT(gbl[0] == n_full_bpr);
        TEUCHOS_ASSERT(gbl[1] == n_4_bpr);
      }
    }

    if (different_maps) {
      auto G = Teuchos::rcp(new Tpetra_CrsGraph(row_map, col_map, g.row_map, g.entries));
      Teuchos::RCP<Tpetra_Map> rng_map, dmn_map;
      { // Make all the maps have different orderings. These additional
        // permutations are not decomposition-independent, unlike every other
        // aspect of these manufactured problems.
        std::vector<GO> my_gids(my_row_gids);
        const auto rank = comm->getRank();
        // Deterministic permutation.
        for (size_t i = 0; i < my_gids.size(); ++i)
          std::swap(my_gids[i],
                    my_gids[((rank + 3*i) + (i % 11)*(i % 7)) % my_gids.size()]);
        rng_map = Teuchos::rcp(new Tpetra_Map(num_gbl_cell, my_gids.data(),
                                              my_gids.size(), idx_base, comm));
        for (size_t i = 0; i < my_gids.size(); ++i)
          std::swap(my_gids[i],
                    my_gids[((rank + 2*i) + (i % 13)*(i % 5)) % my_gids.size()]);
        dmn_map = Teuchos::rcp(new Tpetra_Map(num_gbl_cell, my_gids.data(),
                                              my_gids.size(), idx_base, comm));
      }
      G->fillComplete(dmn_map, rng_map);
      return G;
    } else {
      return Teuchos::rcp(new Tpetra_CrsGraph(row_map, col_map, g, Teuchos::null));
    }
  }

  static void
  get_offdiag_idxs (const StructuredBlock& sb, const Tpetra_CrsGraph& g, const Tpetra_Map& col_map,
                    const Int& lr, const Int& I, const Int& J, const Int& K, Int offdiag_idxs[2]) {
    offdiag_idxs[0] = offdiag_idxs[1] = -1;
    const auto& rowptr = g.getLocalGraphHost().row_map;
    const auto& colidx = g.getLocalGraphHost().entries;
    GO rid_offdiags[2];
    rid_offdiags[0] = rid_offdiags[1] = Teuchos::OrdinalTraits<GO>::invalid();
    if (K > 0) rid_offdiags[0] = sb.ijk2id(I, J, K-1);
    if (K+1 < sb.nk) rid_offdiags[1] = sb.ijk2id(I, J, K+1);
    for (size_t j = rowptr(lr); j < rowptr(lr+1); ++j) {
      const auto cid = col_map.getGlobalElement(colidx(j));
      for (Int k = 0; k < 2; ++k)
        if (cid == rid_offdiags[k]) {
          offdiag_idxs[k] = j - rowptr(lr);
          break;
        }
    }
    if (K > 0) TEUCHOS_ASSERT(offdiag_idxs[0] != -1);
    if (K+1 < sb.nk) TEUCHOS_ASSERT(offdiag_idxs[1] != -1);
  }

  template <typename View> static void
  zero_offdiag_idxs (const Int offdiag_idxs[2], View& blockrow) {
    const Int bs = blockrow.extent(2);
    for (Int k = 0; k < 2; ++k) {
      if (offdiag_idxs[k] == -1) continue;
      for (Int i = 0; i < bs; ++i)
        for (Int j = 0; j < bs; ++j)
          blockrow(offdiag_idxs[k],i,j) = 0;
    }
  }

  // Make a BlockCrsMatrix with interesting, deterministic entries independent of
  // parallel decomposition. If ! tridiag_is_identity, then the matrix is row and
  // column diagonally dominant.
  static Teuchos::RCP<Tpetra_BlockCrsMatrix>
  make_bcrs_matrix (const StructuredBlock& sb, const Teuchos::RCP<const Tpetra_CrsGraph>& gr,
                    const Int bs, const bool tridiag_is_identity = false, const bool block_diag = false) {
    auto mr = Teuchos::rcp(new Tpetra_BlockCrsMatrix(*gr, bs));
    // Raw pointers for threading.
    auto m = mr.get();
    auto g = gr.get();
    const auto& rowptr = g->getLocalGraphHost().row_map;
    const auto& colidx = g->getLocalGraphHost().entries;
    const LO nr = rowptr.extent_int(0) - 1;
    const auto row_map = g->getRowMap().get();
    const auto col_map = g->getColMap().get();
    Int max_bpr = 0; {
      Int max_bpr_lcl = 0;
      for (Int r = 0; r < nr; ++r)
        max_bpr_lcl = std::max<Int>(max_bpr_lcl, rowptr(r+1) - rowptr(r));
      Teuchos::reduceAll(*row_map->getComm(), Teuchos::REDUCE_MAX, 1, &max_bpr_lcl, &max_bpr);
    }
    const Int nthreads =
#ifdef KOKKOS_ENABLE_OPENMP
      omp_get_max_threads()
#else
      1
#endif
      ;
    struct ThreadData {
      Kokkos::View<typename Kokkos::ArithTraits<Scalar>::val_type***, Kokkos::HostSpace> blockrow;
      std::vector<GO> cids;
      ThreadData (const Int max_bpr, const Int bs)
        : blockrow("blockrow", max_bpr, bs, bs), cids(max_bpr)
      {}
    };
    std::vector<ThreadData> tds;
    for (Int tid = 0; tid < nthreads; ++tid)
      tds.push_back(ThreadData(max_bpr, bs));
#ifdef KOKKOS_ENABLE_OPENMP
#   pragma omp parallel for
#endif
    for (Int lr = 0; lr < nr; ++lr) {
      const Int tid =
#ifdef KOKKOS_ENABLE_OPENMP
        omp_get_thread_num()
#else
        0
#endif
        ;
      auto& blockrow = tds[tid].blockrow;
      auto& cids = tds[tid].cids;
      const auto rid = row_map->getGlobalElement(lr);
      Int I, J, K;
      sb.id2ijk(rid, I, J, K);
      const auto n_blocks = rowptr(lr+1) - rowptr(lr);
      Int diag_idx = -1;
      const size_t j0 = rowptr(lr);
      for (size_t j = j0; j < rowptr(lr+1); ++j) {
        const auto cid = col_map->getGlobalElement(colidx(j));
        const auto os = j - j0;
        if (cid == rid) diag_idx = os;
        cids[os] = cid;
      }
      TEUCHOS_ASSERT(diag_idx != -1);
      Int offdiag_idxs[2];
      if (tridiag_is_identity || block_diag)
        get_offdiag_idxs(sb, *g, *col_map, lr, I, J, K, offdiag_idxs);
      else
        offdiag_idxs[0] = offdiag_idxs[1] = -1;
      fill_block_row(I, J, K, cids.data(), n_blocks, diag_idx, max_bpr, blockrow);
      if (tridiag_is_identity)
        for (Int i = 0; i < bs; ++i)
          for (Int j = 0; j < bs; ++j)
            blockrow(diag_idx,i,j) = i == j ? 1 : 0;
      if (tridiag_is_identity || block_diag)
        zero_offdiag_idxs(offdiag_idxs, blockrow);
      for (size_t j = rowptr(lr); j < rowptr(lr+1); ++j) {
        auto block = m->getLocalBlockHostNonConst(lr, colidx(j));
        const auto b = j - rowptr(lr);
        for (Int bi = 0; bi < bs; ++bi)
          for (Int bj = 0; bj < bs; ++bj)
            block(bi,bj) = blockrow(b,bi,bj);
      }
    }
    if ( ! (tridiag_is_identity || (block_diag && bs == 1)))
      make_row_and_col_diag_dominant(*m);
    return mr;
  }

  // Make a multivector with interesting, deterministic entries independent of
  // parallel decomposition.
  template <typename T> static void make_entry (const GO& gid, const LO& col, T& e) {
    e = std::sin(0.1*(1 + gid + col));
  }
#ifdef HAVE_TEUCHOS_COMPLEX
  template <typename T> static void make_entry (const GO& gid, const LO& col, std::complex<T>& e) {
    T re, im;
    make_entry(gid, col  , re);
    make_entry(gid, col+1, im);
    e = std::complex<T>(re, im);
  }
  template <typename T> static void make_entry (const GO& gid, const LO& col, Kokkos::complex<T>& e) {
    std::complex<T> se;
    make_entry(gid, col, se);
    e = se;
  }
#endif

  static Teuchos::RCP<Tpetra_MultiVector> make_multivector (
    const StructuredBlock& sb, const Teuchos::RCP<const Tpetra_BlockCrsMatrix>& m,
    const Int bs, const Int nvec)
  {
    auto mv = Teuchos::rcp(new Tpetra_MultiVector(m->getDomainMap(), nvec));
    auto v = mv->getLocalViewHost(Tpetra::Access::OverwriteAll);
    const auto map = mv->getMap();
    for (GO lid = 0; lid < v.extent_int(0); ++lid)
      for (LO col = 0; col < v.extent_int(1); ++col) {
        const auto gid = map->getGlobalElement(lid);
        make_entry(gid, col, v(lid,col));
      }
    return mv;
  }

  static Magnitude reldif (const Tpetra_MultiVector& x, const Tpetra_MultiVector& y) {
    const Int nvec = x.getNumVectors();
    Tpetra_MultiVector xmy(x, Teuchos::Copy);
    xmy.update(-1, y, 1);
    Teuchos::Array<Magnitude> num_norm2(nvec), den_norm2(nvec);
    xmy.norm2(num_norm2());
    x.norm2(den_norm2());
    Magnitude num = 0, den = 0;
    for (Int i = 0; i < nvec; ++i) {
      num += num_norm2[i]*num_norm2[i];
      den += den_norm2[i]*den_norm2[i];
    }
    const Magnitude rd = std::sqrt(num/den);
    return (std::isnan(rd) || std::isinf(rd)) ? 1 : rd;
  }

  // Check the manufactured test matrix for parallel consistency by comparing a
  // serial MVP with a parallel one.
  static Int
  test_bcrs_matrix (const Teuchos::RCP<const Teuchos::Comm<int> >& pcomm,
                    const StructuredBlock& sb, const StructuredBlockPart& psbp,
                    const Int bs, const bool tridiags_only = false,
                    const bool different_maps = false) {
    // Construct matrix and perform MVP. Return the range map if requested.
    auto mvp = [&] (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const StructuredBlockPart& sbp,
                    const bool want_rng_map, Teuchos::RCP<const Tpetra_Map>& rng_map) {
      auto g = make_crs_graph(comm, sb, sbp, tridiags_only, different_maps);
      auto m = make_bcrs_matrix(sb, g, bs);
      if (want_rng_map) rng_map = m->getRangeMap();
      const Int nvec = 3;
      auto x = make_multivector(sb, m, bs, nvec);
      auto y = Teuchos::rcp(new Tpetra_MultiVector(m->getRangeMap(), nvec));
      m->apply(*x, *y);
      return y;
    };
    // Parallel MVP.
    Teuchos::RCP<const Tpetra_Map> srl_rng_map;
    auto py = mvp(pcomm, psbp, false, srl_rng_map);
    // Serial MVP, and get the range map.
    Teuchos::RCP<Tpetra_MultiVector> sy;
    if (pcomm->getRank() == 0) {
      const auto scomm = Teuchos::rcp(new Teuchos::SerialComm<int>());
      const auto ssbp = make_StructuredBlockPart(sb, 1, 1, 0);
      sy = mvp(scomm, ssbp, true, srl_rng_map);
    }
    // Bring the parallel MVP's solution to one processor, using the serial range
    // map's global index list to set up the target.
    Teuchos::RCP<Tpetra_MultiVector> y0; {
      const GO idx_base = 0;
      const GO num_gbl_cell = sb.ni*sb.nj*sb.nk, N = bs*num_gbl_cell;
      std::vector<GO> my_gids;
      if (pcomm->getRank() == 0) {
        // Only the root rank gets any GIDs, of course.
        my_gids.resize(N);
        const auto& sr_gids = srl_rng_map->getMyGlobalIndices();
        for (GO i = 0; i < N; ++i) my_gids[i] = sr_gids[i];
      }
      const auto allon0 = Teuchos::rcp(new Tpetra_Map(N, my_gids.data(), my_gids.size(),
                                                      idx_base, pcomm));
      y0 = Teuchos::rcp(new Tpetra_MultiVector(allon0, py->getNumVectors()));
      Tpetra_Export e(py->getMap(), allon0);
      y0->doExport(*py, e, Tpetra::REPLACE);
    }
    // Check the difference.
    Magnitude relerr = 0; {
      if (pcomm->getRank() == 0)
        relerr = reldif(*sy, *y0);
      Magnitude grelerr = 0;
      Teuchos::reduceAll(*pcomm, Teuchos::REDUCE_MAX, 1, &relerr, &grelerr);
      relerr = grelerr;
    }
    return relerr > 1e1*std::numeric_limits<Magnitude>::epsilon();
  }
}; // struct BlockCrsMatrixMaker

} // namespace tif_utest

#endif // IFPACK2_UNITTEST_BLOCKCRS_UTIL_HPP
