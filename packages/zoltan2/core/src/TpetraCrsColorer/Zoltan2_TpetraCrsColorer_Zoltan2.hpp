// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#include "Zoltan2_ColoringProblem.hpp"
#include "Zoltan2_ColoringSolution.hpp"
#include "Zoltan2_XpetraCrsMatrixAdapter.hpp"

namespace Zoltan2 {

// Implementation of CrsColorer<> using Zoltan2.  Currently this is a distance-1
// algorithm, which requires forming A^T*A, and only works in serial
template <typename CrsMatrixType> class Zoltan2CrsColorer {
public:
  using matrix_t = CrsMatrixType;
  using graph_t = typename matrix_t::crs_graph_type;
  using node_t = typename matrix_t::node_type;
  using device_t = typename node_t::device_type;
  using list_of_colors_t = Kokkos::View<int *, device_t>;
  using list_of_colors_host_t = typename list_of_colors_t::HostMirror;
  using SC = typename matrix_t::scalar_type;
  using LO = typename matrix_t::local_ordinal_type;
  using GO = typename matrix_t::global_ordinal_type;
  using vector_t = typename Tpetra::Vector<SC, LO, GO, node_t>;

  // Constructor
  Zoltan2CrsColorer(const Teuchos::RCP<matrix_t> &matrix_)
      : matrix(matrix_), graph(Impl::get_graph(matrix_)),
        colorVecLocal_(matrix_->getRowMap()),
        colorVecGlobal_(matrix_->getColMap()) {}

  // Compute coloring data
  void computeColoring(Teuchos::ParameterList &coloring_params, int &num_colors,
                       list_of_colors_host_t &list_of_colors_host,
                       list_of_colors_t &list_of_colors) {
    auto inputMat = this->matrix;

    const auto matrixType = coloring_params.get("matrixType", "Jacobian");
    const auto symmetric = coloring_params.get(
        "symmetric", (matrixType == "Jacobian" ? false : true));

    // Force symmetrize when input matrix is not symmetric (we can change that
    // once we implement Bipartate symmetrization)
    if (!symmetric) {
      // Transpose symmetrize A+A^T
      const auto nzpr = this->matrix->getGlobalMaxNumRowEntries();

      inputMat =
          Teuchos::rcp(new CrsMatrixType(matrix->getRowMap(), nzpr * nzpr));

      Tpetra::MatrixMatrix::Add(*(this->matrix), false, 1.0, *(this->matrix),
                                true, 1.0, inputMat);

      inputMat->fillComplete();
    }

    // Create Zoltan2 coloring problem and solve
    using Z2Adapter_t = Zoltan2::XpetraCrsMatrixAdapter<matrix_t>;
    Z2Adapter_t z2_adapter(inputMat);

    auto z2_params = coloring_params.sublist("Zoltan2");

    // Once we implement
    z2_params.set("color_method", "D2");

    Zoltan2::ColoringProblem<Z2Adapter_t> z2_problem(&z2_adapter, &z2_params);
    z2_problem.solve();

    // Extract colors
    auto z2_solution = z2_problem.getSolution();
    auto local_num_colors = z2_solution->getNumColors();

    auto local_list_of_colors = z2_solution->getColorsRCP();
    const auto len = local_list_of_colors.size();

    // Compute global number of colors
    auto comm = this->graph->getRowMap()->getComm();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &local_num_colors,
                       &num_colors);

    list_of_colors_host_t list_of_colors_tmp(local_list_of_colors.getRawPtr(),
                                             len);
    {
      auto colors_mv =
          colorVecLocal_.getLocalViewHost(Tpetra::Access::ReadWrite);
      auto local_colors = Kokkos::subview(colors_mv, Kokkos::ALL(), 0);

      TEUCHOS_TEST_FOR_EXCEPTION(
          local_colors.extent(0) != list_of_colors_tmp.extent(0),
          std::logic_error, "Incorrect length of color list!");

      for (size_t i = 0; i < local_colors.extent(0); ++i) {
        local_colors(i) = list_of_colors_tmp(i);
      }
    }

    using Import_t = Tpetra::Import<typename matrix_t::local_ordinal_type,
                                    typename matrix_t::global_ordinal_type,
                                    typename matrix_t::node_type>;

    // Create an importer from RowMap to ColMap
    Import_t importer(matrix->getRowMap(), matrix->getColMap());

    colorVecGlobal_.doImport(colorVecLocal_, importer, Tpetra::INSERT);

    {
      auto colors_mv =
          colorVecGlobal_.getLocalViewHost(Tpetra::Access::ReadOnly);
      auto local_colors = Kokkos::subview(colors_mv, Kokkos::ALL(), 0);
      const auto num_cols = this->matrix->getLocalNumCols();

      TEUCHOS_TEST_FOR_EXCEPTION(local_colors.extent(0) != num_cols,
                                 std::logic_error,
                                 "Incorrect length of color list!");

      list_of_colors = list_of_colors_t("list_of_colors", num_cols);
      list_of_colors_host = Kokkos::create_mirror_view(list_of_colors);

      Kokkos::deep_copy(list_of_colors_host, local_colors);
      Kokkos::deep_copy(list_of_colors, list_of_colors_host);
    }
  }

private:
  const Teuchos::RCP<matrix_t> matrix;
  const Teuchos::RCP<const graph_t> graph;
  vector_t colorVecLocal_;
  vector_t colorVecGlobal_;
};

///////////////////////////////////////////////////////////////////////////////
// Specialization of Tpetra::Zoltan2CrsColorer for BlockCrsMatrix
// Zoltan2 does not directly support BlockCrs, so this implementation
// creates a point matrix from the graph of the BlockCrs matrix
template <typename SC, typename LO, typename GO, typename NO>
class Zoltan2CrsColorer<Tpetra::BlockCrsMatrix<SC, LO, GO, NO>> {
public:
  typedef Tpetra::BlockCrsMatrix<SC, LO, GO, NO> matrix_t;
  typedef typename matrix_t::crs_graph_type graph_t;
  typedef typename matrix_t::node_type node_t;
  typedef typename node_t::device_type device_t;
  typedef Kokkos::View<int *, device_t> list_of_colors_t;
  typedef typename list_of_colors_t::HostMirror list_of_colors_host_t;

  // Constructor
  Zoltan2CrsColorer(const Teuchos::RCP<matrix_t> &matrix_)
      : matrix(matrix_), graph(Teuchos::rcp(&(matrix->getCrsGraph()), false)) {}

  // Destructor
  ~Zoltan2CrsColorer() {}

  // Compute coloring data
  void computeColoring(Teuchos::ParameterList &coloring_params, int &num_colors,
                       list_of_colors_host_t &list_of_colors_host,
                       list_of_colors_t &list_of_colors) {
    using point_matrix_t = Tpetra::CrsMatrix<SC, LO, GO, NO>;
    Teuchos::RCP<point_matrix_t> point_matrix =
        Teuchos::rcp(new point_matrix_t(graph));
    point_matrix->setAllToScalar(1.0);
    point_matrix->fillComplete();
    Zoltan2CrsColorer<point_matrix_t> point_colorer(point_matrix);
    point_colorer.computeColoring(coloring_params, num_colors,
                                  list_of_colors_host, list_of_colors);
  }

private:
  const Teuchos::RCP<matrix_t> matrix;
  const Teuchos::RCP<const graph_t> graph;
};
} // namespace Zoltan2
