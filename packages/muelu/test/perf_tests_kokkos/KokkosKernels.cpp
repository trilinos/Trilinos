#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <impl/Kokkos_Timer.hpp>

typedef Kokkos::DefaultExecutionSpace                   execution_space;
typedef typename execution_space::memory_space          memory_space;
typedef Kokkos::Device<execution_space, memory_space>   device_type;

typedef double                                          scalar_type;
typedef int                                             local_ordinal_type;

typedef Kokkos::ArithTraits<scalar_type>                ATS;
typedef typename ATS::mag_type                          magnitude_type;

typedef Kokkos::StaticCrsGraph<local_ordinal_type,
    Kokkos::LayoutLeft, device_type>                    local_graph_type;
typedef typename local_graph_type::row_map_type         row_map_type;
typedef typename local_graph_type::entries_type         entries_type;

typedef Kokkos::View<const local_ordinal_type*,
        device_type>                                    row_type;


typedef Kokkos::CrsMatrix<scalar_type, local_ordinal_type,
    execution_space, void,
    typename local_graph_type::size_type>               local_matrix_type;
typedef Kokkos::View<bool*, device_type>                boundary_nodes_type;

typedef typename local_matrix_type::size_type           size_type;

class LWGraph_kokkos {
public:
  LWGraph_kokkos(const local_graph_type& graph) : graph_(graph) {
    minLocalIndex_ = 0;
    maxLocalIndex_ = graph_.numRows();
  }
  ~LWGraph_kokkos() { }

  //! Return number of graph vertices
  KOKKOS_INLINE_FUNCTION size_type GetNodeNumVertices() const {
    return graph_.numRows();
  }
  //! Return number of graph edges
  KOKKOS_INLINE_FUNCTION size_type GetNodeNumEdges() const {
    return graph_.row_map(GetNodeNumVertices());
  }

  //! Return the list of vertices adjacent to the vertex 'v'.
  KOKKOS_INLINE_FUNCTION row_type getNeighborVertices(local_ordinal_type i) const;

  //! Return true if vertex with local id 'v' is on current process.
  KOKKOS_INLINE_FUNCTION bool isLocalNeighborVertex(local_ordinal_type i) const {
    return i >= minLocalIndex_ && i <= maxLocalIndex_;
  }

  //! Set boolean array indicating which rows correspond to Dirichlet boundaries.
  KOKKOS_INLINE_FUNCTION void SetBoundaryNodeMap(const boundary_nodes_type bndry) {
    dirichletBoundaries_ = bndry;
  }

  //! Returns the maximum number of entries across all rows/columns on this node
  KOKKOS_INLINE_FUNCTION size_type getNodeMaxNumRowEntries () const {
    return maxNumRowEntries_;
  }

  //! Returns map with global ids of boundary nodes.
  KOKKOS_INLINE_FUNCTION const boundary_nodes_type GetBoundaryNodeMap() const {
    return dirichletBoundaries_;
  }

private:

  //! Underlying graph (with label)
  const local_graph_type      graph_;

  //! Boolean array marking Dirichlet rows.
  boundary_nodes_type         dirichletBoundaries_;

  //! Local index boundaries (cached from domain map)
  local_ordinal_type    minLocalIndex_, maxLocalIndex_;
  size_type             maxNumRowEntries_;
};

local_matrix_type kernel_construct(local_ordinal_type numRows) {
  local_ordinal_type numCols         = numRows;
  local_ordinal_type nnz             = 10*numRows;

  local_ordinal_type varianz_nel_row = 0.2*nnz/numRows;
  local_ordinal_type width_row       = 0.01*numRows;

  local_ordinal_type elements_per_row = nnz/numRows;

  local_ordinal_type *rowPtr = new local_ordinal_type[numRows+1];
  rowPtr[0] = 0;
  for (int row = 0; row < numRows; row++) {
    int varianz = (1.0*rand()/INT_MAX-0.5)*varianz_nel_row;

    rowPtr[row+1] = rowPtr[row] + elements_per_row + varianz;
  }
  nnz = rowPtr[numRows];

  local_ordinal_type *colInd = new local_ordinal_type[nnz];
  scalar_type        *values = new scalar_type       [nnz];
  for (int row = 0; row < numRows; row++) {
    for (int k = rowPtr[row]; k < rowPtr[row+1]; k++) {
      int pos = row + (1.0*rand()/INT_MAX-0.5)*width_row;

      if (pos <  0)       pos += numCols;
      if (pos >= numCols) pos -= numCols;
      colInd[k] = pos;
      values[k] = 100.0*rand()/INT_MAX - 50.0;
    }
  }

  return local_matrix_type("A", numRows, numCols, nnz, values, rowPtr, colInd, false/*pad*/);
}

void kernel_coalesce_drop(local_matrix_type A) {
  auto numRows = A.numRows();

  scalar_type eps = 0.05;

  // Stage 0: detect Dirichlet rows
  boundary_nodes_type boundaryNodes("boundaryNodes", numRows);
  Kokkos::parallel_for("MueLu:Utils::DetectDirichletRows", numRows, KOKKOS_LAMBDA(const local_ordinal_type row) {
    auto rowView = A.row (row);
    auto length  = rowView.length;

    boundaryNodes(row) = true;
    for (decltype(length) colID = 0; colID < length; colID++)
      if ((rowView.colidx(colID) != row) && (ATS::magnitude(rowView.value(colID)) > 1e-13)) {
        boundaryNodes(row) = false;
        break;
      }
  });

  // Stage 1: calculate the number of remaining entries per row
  typename entries_type::non_const_type diag("ghosted", numRows);
  typename row_map_type::non_const_type rows("row_map", numRows+1);       // rows(0) = 0 automatically


  local_ordinal_type realnnz = 0;
  Kokkos::parallel_reduce("kernel_cd:stage1_reduce", numRows,
    KOKKOS_LAMBDA(const local_ordinal_type row, local_ordinal_type& nnz) {
      auto rowView = A.row (row);
      auto length  = rowView.length;

      local_ordinal_type rownnz = 0;
      for (decltype(length) colID = 0; colID < length; colID++) {
        local_ordinal_type col = rowView.colidx(colID);

        // Avoid square root by using squared values
        magnitude_type aiiajj = eps*eps * ATS::magnitude(diag(row))*ATS::magnitude(diag(col));                  // eps^2*|a_ii|*|a_jj|
        magnitude_type aij2   = ATS::magnitude(rowView.value(colID)) * ATS::magnitude(rowView.value(colID));    // |a_ij|^2

        if (aij2 > aiiajj || row == col)
          rownnz++;
      }
      rows(row+1) = rownnz;
      nnz += rownnz;
    }, realnnz);

  // parallel_scan (exclusive)
  Kokkos::parallel_scan("kernel_cd:stage1_scan", numRows+1,
    KOKKOS_LAMBDA(const local_ordinal_type i, local_ordinal_type& upd, const bool& final_result) {
      upd += rows(i);
      if (final_result)
        rows(i) = upd;
    });

  // Stage 2: fill in the column indices
  typename boundary_nodes_type::non_const_type bndNodes("boundaryNodes", numRows);
  typename entries_type::non_const_type        cols    ("entries",       realnnz);
  local_ordinal_type numDropped = 0;
  Kokkos::parallel_reduce("kernel_cd:stage2_reduce", numRows,
    KOKKOS_LAMBDA(const local_ordinal_type row, local_ordinal_type& dropped) {
      auto rowView = A.row (row);
      auto length = rowView.length;

      local_ordinal_type rownnz = 0;
      for (decltype(length) colID = 0; colID < length; colID++) {
        local_ordinal_type col = rowView.colidx(colID);

        // Avoid square root by using squared values
        magnitude_type aiiajj = eps*eps * ATS::magnitude(diag(row))*ATS::magnitude(diag(col));                  // eps^2*|a_ii|*|a_jj|
        magnitude_type aij2   = ATS::magnitude(rowView.value(colID)) * ATS::magnitude(rowView.value(colID));    // |a_ij|^2

        if (aij2 > aiiajj || row == col) {
          cols(rows(row) + rownnz) = col;
          rownnz++;
        } else {
          dropped++;
        }
        if (rownnz == 1) {
          // If the only element remaining after filtering is diagonal, mark node as boundary
          // FIXME: this should really be replaced by the following
          //    if (indices.size() == 1 && indices[0] == row)
          //        boundaryNodes[row] = true;
          // We do not do it this way now because there is no framework for distinguishing isolated
          // and boundary nodes in the aggregation algorithms
          bndNodes(row) = true;
        }
      }
    }, numDropped);
  boundaryNodes = bndNodes;

  local_graph_type kokkosGraph(cols, rows);

  LWGraph_kokkos graph(kokkosGraph);
  graph.SetBoundaryNodeMap(boundaryNodes);
}

int main(int argc, char **argv) {
  int n    = 100000;
  int loop = 10;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) { n    = atoi(argv[++i]); continue; }
    if (strcmp(argv[i], "-l") == 0) { loop = atoi(argv[++i]); continue; }
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      std::cout << "./MueLu_KokkosKernels.exe [-n <matrix_size>] [-l <number_of_loops>]" << std::endl;
      return 0;
    }
  }

  Kokkos::initialize(argc, argv);

  srand(13721);

  local_matrix_type A = kernel_construct(n);

  execution_space::fence();
  Kokkos::Impl::Timer timer;

  for (int i = 0; i < loop; i++)
    kernel_coalesce_drop(A);

  double kernel_time = timer.seconds();

  execution_space::fence();

  printf("kernel_coalesce_drop: %.2e (s)\n", kernel_time / loop);

  execution_space::finalize();
}
