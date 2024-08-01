// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_TestingUtilities.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_Core.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <Teuchos_UnitTestHarness.hpp>

namespace { // (anonymous)

using Teuchos::RCP;
using Teuchos::Comm;
using Tpetra::createNonContigMapWithNode;

template<int NumElemNodes, class LO, class GO, class NT>
class MeshInfo {
public:
  RCP<const Tpetra::Map<LO,GO,NT> > uniqueMap;
  RCP<const Tpetra::Map<LO,GO,NT> > overlapMap;
  std::vector<std::vector<GO> > element2node;

  int num_elem_shared (const GO idof, const GO jdof) const {
    int res = 0;
    for (const auto& elem : element2node) {
      if (vec_has_dof(idof) && vec_has_dof(jdof)) {
        ++res;
      }
    }
    return res;
  }

private:
  bool vec_has_dof (const std::vector<GO>& dofs, const GO dof) {
    return std::find(dofs.begin(),dofs.end(),dof)!=dofs.end();
  }
};

template<class LO, class GO, class NT>
void generate_fem2d_q1_graph(size_t numCells1D, RCP<const Comm<int> > comm , MeshInfo<4,LO,GO,NT> & mesh) {

  // We assume a NxN Q1 fem grid. Cells are partitioned among ranks.
  // Cell/Nodes ids are as follows (assuming N=2):
  //
  //  0---1---2
  //  | 0 | 1 |
  //  3---4---5
  //  | 2 | 3 |
  //  6---7---8
  //
  // The order in which we store the ids within a cell is irrelevant (so long as consistent),
  // since the local pattern is full anyways.
  //
  // NOTE: this routine sets overlapMap NOT to be the "col" map. Instead, since we partition
  //       cells, it contains all nodes belonging to local cells. If all cells are "fully" owned,
  //       then locally overlapMap=uniqueMap, while the "col" map of the graph will ALWAYS
  //       have something more (the halo).

  const size_t rank    = comm->getRank();
  const size_t numProc = comm->getSize();

  const size_t numGlobalCells = numCells1D*numCells1D;
  const size_t numNodes1D     = numCells1D+1;

  size_t numMyCells  = numGlobalCells / numProc;
  size_t remainder   = numGlobalCells % numProc;
  size_t myCellStart = numMyCells*rank;
  if (rank < remainder) {
    ++numMyCells;
    myCellStart += rank;
  } else {
    myCellStart += remainder;
  }

  // Add also repeated GIDs, we'll remove duplicates later
  Teuchos::Array<GO> ovNodes;
  mesh.element2node.resize(numMyCells);
  for (size_t cell=0; cell<numMyCells; ++cell) {
    //  0---1---2
    //  | 0 | 1 |
    //  3---4---5
    //  | 2 | 3 |
    //  6---7---8

    auto cellId = myCellStart+cell;

    auto icell = cellId / numCells1D;
    auto jcell = cellId % numCells1D;

    auto offset = icell*numNodes1D;

    // Store for the map
    ovNodes.append(offset+jcell);
    ovNodes.append(offset+jcell+1);
    ovNodes.append(offset+jcell+numNodes1D);
    ovNodes.append(offset+jcell+numNodes1D+1);

    // Store for the assembly
    mesh.element2node[cell].resize(4);
    mesh.element2node[cell][0] = offset + jcell;
    mesh.element2node[cell][1] = offset + jcell + 1;
    mesh.element2node[cell][2] = offset + jcell + numNodes1D;
    mesh.element2node[cell][3] = offset + jcell + numNodes1D + 1;
  }

  // Remove duplicates. Note: std::unique needs consecutive duplicates, so sort first
  auto start = ovNodes.data();
  auto end   = ovNodes.data()+ovNodes.size();
  std::sort(start,end);
  auto new_end = std::unique(start,end);
  auto new_size = std::distance(start, new_end);
  ovNodes.resize(new_size);

  // Create overlap map, and use Tpetra utility to create the unique one
  mesh.overlapMap = createNonContigMapWithNode<LO,GO,NT>(ovNodes(),comm);
  mesh.uniqueMap  = createOneToOne(mesh.overlapMap);
}

template<typename LO, typename GO, typename NT>
Teuchos::RCP<Tpetra::CrsGraph<LO,GO,NT>>
generate_crs_graph (const MeshInfo<4,LO,GO,NT>& mesh)
{
  using CG = Tpetra::CrsGraph<LO,GO,NT>;

  Teuchos::RCP<CG> g(new CG(mesh.uniqueMap,9));
  for (const auto& elem_dofs : mesh.element2node) {
    for (const GO gid_i : elem_dofs) {
      for (const GO gid_j : elem_dofs) {
        g->insertGlobalIndices(gid_i,1,&gid_j);
      }
    }
  }
  g->fillComplete();

  return g;
}

template<typename LO, typename GO, typename NT>
Teuchos::RCP<Tpetra::FECrsGraph<LO,GO,NT>>
generate_fecrs_graph (const MeshInfo<4,LO,GO,NT>& mesh)
{
  using FEG = Tpetra::FECrsGraph<LO,GO,NT>;

  Teuchos::RCP<FEG> feg(new FEG(mesh.uniqueMap,mesh.overlapMap,9,mesh.overlapMap));
  feg->beginAssembly();
  for (const auto& elem_dofs : mesh.element2node) {
    for (const GO gid_i : elem_dofs) {
      for (const GO gid_j : elem_dofs) {
        feg->insertGlobalIndices(gid_i,1,&gid_j);
      }
    }
  }
  feg->endAssembly();

  return feg;
}

// Builds a FECrs matrix with pattern generated from the above routine.
// Loops on cells to fill the matrix. On each cell, add -1 to off-diagonal
// entries, and N to diag entries, where N is the number of local nodes.
// Note that such matrix is strictly diagonally dominant.
template<typename ST, typename LO, typename GO, typename NT>
void
fill_matrices (Tpetra::FECrsMatrix<ST,LO,GO,NT>& fe_mat,
               Tpetra::CrsMatrix<ST,LO,GO,NT>& mat,
               const MeshInfo<4,LO,GO,NT>& mesh)
{
  const ST zero = Teuchos::ScalarTraits<ST>::zero();

  fe_mat.beginAssembly();
  mat.resumeFill();

  fe_mat.setAllToScalar(zero);
  mat.setAllToScalar(zero);

  Teuchos::Array<GO> col(1);
  Teuchos::Array<ST> val(1);
  for (const auto& elem_dofs : mesh.element2node) {
    for (const GO& idof : elem_dofs) {
      for (const GO& jdof : elem_dofs) {
        col[0] = jdof;
        val[0] = Teuchos::ScalarTraits<ST>::random();

        fe_mat.sumIntoGlobalValues (idof, col(), val());
           mat.sumIntoGlobalValues (idof, col(), val());
      }
    }
  }
  fe_mat.endAssembly();
  mat.fillComplete();
}

template<typename ST, typename LO, typename GO, typename NT>
bool compare_matrices (const Tpetra::CrsMatrix<ST,LO,GO,NT>& A,
                       const Tpetra::CrsMatrix<ST,LO,GO,NT>& B,
                       Teuchos::FancyOStream &out)
{
  using TST = Teuchos::ScalarTraits<ST>;
  using MT = typename TST::magnitudeType;
  auto eps = 1000*Teuchos::ScalarTraits<MT>::eps();

  // They should have the same row/range/domain maps
  if (!A.getRowMap()->isSameAs(*B.getRowMap())) {
    out<<"Compare: RowMap failed.\n";
    return false;
  }
  if (!A.getRangeMap()->isSameAs(*B.getRangeMap())) {
    out<<"Compare: RangeMap failed.\n";
    return false;
  }
  if (!A.getDomainMap()->isSameAs(*B.getDomainMap())) {
    out<<"Compare: DomainMap failed.\n";
    return false;
  }

  // Now we can test the equality of the rows, in a mathematical way.
  // We do not care about the order in which entries appear, only the "mathematical" object.
  const auto& gA = *A.getGraph();
  const auto& gB = *B.getGraph();
  const LO num_my_rows = gA.getLocalNumRows();
  if (num_my_rows!=static_cast<LO>(gB.getLocalNumRows())) {
    out << "Compare: number of local rows differ on some MPI rank: "
        << num_my_rows << " vs " << gB.getLocalNumRows() << ".\n";
    return false;
  }

  typedef typename Tpetra::CrsMatrix<ST,LO,GO,NT> crs_matrix_type;
  auto findLID = [](
       const typename crs_matrix_type::local_inds_host_view_type& lids,
       const LO lid) -> int {
    auto it = std::find(lids.data(),lids.data()+lids.extent(0),lid);
    if (it==lids.data()+lids.extent(0)) {
      return -1;
    } else {
      return std::distance(lids.data(),it);
    }
  };

  typename crs_matrix_type::values_host_view_type  valsA, valsB;
  typename crs_matrix_type::local_inds_host_view_type colsA, colsB;
  const LO invLO = Teuchos::OrdinalTraits<LO>::invalid();
  const auto& colMapA = *gA.getColMap();
  const auto& colMapB = *gB.getColMap();

  Kokkos::fence(); // must protect UVM access

  for (LO irow=0; irow<num_my_rows; ++irow) {
    const GO grow = gA.getRowMap()->getGlobalElement(irow);

    // Extract rows
    A.getLocalRowView(irow, colsA, valsA);
    B.getLocalRowView(irow, colsB, valsB);

    // If different row sizes, then the two matrices are different.
    auto numEntries = colsA.size();
    if (numEntries!=colsB.size()) {
      out << "Compare: global row " << grow << " has different lengths.\n";
      return false;
    }

    // Loop over rows entries
    for (size_t j=0; j<numEntries; ++j) {
      const LO lidA = colsA[j];
      const GO gid = colMapA.getGlobalElement(lidA);
      const LO lidB = colMapB.getLocalElement(gid);

      // If B does not have this GID in its column map, then the two matrices are different.
      if (lidB==invLO) {
        out << "Compare: col maps store different global indices.\n";
        return false;
      }

      // If B does not have this GID in this row, then the two matrices are different.
      int pos = findLID(colsB,lidB);
      if (pos==-1) {
        out << "Compare: global row " << grow << " has different global indicess.\n";
        return false;
      }

      // Finally, check the numerical values

      if( TST::magnitude(valsB[pos]-valsA[j]) > eps) {
        out << "Compare: global row " << grow << " has different values.\n";
        return false;
      }
    }
  }

  // If we got this far, then none of the negative checks happened,
  // so the matrices are (locally) truly identical (from the mathematical point of view)
  return true;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL (Tpetra_MatMat, FECrsMatrix, SC, LO, GO, NT)
{
  using FEMAT = typename Tpetra::FECrsMatrix<SC,LO,GO,NT>;
  using   MAT = typename Tpetra::CrsMatrix<SC,LO,GO,NT>;

  // get a comm
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  // Generate a mesh
  const int numCells1D = 4;
  MeshInfo<4,LO,GO,NT> mesh;
  generate_fem2d_q1_graph(numCells1D,comm,mesh);

  auto fe_graph = generate_fecrs_graph(mesh);
  auto    graph = generate_crs_graph(mesh);

  FEMAT feA(fe_graph);
    MAT   A(   graph);
  fill_matrices(feA,A,mesh);
  success = compare_matrices(A,feA,out);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success);

  FEMAT feB(fe_graph);
    MAT   B(   graph);
  fill_matrices(feB,B,mesh);
  success = compare_matrices(B,feB,out);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success);

  for (bool transA : {false, true}) {
    for (bool transB : {false, true}) {
      Teuchos::RCP<Teuchos::ParameterList> params1(new Teuchos::ParameterList());
      Teuchos::RCP<Teuchos::ParameterList> params2(new Teuchos::ParameterList());
      params2->set("MM_TAFC_OptimizationCoreCount",1);
      for (auto params : {params1, params2}) {

        // A and feA should have the same row map, so pick one.
        auto C_row_map = transA ? feA.getDomainMap() : feA.getRangeMap();

        // For the test, use a ridicolously large upper bound for the nnz per row
        MAT feC(C_row_map,feA.getGraph()->getGlobalNumEntries());
        MAT   C(C_row_map,  A.getGraph()->getGlobalNumEntries());

        Tpetra::MatrixMatrix::Multiply(feA, transA, feB, transB, feC, true, "", params);
        Tpetra::MatrixMatrix::Multiply(  A, transA,   B, transB,   C, true, "", params);

        success = compare_matrices(C,feC,out);

        TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success);
      }
    }
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )			\
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, FECrsMatrix, SC, LO, GO, NT)

  TPETRA_ETI_MANGLING_TYPEDEFS()

// FIXME_SYCL
#if !defined(KOKKOS_ENABLE_SYCL) || 1
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )
#endif

} // anonymous namespace
