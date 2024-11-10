// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <utility>
#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_FECrsGraph.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Assembly_Helpers.hpp"
#include "Tpetra_Details_getNumDiags.hpp"


namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Tpetra::createContigMapWithNode;
using Tpetra::createNonContigMapWithNode;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::outArg;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using std::endl;
typedef Tpetra::global_size_t GST;

template<class RP, class CI>
void print_graph(int rank, const char *prefix, RP rowptr, CI colind) {
  printf("[%d] %s entries = ",rank,prefix);
  for(size_t i=0; i<rowptr.extent(0)-1; i++) {
    printf("( ");
    for(size_t j=rowptr[i]; j<rowptr[i+1]; j++)
      printf("%2d ",colind[j]);
    printf(") ");
    printf("\n");
  }
}

template<class LO, class GO, class Node>
bool compare_final_graph_structure(Teuchos::FancyOStream &out,Tpetra::CrsGraph<LO,GO,Node> & g1, Tpetra::CrsGraph<LO,GO,Node> & g2) {
  using std::endl;
  if (!g1.isFillComplete() || !g2.isFillComplete()) {out<<"Compare: FillComplete failed"<<endl;return false;}
  if (!g1.getRangeMap()->isSameAs(*g2.getRangeMap())) {out<<"Compare: RangeMap failed"<<endl;return false;}
  if (!g1.getRowMap()->isSameAs(*g2.getRowMap())) {out<<"Compare: RowMap failed"<<endl;return false;}
  if (!g1.getDomainMap()->isSameAs(*g2.getDomainMap())) {out<<"Compare: DomainMap failed"<<endl;return false;}
  if (!g1.getColMap()->isSameAs(*g2.getColMap())) {out<<"Compare: ColMap failed"<<endl;g1.describe(out,Teuchos::VERB_EXTREME);g2.describe(out,Teuchos::VERB_EXTREME);return false;}

  auto rowptr1 = g1.getLocalGraphHost().row_map;
  auto rowptr2 = g2.getLocalGraphHost().row_map;

  auto colind1 = g1.getLocalGraphHost().entries;
  auto colind2 = g2.getLocalGraphHost().entries;

  if (rowptr1.extent(0) != rowptr2.extent(0)) {out<<"Compare: rowptr extent failed"<<endl;return false;}
  if (colind1.extent(0) != colind2.extent(0)) {out<<"Compare: colind extent failed: "<<colind1.extent(0)<<" vs "<<colind2.extent(0)<<endl;
    int rank = g1.getRowMap()->getComm()->getRank();
    print_graph(rank,"G1",rowptr1,colind1);
    print_graph(rank,"G2",rowptr2,colind2);
    return false;
}

  bool success=true;
  TEST_COMPARE_ARRAYS(rowptr1,rowptr2);
  if (!success) {out<<"Compare: rowptr match failed"<<endl;return false;}

  TEST_COMPARE_ARRAYS(colind1,colind2);
  if (!success) {out<<"Compare: colind match failed"<<endl;return false;}


  return true;
}

// This routine checks that two graphs are the same "up to permutation of indices inside rows".
// This means that the two graphs must have the same row/range/domain maps, the same number
// of rows, and the same GLOBAL indices in each row. The order in which the indices appear
// in the row, as well as the local index they are given in the column map is irrelevant.
template<class LO, class GO, class Node>
bool compare_final_graph_structure_relaxed(Teuchos::FancyOStream &out,
                                        const Tpetra::CrsGraph<LO,GO,Node> & g1,
                                        const Tpetra::CrsGraph<LO,GO,Node> & g2) {
  // Make sure we finished filling the two graphs
  if (!g1.isFillComplete() || !g2.isFillComplete()) {
    out << "Compare: FillComplete failed.\n";
    return false;
  }

  // Range/domain/row maps *must* be the same
  if (!g1.getRangeMap()->isSameAs(*g2.getRangeMap())) {
    out<<"Compare: RangeMap failed.\n";
    return false;
  }
  if (!g1.getRowMap()->isSameAs(*g2.getRowMap())) {
    out << "Compare: RowMap failed.\n";
    return false;
  }
  if (!g1.getDomainMap()->isSameAs(*g2.getDomainMap())) {
    out << "Compare: DomainMap failed.\n";
    return false;
  }

  const LO num_my_rows = g1.getLocalNumRows();
  if (num_my_rows!=static_cast<LO>(g2.getLocalNumRows())) {
    out << "Compare: number of local rows differ on some MPI rank: "
        << num_my_rows << " vs " << g2.getLocalNumRows() << ".\n";
    return false;
  }

  typedef typename Tpetra::CrsGraph<LO,GO,Node>::local_inds_host_view_type
                   lcl_ind_type;

  auto hasLID = [](const lcl_ind_type & lids, const LO lid) -> bool {
    auto it = std::find(lids.data(),lids.data()+lids.extent(0),lid);
    return it!=lids.data()+lids.extent(0);
  };

  lcl_ind_type cols1, cols2;
  const LO invLO = Teuchos::OrdinalTraits<LO>::invalid();
  const auto& colMap1 = *g1.getColMap();
  const auto& colMap2 = *g2.getColMap();
  for (LO irow=0; irow<num_my_rows; ++irow) {
    const GO grow = g1.getRowMap()->getGlobalElement(irow);

    // Extract rows
    g1.getLocalRowView(irow, cols1);
    g2.getLocalRowView(irow, cols2);

    // If different row sizes, then the two graphs are different.
    auto numEntries = cols1.size();
    if (numEntries!=cols2.size()) {
      out << "Compare: global row " << grow << " has different lengths.\n";
      return false;
    }

    // Loop over rows indices
    for (decltype(numEntries) j=0; j<numEntries; ++j) {
      const LO lid1 = cols1[j];
      const GO gid = colMap1.getGlobalElement(lid1);
      const LO lid2 = colMap2.getLocalElement(gid);

      // If g2 does not have this GID in its column map, then the two graphs are different.
      if (lid2==invLO) {
        out << "Compare: col maps store different global indices.\n";
        return false;
      }

      // If g2 does not have this GID in this row, then the two graphs are different.
      if (!hasLID(cols2,lid2)) {
        out << "Compare: global row " << grow << " has different global indicess.\n";
        return false;
      }
    }
  }

  // If we got this far, then none of the negative checks happened,
  // so the graphs are (locally) truly identical (from the mathematical point of view)
  return true;
}


template<int NumElemNodes, class LO, class GO, class Node>
class GraphPack {
public:
  RCP<const Tpetra::Map<LO,GO,Node> > uniqueMap;
  RCP<const Tpetra::Map<LO,GO,Node> > overlapMap;
  std::vector<std::vector<GO> > element2node;

  typedef Kokkos::View<LO*[NumElemNodes], Kokkos::LayoutLeft, typename Node::device_type > k_element2node_type;
  k_element2node_type k_element2node;

  void print(int rank, std::ostream & out) {
    using std::endl;
    out << "["<<rank<<"] Unique Map  : ";
    for(size_t i=0; i<uniqueMap->getLocalNumElements(); i++)
      out << uniqueMap->getGlobalElement(i) << " ";
    out<<endl;

    out << "["<<rank<<"] Overlap Map : ";
    for(size_t i=0; i<overlapMap->getLocalNumElements(); i++)
      out << overlapMap->getGlobalElement(i) << " ";
    out<<endl;

    out << "["<<rank<<"] element2node: ";
    for(size_t i=0; i<(size_t)element2node.size(); i++) {
      out <<"(";
      for(size_t j=0; j<(size_t)element2node[i].size(); j++)
        out<<element2node[i][j] << " ";
      out<<") ";
    }
    out<<endl;
  }
};


template<class LO, class GO, class Node>
void generate_fem1d_graph(size_t numLocalNodes, RCP<const Comm<int> > comm , GraphPack<2,LO,GO,Node> & pack) {
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
  int rank    = comm->getRank();
  int numProc = comm->getSize();
  size_t numOverlapNodes = numLocalNodes; if(rank!=numProc-1) numOverlapNodes++;  if(rank!=0) numOverlapNodes++;
  size_t numLocalElements = (rank == numProc-1) ? numLocalNodes -1 : numLocalNodes;
  //  printf("CMS numOverlapNodes = %d numLocalElements = %d\n",numOverlapNodes,numLocalElements);

  pack.uniqueMap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocalNodes,comm);

  Teuchos::Array<GO> overlapIndices(numOverlapNodes);
  for(size_t i=0; i<numLocalNodes; i++) {
    overlapIndices[i] = pack.uniqueMap->getGlobalElement(i);
  }
  size_t last = numLocalNodes;
  if(rank != 0)           {overlapIndices[last] = overlapIndices[0] - 1; last++;}
  if(rank != numProc -1)  {overlapIndices[last] = overlapIndices[numLocalNodes-1] + 1; last++;}

  pack.overlapMap = rcp(new Tpetra::Map<LO,GO,Node>(INVALID,overlapIndices,0,comm));

  pack.element2node.resize(numLocalElements);
  for(size_t i=0; i<numLocalElements; i++) {
    pack.element2node[i].resize(2);
    pack.element2node[i][0] = pack.uniqueMap->getGlobalElement(i);
    pack.element2node[i][1] = pack.uniqueMap->getGlobalElement(i) + 1;
  }

  // Kokkos version of the element2node array
  Kokkos::resize(pack.k_element2node,numLocalElements);
  auto k_e2n = pack.k_element2node;
  auto l_umap = pack.uniqueMap->getLocalMap();
  Kokkos::parallel_for(Kokkos::RangePolicy<typename Node::execution_space>(0,numLocalElements),KOKKOS_LAMBDA(const size_t i) {
      GO gid = l_umap.getGlobalElement(i);
      k_e2n(i,0) = gid;
      k_e2n(i,1) = gid+1;
    });
  Kokkos::fence();
}

template<class LO, class GO, class Node>
void generate_fem2d_q1_graph(size_t numCells1D, RCP<const Comm<int> > comm , GraphPack<4,LO,GO,Node> & pack) {

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
  pack.element2node.resize(numMyCells);
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
    pack.element2node[cell].resize(4);
    pack.element2node[cell][0] = offset + jcell;
    pack.element2node[cell][1] = offset + jcell + 1;
    pack.element2node[cell][2] = offset + jcell + numNodes1D;
    pack.element2node[cell][3] = offset + jcell + numNodes1D + 1;
  }

  // Remove duplicates. Note: std::unique needs consecutive duplicates, so sort first
  auto start = ovNodes.data();
  auto end   = ovNodes.data()+ovNodes.size();
  std::sort(start,end);
  auto new_end = std::unique(start,end);
  auto new_size = std::distance(start, new_end);
  ovNodes.resize(new_size);

  // Create overlap map, and use Tpetra utility to create the unique one
  pack.overlapMap = createNonContigMapWithNode<LO,GO,Node>(ovNodes(),comm);
  pack.uniqueMap  = createOneToOne(pack.overlapMap);

  // Kokkos version of the element2node array
  Kokkos::resize(pack.k_element2node,numMyCells);
  auto k_e2n = pack.k_element2node;
  Kokkos::parallel_for(Kokkos::RangePolicy<typename Node::execution_space>(0,numMyCells),
                       KOKKOS_LAMBDA(const size_t cell) {
      auto cellId = myCellStart+cell;

      auto icell = cellId / numCells1D;
      auto jcell = cellId % numCells1D;

      auto offset = icell*numNodes1D;

      k_e2n(icell,0) = offset + jcell;
      k_e2n(icell,1) = offset + jcell + 1;
      k_e2n(icell,2) = offset + jcell + numNodes1D;
      k_e2n(icell,3) = offset + jcell + numNodes1D + 1;
    });
  Kokkos::fence();
}

//
// UNIT TESTS
//






////
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FECrsGraph, Diagonal, LO, GO, Node )
{
    typedef Tpetra::FECrsGraph<LO,GO,Node> FEG;
    typedef Tpetra::CrsGraph<LO,GO,Node> CG;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    // create a Map
    const size_t numLocal = 10;
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);


    // Trivial test that makes sure a diagonal graph can be built
    CG g1(map,1);
    FEG g2(map,map,1);

    Tpetra::beginAssembly(g2);
    for(size_t i=0; i<numLocal; i++) {
      GO gid = map->getGlobalElement(i);
      g1.insertGlobalIndices(gid,1,&gid);
      g2.insertGlobalIndices(gid,1,&gid);
    }
    Tpetra::endAssembly(g2);
    g1.fillComplete();

    success = compare_final_graph_structure(out,g1,g2) && success;
    TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FECrsGraph, Diagonal_LocalIndex, LO, GO, Node )
{
    typedef Tpetra::FECrsGraph<LO,GO,Node> FEG;
    typedef Tpetra::CrsGraph<LO,GO,Node> CG;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    // create a Map
    const size_t numLocal = 10;
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);


    // Trivial test that makes sure a diagonal graph can be built
    CG g1(map,1);
    FEG g2(map,map,map,1);

    Tpetra::beginAssembly(g2);
    for(size_t i=0; i<numLocal; i++) {
      GO gid = map->getGlobalElement(i);
      LO lid = (LO) i;
      g1.insertGlobalIndices(gid,1,&gid);
      g2.insertLocalIndices(lid,1,&lid);
    }
    Tpetra::endAssembly(g2);
    g1.fillComplete();

    success = compare_final_graph_structure(out,g1,g2);
    TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}


////
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FECrsGraph, Assemble1D, LO, GO, Node )
{
  typedef Tpetra::FECrsGraph<LO,GO,Node> FEG;
  typedef Tpetra::CrsGraph<LO,GO,Node> CG;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();

  // create a Map
  const size_t numLocal = 10;

  // Generate a mesh
  GraphPack<2,LO,GO,Node> pack;
  generate_fem1d_graph(numLocal,comm,pack);
  //pack.print(comm->getRank(),std::cout);

  // Comparative assembly
  // FIXME: We should be able to get away with 3 here, but we need 4 since duplicates are
  // not being handled correctly.
  CG g1(pack.uniqueMap,4);
  FEG g2(pack.uniqueMap,pack.overlapMap,4);

  g2.beginAssembly();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        //        printf("Inserting (%d,%d)\n",gid_j,gid_k);
        g1.insertGlobalIndices(gid_j,1,&gid_k);
        g2.insertGlobalIndices(gid_j,1,&gid_k);
      }
    }
  }
  g1.fillComplete();
  g2.endAssembly();

  success = compare_final_graph_structure(out,g1,g2);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}


////
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FECrsGraph, Assemble1D_LocalIndex, LO, GO, Node )
{
  typedef Tpetra::FECrsGraph<LO,GO,Node> FEG;
  typedef Tpetra::CrsGraph<LO,GO,Node> CG;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();

  // create a Map
  const size_t numLocal = 10;

  // Generate a mesh
  GraphPack<2,LO,GO,Node> pack;
  generate_fem1d_graph(numLocal,comm,pack);
  //  pack.print(comm->getRank(),std::cout);fflush(stdout);

  // Comparative assembly
  // FIXME: We should be able to get away with 3 here, but we need 4 since duplicates are
  // not being handled correctly.

  CG g1(pack.uniqueMap,4);
  FEG g2(pack.uniqueMap,pack.overlapMap,pack.overlapMap,4);

  g2.beginAssembly();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      LO lid_j = pack.overlapMap->getLocalElement(gid_j);
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        LO lid_k = pack.overlapMap->getLocalElement(gid_k);
        //        printf("[%d] Inserting gid (%d,%d) lid (%d,%d)\n",comm->getRank(),gid_j,gid_k,lid_j,lid_k);fflush(stdout);
        g1.insertGlobalIndices(gid_j,1,&gid_k);
        g2.insertLocalIndices(lid_j,1,&lid_k);
      }
    }
  }
  g1.fillComplete();
  g2.endAssembly();

  success = compare_final_graph_structure(out,g1,g2);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FECrsGraph, Assemble2D_OPSDomain, LO, GO, Node )
{
  typedef Tpetra::FECrsGraph<LO,GO,Node> FEG;
  typedef Tpetra::CrsGraph<LO,GO,Node> CG;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();

  const int numCells1D = 4;

  // Generate a mesh
  GraphPack<4,LO,GO,Node> pack;
  generate_fem2d_q1_graph(numCells1D,comm,pack);

  // Comparative assembly

  CG  g0(pack.overlapMap,9);
  CG  g1(pack.uniqueMap,9);
  FEG g2(pack.uniqueMap,pack.overlapMap,9);
  FEG g3(pack.uniqueMap,pack.overlapMap,9,pack.overlapMap);

  g2.beginAssembly();
  g3.beginAssembly();
  for (const auto& cell_dofs : pack.element2node) {
    for (const auto& gid_i : cell_dofs) {
      for (const auto& gid_j : cell_dofs) {
        g0.insertGlobalIndices(gid_i,1,&gid_j);
        g1.insertGlobalIndices(gid_i,1,&gid_j);
        g2.insertGlobalIndices(gid_i,1,&gid_j);
        g3.insertGlobalIndices(gid_i,1,&gid_j);
      }
    }
  }
  g0.fillComplete();
  g1.fillComplete();
  g2.endAssembly();
  g3.endAssembly();
  Tpetra::Import<LO,GO,Node> import(pack.uniqueMap,pack.overlapMap);

  Teuchos::FancyOStream myout(Teuchos::rcpFromRef(std::cout));
  CG  g11(pack.uniqueMap,9);
  g11.doExport(g0,import,Tpetra::INSERT);
  g11.fillComplete(pack.uniqueMap,pack.uniqueMap);
  success = compare_final_graph_structure(out,g11,g1);
  TPETRA_GLOBAL_SUCCESS_CHECK(myout,comm,success)

  success = compare_final_graph_structure(out,g1,g2) && success;
  TPETRA_GLOBAL_SUCCESS_CHECK(myout,comm,success)

  // Passing a domain map different from the uniqueMap should cause a rearrangement of
  // off-rank indices in the column map. Other than that, the graphs should still coincide.
  success = compare_final_graph_structure_relaxed(out,g1,g3) && success;
  TPETRA_GLOBAL_SUCCESS_CHECK(myout,comm,success)

  // We can check that g3 is exactly what we would get if g1's col map had been
  // built with a overlapMap locally fitted to it. We can recreate that scenario by
  // a) use overlapMap as domainMap during fillComplete (this will generate the
  // correct colMap), and b) reset domainMap to be the unique map (so that the
  // checks in compare_final_graph_structure do not fail b/c of the domain map).
  CG  g12(pack.uniqueMap,9);
  g12.doExport(g0,import,Tpetra::INSERT);
  g12.fillComplete(pack.overlapMap,pack.uniqueMap);
  auto importer = g12.getImporter();
  g12.replaceDomainMapAndImporter(pack.uniqueMap,importer);

  success = compare_final_graph_structure(myout,g12,g3) && success;
  TPETRA_GLOBAL_SUCCESS_CHECK(myout,comm,success)
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FECrsGraph, Diagonal, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FECrsGraph, Diagonal_LocalIndex, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FECrsGraph, Assemble1D, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FECrsGraph, Assemble1D_LocalIndex, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FECrsGraph, Assemble2D_OPSDomain, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // end namespace (anonymous)
