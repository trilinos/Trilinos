// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CLASSICALMAPFACTORY_DEF_HPP_
#define MUELU_CLASSICALMAPFACTORY_DEF_HPP_

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#endif

#include <Xpetra_Vector.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_ClassicalMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_LWGraph.hpp"

#ifdef HAVE_MUELU_ZOLTAN2
#include "MueLu_Zoltan2GraphAdapter.hpp"
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_ColoringSolution.hpp>

#endif

#include "MueLu_LWGraph_kokkos.hpp"
#include <KokkosGraph_Distance1ColorHandle.hpp>
#include <KokkosGraph_Distance1Color.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: deterministic");
  SET_VALID_ENTRY("aggregation: coloring algorithm");
  SET_VALID_ENTRY("aggregation: coloring: use color graph");
#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase> >("Graph", null, "Generating factory of the graph");
  validParamList->set<RCP<const FactoryBase> >("Coloring Graph", null, "Generating factory of the graph");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "UnAmalgamationInfo");
  Input(currentLevel, "Graph");

  const ParameterList& pL = GetParameterList();
  bool use_color_graph    = pL.get<bool>("aggregation: coloring: use color graph");
  if (use_color_graph)
    Input(currentLevel, "Coloring Graph");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);
  const ParameterList& pL = GetParameterList();
  RCP<const Matrix> A     = Get<RCP<Matrix> >(currentLevel, "A");

  RCP<const LWGraph> graph;
  bool use_color_graph = pL.get<bool>("aggregation: coloring: use color graph");
  if (use_color_graph)
    graph = Get<RCP<LWGraph> >(currentLevel, "Coloring Graph");
  else
    graph = Get<RCP<LWGraph> >(currentLevel, "Graph");

  /* ============================================================= */
  /* Phase 1 : Compute an initial MIS                              */
  /* ============================================================= */
  ArrayRCP<LO> myColors;
  LO numColors = 0;

  RCP<LocalOrdinalVector> fc_splitting;
  std::string coloringAlgo = pL.get<std::string>("aggregation: coloring algorithm");

  // Switch to Zoltan2 if we're parallel and Tpetra (and not file)
#ifdef HAVE_MUELU_ZOLTAN2
  int numProcs = A->getRowMap()->getComm()->getSize();
  if (coloringAlgo != "file" && numProcs > 1 && graph->GetDomainMap()->lib() == Xpetra::UseTpetra)
    coloringAlgo = "Zoltan2";
#endif

    //#define CMS_DUMP
#ifdef CMS_DUMP
  {
    int rank = graph->GetDomainMap()->getComm()->getRank();

    printf("[%d,%d] graph local size = %dx%d\n", rank, currentLevel.GetLevelID(), (int)graph->GetDomainMap()->getLocalNumElements(), (int)graph->GetImportMap()->getLocalNumElements());

    std::ofstream ofs(std::string("m_dropped_graph_") + std::to_string(currentLevel.GetLevelID()) + std::string("_") + std::to_string(rank) + std::string(".dat"), std::ofstream::out);
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(ofs));
    graph->print(*fancy, Debug);
  }
  {
    A->getRowMap()->getComm()->barrier();
  }

#endif

  // Switch to MIS if we're in Epetra (and not file)
  if (coloringAlgo != "file" && graph->GetDomainMap()->lib() == Xpetra::UseEpetra)
    coloringAlgo = "MIS";

  if (coloringAlgo == "file") {
    // Read the CF splitting from disk
    // NOTE: For interoperability reasons, this is dependent on the point_type enum not changing
    std::string map_file   = std::string("map_fcsplitting_") + std::to_string(currentLevel.GetLevelID()) + std::string(".m");
    std::string color_file = std::string("fcsplitting_") + std::to_string(currentLevel.GetLevelID()) + std::string(".m");

    FILE* mapfile               = fopen(map_file.c_str(), "r");
    using real_type             = typename Teuchos::ScalarTraits<SC>::magnitudeType;
    using RealValuedMultiVector = typename Xpetra::MultiVector<real_type, LO, GO, NO>;
    RCP<RealValuedMultiVector> mv;

    GetOStream(Statistics1) << "Reading FC splitting from " << color_file << ", using map file " << map_file << ". On rank " << A->getRowMap()->getComm()->getRank() << " local size is " << A->getRowMap()->getLocalNumElements() << std::endl;
    if (mapfile) {
      fclose(mapfile);
      RCP<const Map> colorMap = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMap(map_file, A->getRowMap()->lib(), A->getRowMap()->getComm());
      TEUCHOS_TEST_FOR_EXCEPTION(!colorMap->isCompatible(*A->getRowMap()), std::invalid_argument, "Coloring on disk has incompatible map with A");

      mv = Xpetra::IO<real_type, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(color_file, colorMap);
    } else {
      // Use A's rowmap and hope it matches
      mv = Xpetra::IO<real_type, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(color_file, A->getRowMap());
    }
    TEUCHOS_TEST_FOR_EXCEPTION(mv.is_null(), std::invalid_argument, "Coloring on disk cannot be read");
    fc_splitting = LocalOrdinalVectorFactory::Build(A->getRowMap());
    TEUCHOS_TEST_FOR_EXCEPTION(mv->getLocalLength() != fc_splitting->getLocalLength(), std::invalid_argument, "Coloring map mismatch");

    // Overlay the Dirichlet Points (and copy out the rest)
    auto boundaryNodes                = graph->GetBoundaryNodeMap();
    ArrayRCP<const real_type> mv_data = mv->getData(0);
    ArrayRCP<LO> fc_data              = fc_splitting->getDataNonConst(0);
    for (LO i = 0; i < (LO)fc_data.size(); i++) {
      if (boundaryNodes[i])
        fc_data[i] = DIRICHLET_PT;
      else
        fc_data[i] = Teuchos::as<LO>(mv_data[i]);
    }
  }
#ifdef HAVE_MUELU_ZOLTAN2
  else if (coloringAlgo.find("Zoltan2") != std::string::npos && graph->GetDomainMap()->lib() == Xpetra::UseTpetra) {
    SubFactoryMonitor sfm(*this, "DistributedGraphColoring", currentLevel);
    DoDistributedGraphColoring(graph, myColors, numColors);
  }
#endif
  else if (coloringAlgo == "MIS" || graph->GetDomainMap()->lib() == Xpetra::UseTpetra) {
    SubFactoryMonitor sfm(*this, "MIS", currentLevel);
    TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getComm()->getSize() != 1, std::invalid_argument, "MIS on more than 1 MPI rank is not supported");
    DoMISNaive(*graph, myColors, numColors);
  } else {
    SubFactoryMonitor sfm(*this, "GraphColoring", currentLevel);
    TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getComm()->getSize() != 1, std::invalid_argument, "KokkosKernels graph coloring on more than 1 MPI rank is not supported");
    DoGraphColoring(*graph, myColors, numColors);
  }

#ifdef CMS_DUMP
  {
    int rank = graph->GetDomainMap()->getComm()->getRank();

    printf("[%d,%d] num colors %d\n", rank, currentLevel.GetLevelID(), numColors);

    std::ofstream ofs(std::string("m_colors_") + std::to_string(currentLevel.GetLevelID()) + std::string("_") + std::to_string(rank) + std::string(".dat"), std::ofstream::out);
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(ofs));
    *fancy << myColors();
  }
  {
    A->getRowMap()->getComm()->barrier();
  }

#endif

  /* ============================================================= */
  /* Phase 2 : Mark the C-Points                                   */
  /* ============================================================= */
  LO num_c_points = 0, num_d_points = 0, num_f_points = 0;
  if (fc_splitting.is_null()) {
    // We just have a coloring, so we need to generate a splitting
    auto boundaryNodes       = graph->GetBoundaryNodeMap();
    fc_splitting             = LocalOrdinalVectorFactory::Build(A->getRowMap());
    ArrayRCP<LO> myPointType = fc_splitting->getDataNonConst(0);
    for (LO i = 0; i < (LO)myColors.size(); i++) {
      if (boundaryNodes[i]) {
        myPointType[i] = DIRICHLET_PT;
        num_d_points++;
      } else if ((LO)myColors[i] == 1) {
        myPointType[i] = C_PT;
        num_c_points++;
      } else
        myPointType[i] = F_PT;
    }
    num_f_points = (LO)myColors.size() - num_d_points - num_c_points;
  } else {
    // If we read the splitting off disk, we just need to count
    ArrayRCP<LO> myPointType = fc_splitting->getDataNonConst(0);

    for (LO i = 0; i < (LO)myPointType.size(); i++) {
      if (myPointType[i] == DIRICHLET_PT)
        num_d_points++;
      else if (myPointType[i] == C_PT)
        num_c_points++;
    }
    num_f_points = (LO)myPointType.size() - num_d_points - num_c_points;
  }

  /* Output statistics on c/f/d points */
  if (GetVerbLevel() & Statistics1) {
    // NOTE: We batch the communication here
    GO l_counts[] = {(GO)num_c_points, (GO)num_f_points, (GO)num_d_points};
    GO g_counts[3];

    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 3, l_counts, g_counts);
    GetOStream(Statistics1) << "ClassicalMapFactory(" << coloringAlgo << "): C/F/D = " << g_counts[0] << "/" << g_counts[1] << "/" << g_counts[2] << std::endl;
  }

  /* Generate the Coarse map */
  RCP<const Map> coarseMap;
  {
    SubFactoryMonitor sfm(*this, "Coarse Map", currentLevel);
    GenerateCoarseMap(*A->getRowMap(), num_c_points, coarseMap);
  }

  Set(currentLevel, "FC Splitting", fc_splitting);
  Set(currentLevel, "CoarseMap", coarseMap);
}

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GenerateCoarseMap(const Map& fineMap, LO num_c_points, RCP<const Map>& coarseMap) const {
  // FIXME: Assumes scalar PDE
  std::vector<size_t> stridingInfo_(1);
  stridingInfo_[0]   = 1;
  GO domainGIDOffset = 0;

  coarseMap = StridedMapFactory::Build(fineMap.lib(),
                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                       num_c_points,
                                       fineMap.getIndexBase(),
                                       stridingInfo_,
                                       fineMap.getComm(),
                                       domainGIDOffset);
}

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DoGraphColoring(const LWGraph& graph, ArrayRCP<LO>& myColors_out, LO& numColors) const {
  const ParameterList& pL = GetParameterList();
  using graph_t           = typename LWGraph_kokkos::local_graph_type;
  using KernelHandle      = KokkosKernels::Experimental::
      KokkosKernelsHandle<typename graph_t::row_map_type::value_type,
                          typename graph_t::entries_type::value_type,
                          typename graph_t::entries_type::value_type,
                          typename graph_t::device_type::execution_space,
                          typename graph_t::device_type::memory_space,
                          typename graph_t::device_type::memory_space>;
  KernelHandle kh;

  // Leave gc algorithm choice as the default
  kh.create_graph_coloring_handle();

  // Get the distance-1 graph coloring handle
  auto coloringHandle = kh.get_graph_coloring_handle();

  // Set the distance-1 coloring algorithm to use
  if (pL.get<bool>("aggregation: deterministic") == true) {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_SERIAL);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "serial") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_SERIAL);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "vertex based") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_VB);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "vertex based bit array") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_VBBIT);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based bit array" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "vertex based color set") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_VBCS);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based color set" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "vertex based deterministic") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_VBD);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based deterministic" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "vertex based deterministic bit array") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_VBDBIT);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based deterministic bit array" << std::endl;
  } else if (pL.get<std::string>("aggregation: coloring algorithm") == "edge based") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_EB);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: edge based" << std::endl;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unrecognized distance 1 coloring algorithm");
  }

  // Create device views for graph rowptrs/colinds
  size_t numRows = graph.GetNodeNumVertices();
  // auto graphLWK  = dynamic_cast<const LWGraph_kokkos*>(&graph);
  auto graphLW = dynamic_cast<const LWGraph*>(&graph);
  TEUCHOS_TEST_FOR_EXCEPTION(!graphLW, std::invalid_argument, "Graph is not a LWGraph object");
  // Run d1 graph coloring
  // Assume that the graph is symmetric so row map/entries and col map/entries are the same

  // if (graphLWK) {
  //   KokkosGraph::Experimental::graph_color(&kh,
  //                                          numRows,
  //                                          numRows,  // FIXME: This should be the number of columns
  //                                          graphLWK->getRowPtrs(),
  //                                          graphLWK->getEntries(),
  //                                          true);
  // } else
  if (graphLW) {
    auto rowptrs = graphLW->getRowPtrs();
    auto entries = graphLW->getEntries();
    KokkosGraph::Experimental::graph_color(&kh,
                                           numRows,
                                           numRows,  // FIXME: This should be the number of columns
                                           rowptrs,
                                           entries,
                                           true);
  }

  // Extract the colors and store them in the aggregates
  auto myColors_d = coloringHandle->get_vertex_colors();
  numColors       = static_cast<LO>(coloringHandle->get_num_colors());

  // Copy back to host
  auto myColors_h = Kokkos::create_mirror_view(myColors_d);
  myColors_out.resize(myColors_h.size());
  Kokkos::View<LO*, Kokkos::LayoutLeft, Kokkos::HostSpace> myColors_v(&myColors_out[0], myColors_h.size());
  Kokkos::deep_copy(myColors_v, myColors_h);

  // clean up coloring handle
  kh.destroy_graph_coloring_handle();

}  // end DoGraphColoring

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DoMISNaive(const LWGraph& graph, ArrayRCP<LO>& myColors, LO& numColors) const {
  // This is a fall-back routine for when we don't have Kokkos or when it isn't initialized
  // We just do greedy MIS because this is easy to write.

  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  LO MIS        = Teuchos::ScalarTraits<LO>::one();

  // FIXME: Not efficient
  myColors.resize(0);
  myColors.resize(graph.GetNodeNumVertices(), LO_INVALID);
  auto boundaryNodes = graph.GetBoundaryNodeMap();
  LO Nrows           = (LO)graph.GetNodeNumVertices();

  for (LO row = 0; row < Nrows; row++) {
    if (boundaryNodes[row])
      continue;
    auto indices              = graph.getNeighborVertices(row);
    bool has_colored_neighbor = false;
    for (LO j = 0; !has_colored_neighbor && j < (LO)indices.length; j++) {
      // FIXME: This does not handle ghosting correctly
      if (myColors[indices(j)] == MIS)
        has_colored_neighbor = true;
    }
    if (!has_colored_neighbor)
      myColors[row] = MIS;
  }
  numColors = 1;
}

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DoDistributedGraphColoring(RCP<const LWGraph>& graph, ArrayRCP<LO>& myColors_out, LO& numColors) const {
#ifdef HAVE_MUELU_ZOLTAN2
  //  const ParameterList& pL = GetParameterList();
  Teuchos::ParameterList params;
  params.set("color_choice", "FirstFit");
  params.set("color_method", "D1");
  //  params.set("color_choice", colorMethod);
  //  params.set("color_method", colorAlg);
  //  params.set("verbose", verbose);
  //  params.set("serial_threshold",serialThreshold);
  // params.set("recolor_degrees",recolorDegrees);

  // Do the coloring via Zoltan2
  using GraphAdapter = MueLuGraphBaseAdapter<LWGraph>;
  GraphAdapter z_adapter(graph);

  // We need to provide the MPI Comm, or else we wind up using the default (eep!)
  Zoltan2::ColoringProblem<GraphAdapter> problem(&z_adapter, &params, graph->GetDomainMap()->getComm());
  problem.solve();
  Zoltan2::ColoringSolution<GraphAdapter>* soln = problem.getSolution();
  ArrayRCP<int> colors                          = soln->getColorsRCP();
  numColors                                     = (LO)soln->getNumColors();

  // Assign the Array RCP or Copy Out
  // FIXME:  This probably won't work if LO!=int
  if (std::is_same<LO, int>::value)
    myColors_out = colors;
  else {
    myColors_out.resize(colors.size());
    for (LO i = 0; i < (LO)myColors_out.size(); i++)
      myColors_out[i] = (LO)colors[i];
  }

  /*

  printf("CMS: numColors = %d\ncolors = ",numColors);
  for(int i=0;i<colors.size(); i++)
    printf("%d ",colors[i]);
  printf("\n");

  */
#endif  // ifdef HAVE_MUELU_ZOLTAN2
}

}  // namespace MueLu

#endif /* MUELU_CLASSICALMAPFACTORY_DEF_HPP_ */
