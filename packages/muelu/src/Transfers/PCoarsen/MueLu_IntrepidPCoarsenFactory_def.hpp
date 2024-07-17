// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_IPCFACTORY_DEF_HPP
#define MUELU_IPCFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_IO.hpp>
#include <sstream>
#include <algorithm>

#include "MueLu_IntrepidPCoarsenFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

#include "Teuchos_ScalarTraits.hpp"

// Intrepid Headers

// Intrepid_HGRAD_HEX_C1_FEM.hpp
// Intrepid_HGRAD_HEX_C2_FEM.hpp
// Intrepid_HGRAD_HEX_Cn_FEM.hpp
// Intrepid_HGRAD_HEX_I2_FEM.hpp
#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
// Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp
// Intrepid_HGRAD_POLY_C1_FEM.hpp
// Intrepid_HGRAD_PYR_C1_FEM.hpp
// Intrepid_HGRAD_PYR_I2_FEM.hpp
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
//#include Intrepid_HGRAD_QUAD_C2_FEM.hpp
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
// Intrepid_HGRAD_TET_C1_FEM.hpp
// Intrepid_HGRAD_TET_C2_FEM.hpp
// Intrepid_HGRAD_TET_Cn_FEM.hpp
// Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp
// Intrepid_HGRAD_TET_COMP12_FEM.hpp
// Intrepid_HGRAD_TRI_C1_FEM.hpp
// Intrepid_HGRAD_TRI_C2_FEM.hpp
// Intrepid_HGRAD_TRI_Cn_FEM.hpp
// Intrepid_HGRAD_TRI_Cn_FEM_ORTH.hpp
// Intrepid_HGRAD_WEDGE_C1_FEM.hpp
// Intrepid_HGRAD_WEDGE_C2_FEM.hpp
// Intrepid_HGRAD_WEDGE_I2_FEM.hpp

// Helper Macro to avoid "unrequested" warnings
#define MUELU_LEVEL_SET_IF_REQUESTED_OR_KEPT(level, ename, entry)                                              \
  {                                                                                                            \
    if (level.IsRequested(ename, this) || level.GetKeepFlag(ename, this) != 0) this->Set(level, ename, entry); \
  }

namespace MueLu {

/*********************************************************************************************************/
namespace MueLuIntrepid {
inline std::string tolower(const std::string &str) {
  std::string data(str);
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  return data;
}

/*********************************************************************************************************/
template <class Basis, class LOFieldContainer, class LocalOrdinal, class GlobalOrdinal, class Node>
void FindGeometricSeedOrdinals(Teuchos::RCP<Basis> basis, const LOFieldContainer &elementToNodeMap,
                               std::vector<std::vector<LocalOrdinal>> &seeds,
                               const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> &rowMap,
                               const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> &columnMap) {
  // For each subcell represented by the elements in elementToNodeMap, we want to identify a globally
  // unique degree of freedom.  Because the other "seed" interfaces in MueLu expect a local ordinal, we
  // store local ordinals in the resulting seeds container.

  // The approach is as follows.  For each element, we iterate through the subcells of the domain topology.
  // We determine which, if any, of these has the lowest global ID owned that is locally owned.  We then insert
  // the local ID corresponding to this in a vector<set<int>> container whose outer index is the spatial dimension
  // of the subcell.  The set lets us conveniently enforce uniqueness of the stored LIDs.

  shards::CellTopology cellTopo = basis->getBaseCellTopology();
  int spaceDim                  = cellTopo.getDimension();
  seeds.clear();
  seeds.resize(spaceDim + 1);
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;

  LocalOrdinal lo_invalid  = Teuchos::OrdinalTraits<LO>::invalid();
  GlobalOrdinal go_invalid = Teuchos::OrdinalTraits<GO>::invalid();

  std::vector<std::set<LocalOrdinal>> seedSets(spaceDim + 1);

  int numCells               = elementToNodeMap.extent(0);
  auto elementToNodeMap_host = Kokkos::create_mirror_view(elementToNodeMap);
  Kokkos::deep_copy(elementToNodeMap_host, elementToNodeMap);
  for (int cellOrdinal = 0; cellOrdinal < numCells; cellOrdinal++) {
    for (int d = 0; d <= spaceDim; d++) {
      int subcellCount = cellTopo.getSubcellCount(d);
      for (int subcord = 0; subcord < subcellCount; subcord++) {
        int dofCount = basis->getDofCount(d, subcord);
        if (dofCount == 0) continue;
        // otherwise, we want to insert the LID corresponding to the least globalID that is locally owned
        GO leastGlobalDofOrdinal     = go_invalid;
        LO LID_leastGlobalDofOrdinal = lo_invalid;
        for (int basisOrdinalOrdinal = 0; basisOrdinalOrdinal < dofCount; basisOrdinalOrdinal++) {
          int basisOrdinal = basis->getDofOrdinal(d, subcord, basisOrdinalOrdinal);
          int colLID       = elementToNodeMap_host(cellOrdinal, basisOrdinal);
          if (colLID != Teuchos::OrdinalTraits<LO>::invalid()) {
            GlobalOrdinal colGID = columnMap.getGlobalElement(colLID);
            LocalOrdinal rowLID  = rowMap.getLocalElement(colGID);
            if (rowLID != lo_invalid) {
              if ((leastGlobalDofOrdinal == go_invalid) || (colGID < leastGlobalDofOrdinal)) {
                // replace with rowLID
                leastGlobalDofOrdinal     = colGID;
                LID_leastGlobalDofOrdinal = rowLID;
              }
            }
          }
        }
        if (leastGlobalDofOrdinal != go_invalid) {
          seedSets[d].insert(LID_leastGlobalDofOrdinal);
        }
      }
    }
  }
  for (int d = 0; d <= spaceDim; d++) {
    seeds[d] = std::vector<LocalOrdinal>(seedSets[d].begin(), seedSets[d].end());
  }
}

/*********************************************************************************************************/
// Syntax [HGRAD|HCURL|HDIV][_| ][HEX|LINE|POLY|PYR|QUAD|TET|TRI|WEDGE][_| ][C|I][1|2|n]
// Inputs:
//  name - name of the intrepid basis to generate
// Outputs:
//  degree - order of resulting discretization
//  return value - Intrepid2 basis correspionding to the name
template <class Scalar, class KokkosExecutionSpace>
Teuchos::RCP<Intrepid2::Basis<KokkosExecutionSpace, Scalar, Scalar>> BasisFactory(const std::string &name, int &degree) {
  using std::string;
  using Teuchos::rcp;
  string myerror("IntrepidBasisFactory: cannot parse string name '" + name + "'");

  // Syntax [HGRAD|HCURL|HDIV][_| ][HEX|LINE|POLY|PYR|QUAD|TET|TRI|WEDGE][_| ][C|I][1|2|n]

  // Get the derivative type
  size_t pos1 = name.find_first_of(" _");
  if (pos1 == 0) throw std::runtime_error(myerror);
  string deriv = tolower(name.substr(0, pos1));
  if (deriv != "hgrad" && deriv != "hcurl" && deriv != "hdiv") throw std::runtime_error(myerror);

  // Get the element type
  pos1++;
  size_t pos2 = name.find_first_of(" _", pos1);
  if (pos2 == 0) throw std::runtime_error(myerror);
  string el = tolower(name.substr(pos1, pos2 - pos1));
  if (el != "hex" && el != "line" && el != "poly" && el != "pyr" && el != "quad" && el != "tet" && el != "tri" && el != "wedge") throw std::runtime_error(myerror);

  // Get the polynomial type
  pos2++;
  string poly = tolower(name.substr(pos2, 1));
  if (poly != "c" && poly != "i") throw std::runtime_error(myerror);

  // Get the degree
  pos2++;
  degree = std::stoi(name.substr(pos2, 1));
  if (degree <= 0) throw std::runtime_error(myerror);

  // FIXME LATER: Allow for alternative point types for Kirby elements
  if (deriv == "hgrad" && el == "quad" && poly == "c") {
    if (degree == 1)
      return rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<KokkosExecutionSpace, Scalar, Scalar>());
    else
      return rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<KokkosExecutionSpace, Scalar, Scalar>(degree, Intrepid2::POINTTYPE_EQUISPACED));
  } else if (deriv == "hgrad" && el == "line" && poly == "c") {
    if (degree == 1)
      return rcp(new Intrepid2::Basis_HGRAD_LINE_C1_FEM<KokkosExecutionSpace, Scalar, Scalar>());
    else
      return rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<KokkosExecutionSpace, Scalar, Scalar>(degree, Intrepid2::POINTTYPE_EQUISPACED));
  }

  // Error out
  throw std::runtime_error(myerror);
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

/*********************************************************************************************************/
// Gets the "lo" nodes nested into a "hi" basis.  Only works on quads and lines for a lo basis of p=1
// Inputs:
//  hi_basis - Higher order Basis
// Outputs:
//  lo_node_in_hi   - std::vector<size_t> of size lo dofs in the reference element, which describes the coindcident hi dots
//  hi_DofCoords    - FC<Scalar> of size (#hi dofs, dim) with the coordinate locations of the hi dofs on the reference element
template <class Scalar, class KokkosDeviceType>
void IntrepidGetP1NodeInHi(const Teuchos::RCP<Intrepid2::Basis<typename KokkosDeviceType::execution_space, Scalar, Scalar>> &hi_basis,
                           std::vector<size_t> &lo_node_in_hi,
                           Kokkos::DynRankView<Scalar, KokkosDeviceType> &hi_DofCoords) {
  typedef typename KokkosDeviceType::execution_space KokkosExecutionSpace;
  // Figure out which unknowns in hi_basis correspond to nodes on lo_basis. This varies by element type.
  size_t degree = hi_basis->getDegree();
  lo_node_in_hi.resize(0);

  if (!rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<KokkosExecutionSpace, Scalar, Scalar>>(hi_basis).is_null()) {
    // HGRAD QUAD Cn: Numbering as per the Kirby convention (straight across, bottom to top)
    lo_node_in_hi.insert(lo_node_in_hi.end(), {0, degree, (degree + 1) * (degree + 1) - 1, degree * (degree + 1)});
  } else if (!rcp_dynamic_cast<Intrepid2::Basis_HGRAD_LINE_Cn_FEM<KokkosExecutionSpace, Scalar, Scalar>>(hi_basis).is_null()) {
    // HGRAD LINE Cn: Numbering as per the Kirby convention (straight across)
    lo_node_in_hi.insert(lo_node_in_hi.end(), {0, degree});
  } else
    throw std::runtime_error("IntrepidPCoarsenFactory: Unknown element type");

  // Get coordinates of the hi_basis dof's
  Kokkos::resize(hi_DofCoords, hi_basis->getCardinality(), hi_basis->getBaseCellTopology().getDimension());
  hi_basis->getDofCoords(hi_DofCoords);
}

/*********************************************************************************************************/
// Given a list of candidates picks a definitive list of "representative" higher order nodes for each lo order node via the "smallest GID" rule
// Input:
//  representative_node_candidates - std::vector<std::vector<size_t> > of lists of "representative candidate" hi dofs for each lo dof
//  hi_elemToNode   - FC<LO> containing the high order element-to-node map
//  hi_columnMap    - Column map of the higher order matrix
// Output:
//  lo_elemToHiRepresentativeNode - FC<LO> of size (# elements, # lo dofs per element) listing the hi unknown chosen as the single representative for each lo unknown for counting purposes
template <class LocalOrdinal, class GlobalOrdinal, class Node, class LOFieldContainer>
void GenerateLoNodeInHiViaGIDs(const std::vector<std::vector<size_t>> &candidates, const LOFieldContainer &hi_elemToNode,
                               RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &hi_columnMap,
                               LOFieldContainer &lo_elemToHiRepresentativeNode) {
  typedef GlobalOrdinal GO;

  // Given: A set of "candidate" hi-DOFs to serve as the "representative" DOF for each lo-DOF on the reference element.
  // Algorithm:  For each element, we choose the lowest GID of the candidates for each DOF to generate the lo_elemToHiRepresentativeNode map

  size_t numElem   = hi_elemToNode.extent(0);
  size_t lo_nperel = candidates.size();
  Kokkos::resize(lo_elemToHiRepresentativeNode, numElem, lo_nperel);

  auto lo_elemToHiRepresentativeNode_host = Kokkos::create_mirror_view(lo_elemToHiRepresentativeNode);
  auto hi_elemToNode_host                 = Kokkos::create_mirror_view(hi_elemToNode);
  Kokkos::deep_copy(hi_elemToNode_host, hi_elemToNode);
  for (size_t i = 0; i < numElem; i++)
    for (size_t j = 0; j < lo_nperel; j++) {
      if (candidates[j].size() == 1)
        lo_elemToHiRepresentativeNode_host(i, j) = hi_elemToNode_host(i, candidates[j][0]);
      else {
        // First we get the GIDs for each candidate
        std::vector<GO> GID(candidates[j].size());
        for (size_t k = 0; k < (size_t)candidates[j].size(); k++)
          GID[k] = hi_columnMap->getGlobalElement(hi_elemToNode_host(i, candidates[j][k]));

        // Find the one with smallest GID
        size_t which = std::distance(GID.begin(), std::min_element(GID.begin(), GID.end()));

        // Record this
        lo_elemToHiRepresentativeNode_host(i, j) = hi_elemToNode_host(i, candidates[j][which]);
      }
    }
  Kokkos::deep_copy(lo_elemToHiRepresentativeNode, lo_elemToHiRepresentativeNode_host);
}

/*********************************************************************************************************/
// Inputs:
//  hi_elemToNode   - FC<LO> containing the high order element-to-node map
//  hi_nodeIsOwned  - std::vector<bool> of size hi's column map, which described hi node ownership
//  lo_elemToHiRepresentativeNode - FC<LO> of size (# elements, # lo dofs per element) listing the hi unknown chosen as the single representative for each lo unknown for counting purposes
// Outputs:
//  lo_elemToNode   - FC<LO> containing the low order element-to-node map.
//  lo_nodeIsOwned  - std::vector<bool> of size lo's (future) column map, which described lo node ownership
//  hi_to_lo_map    - std::vector<LO> of size equal to hi's column map, which contains the lo id each hi idea maps to (or invalid if it doesn't)
//  lo_numOwnedNodes- Number of lo owned nodes
template <class LocalOrdinal, class LOFieldContainer>
void BuildLoElemToNodeViaRepresentatives(const LOFieldContainer &hi_elemToNode,
                                         const std::vector<bool> &hi_nodeIsOwned,
                                         const LOFieldContainer &lo_elemToHiRepresentativeNode,
                                         LOFieldContainer &lo_elemToNode,
                                         std::vector<bool> &lo_nodeIsOwned,
                                         std::vector<LocalOrdinal> &hi_to_lo_map,
                                         int &lo_numOwnedNodes) {
  typedef LocalOrdinal LO;
  using Teuchos::RCP;
  //  printf("CMS:BuildLoElemToNodeViaRepresentatives: hi_elemToNode.rank() = %d hi_elemToNode.size() = %d\n",hi_elemToNode.rank(), hi_elemToNode.size());
  size_t numElem     = hi_elemToNode.extent(0);
  size_t hi_numNodes = hi_nodeIsOwned.size();
  size_t lo_nperel   = lo_elemToHiRepresentativeNode.extent(1);
  Kokkos::resize(lo_elemToNode, numElem, lo_nperel);

  // Start by flagginc the representative nodes
  auto lo_elemToHiRepresentativeNode_host = Kokkos::create_mirror_view(lo_elemToHiRepresentativeNode);
  Kokkos::deep_copy(lo_elemToHiRepresentativeNode_host, lo_elemToHiRepresentativeNode);
  std::vector<bool> is_low_order(hi_numNodes, false);
  for (size_t i = 0; i < numElem; i++)
    for (size_t j = 0; j < lo_nperel; j++) {
      LO id            = lo_elemToHiRepresentativeNode_host(i, j);
      is_low_order[id] = true;  // This can overwrite and that is OK.
    }

  // Count the number of lo owned nodes, generating a local index for lo nodes
  lo_numOwnedNodes   = 0;
  size_t lo_numNodes = 0;
  hi_to_lo_map.resize(hi_numNodes, Teuchos::OrdinalTraits<LO>::invalid());

  for (size_t i = 0; i < hi_numNodes; i++)
    if (is_low_order[i]) {
      hi_to_lo_map[i] = lo_numNodes;
      lo_numNodes++;
      if (hi_nodeIsOwned[i]) lo_numOwnedNodes++;
    }

  // Flag the owned lo nodes
  lo_nodeIsOwned.resize(lo_numNodes, false);
  for (size_t i = 0; i < hi_numNodes; i++) {
    if (is_low_order[i] && hi_nodeIsOwned[i])
      lo_nodeIsOwned[hi_to_lo_map[i]] = true;
  }

  // Translate lo_elemToNode to a lo local index
  auto lo_elemToNode_host = Kokkos::create_mirror_view(lo_elemToNode);
  for (size_t i = 0; i < numElem; i++)
    for (size_t j = 0; j < lo_nperel; j++)
      lo_elemToNode_host(i, j) = hi_to_lo_map[lo_elemToHiRepresentativeNode_host(i, j)];

  // Check for the [E|T]petra column map ordering property, namely LIDs for owned nodes should all appear first.
  // Since we're injecting from the higher-order mesh, it should be true, but we should add an error check & throw in case.
  bool map_ordering_test_passed = true;
  for (size_t i = 0; i < lo_numNodes - 1; i++)
    if (!lo_nodeIsOwned[i] && lo_nodeIsOwned[i + 1])
      map_ordering_test_passed = false;

  if (!map_ordering_test_passed)
    throw std::runtime_error("MueLu::MueLuIntrepid::BuildLoElemToNodeViaRepresentatives failed map ordering test");
  Kokkos::deep_copy(lo_elemToNode, lo_elemToNode_host);
}

/*********************************************************************************************************/
// Inputs:
//  hi_elemToNode   - FC<LO> containing the high order element-to-node map
//  hi_nodeIsOwned  - std::vector<bool> of size hi's column map, which described hi node ownership
//  lo_node_in_hi   - std::vector<size_t> of size lo dofs in the reference element, which describes the coindcident hi dots
//  hi_isDirichlet  - ArrayView<int> of size of hi's column map, which has a 1 if the unknown is Dirichlet and a 0 if it isn't.
// Outputs:
//  lo_elemToNode   - FC<LO> containing the low order element-to-node map.
//  lo_nodeIsOwned  - std::vector<bool> of size lo's (future) column map, which described lo node ownership
//  hi_to_lo_map    - std::vector<LO> of size equal to hi's column map, which contains the lo id each hi idea maps to (or invalid if it doesn't)
//  lo_numOwnedNodes- Number of lo owned nodes
template <class LocalOrdinal, class LOFieldContainer>
void BuildLoElemToNode(const LOFieldContainer &hi_elemToNode,
                       const std::vector<bool> &hi_nodeIsOwned,
                       const std::vector<size_t> &lo_node_in_hi,
                       const Teuchos::ArrayRCP<const int> &hi_isDirichlet,
                       LOFieldContainer &lo_elemToNode,
                       std::vector<bool> &lo_nodeIsOwned,
                       std::vector<LocalOrdinal> &hi_to_lo_map,
                       int &lo_numOwnedNodes) {
  typedef LocalOrdinal LO;
  using Teuchos::RCP;
  LocalOrdinal LOINVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  //  printf("CMS:BuildLoElemToNode: hi_elemToNode.rank() = %d hi_elemToNode.size() = %d\n",hi_elemToNode.rank(), hi_elemToNode.size());

  size_t numElem     = hi_elemToNode.extent(0);
  size_t hi_numNodes = hi_nodeIsOwned.size();

  size_t lo_nperel = lo_node_in_hi.size();
  Kokkos::resize(lo_elemToNode, numElem, lo_nperel);

  // Build lo_elemToNode (in the hi local index ordering) and flag owned ones
  std::vector<bool> is_low_order(hi_numNodes, false);
  auto hi_elemToNode_host = Kokkos::create_mirror_view(hi_elemToNode);
  Kokkos::deep_copy(hi_elemToNode_host, hi_elemToNode);
  auto lo_elemToNode_host = Kokkos::create_mirror_view(lo_elemToNode);
  for (size_t i = 0; i < numElem; i++)
    for (size_t j = 0; j < lo_nperel; j++) {
      LO lid = hi_elemToNode_host(i, lo_node_in_hi[j]);

      // Remove Dirichlet
      if (hi_isDirichlet[lid])
        lo_elemToNode_host(i, j) = LOINVALID;
      else {
        lo_elemToNode_host(i, j)                              = lid;
        is_low_order[hi_elemToNode_host(i, lo_node_in_hi[j])] = true;  // This can overwrite and that is OK.
      }
    }

  // Count the number of lo owned nodes, generating a local index for lo nodes
  lo_numOwnedNodes   = 0;
  size_t lo_numNodes = 0;
  hi_to_lo_map.resize(hi_numNodes, Teuchos::OrdinalTraits<LO>::invalid());

  for (size_t i = 0; i < hi_numNodes; i++)
    if (is_low_order[i]) {
      hi_to_lo_map[i] = lo_numNodes;
      lo_numNodes++;
      if (hi_nodeIsOwned[i]) lo_numOwnedNodes++;
    }

  // Flag the owned lo nodes
  lo_nodeIsOwned.resize(lo_numNodes, false);
  for (size_t i = 0; i < hi_numNodes; i++) {
    if (is_low_order[i] && hi_nodeIsOwned[i])
      lo_nodeIsOwned[hi_to_lo_map[i]] = true;
  }

  // Translate lo_elemToNode to a lo local index
  for (size_t i = 0; i < numElem; i++)
    for (size_t j = 0; j < lo_nperel; j++) {
      if (lo_elemToNode_host(i, j) != LOINVALID)
        lo_elemToNode_host(i, j) = hi_to_lo_map[lo_elemToNode_host(i, j)];
    }
  Kokkos::deep_copy(lo_elemToNode, lo_elemToNode_host);

  // Check for the [E|T]petra column map ordering property, namely LIDs for owned nodes should all appear first.
  // Since we're injecting from the higher-order mesh, it should be true, but we should add an error check & throw in case.
  bool map_ordering_test_passed = true;
  for (size_t i = 0; i < lo_numNodes - 1; i++)
    if (!lo_nodeIsOwned[i] && lo_nodeIsOwned[i + 1])
      map_ordering_test_passed = false;

  if (!map_ordering_test_passed)
    throw std::runtime_error("MueLu::MueLuIntrepid::BuildLoElemToNode failed map ordering test");
}

/*********************************************************************************************************/
// Generates the lo_columnMap
// Input:
//  hi_importer        - Importer from the hi matrix
//  hi_to_lo_map       - std::vector<LO> of size equal to hi's column map, which contains the lo id each hi idea maps to (or invalid if it doesn't)
//  lo_DomainMap       - Domain map for the lo matrix
//  lo_columnMapLength - Number of local columns in the lo column map
// Output:
//  lo_columnMap       - Column map of the lower order matrix
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void GenerateColMapFromImport(const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> &hi_importer, const std::vector<LocalOrdinal> &hi_to_lo_map, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> &lo_domainMap, const size_t &lo_columnMapLength, RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &lo_columnMap) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;
  typedef Xpetra::Map<LO, GO, NO> Map;
  typedef Xpetra::Vector<GO, LO, GO, NO> GOVector;

  GO go_invalid = Teuchos::OrdinalTraits<GO>::invalid();
  LO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();

  RCP<const Map> hi_domainMap = hi_importer.getSourceMap();
  RCP<const Map> hi_columnMap = hi_importer.getTargetMap();
  // Figure out the GIDs of my non-owned P1 nodes
  // HOW: We can build a GOVector(domainMap) and fill the values with either invalid() or the P1 domainMap.GID() for that guy.
  // Then we can use A's importer to get a GOVector(colMap) with that information.

  // NOTE: This assumes rowMap==colMap and [E|T]petra ordering of all the locals first in the colMap
  RCP<GOVector> dvec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(hi_domainMap);
  {
    ArrayRCP<GO> dvec_data = dvec->getDataNonConst(0);
    for (size_t i = 0; i < hi_domainMap->getLocalNumElements(); i++) {
      if (hi_to_lo_map[i] != lo_invalid)
        dvec_data[i] = lo_domainMap.getGlobalElement(hi_to_lo_map[i]);
      else
        dvec_data[i] = go_invalid;
    }
  }

  RCP<GOVector> cvec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(hi_columnMap, true);
  cvec->doImport(*dvec, hi_importer, Xpetra::ADD);

  // Generate the lo_columnMap
  // HOW: We can use the local hi_to_lo_map from the GID's in cvec to generate the non-contiguous colmap ids.
  Array<GO> lo_col_data(lo_columnMapLength);
  {
    ArrayRCP<GO> cvec_data = cvec->getDataNonConst(0);
    for (size_t i = 0, idx = 0; i < hi_columnMap->getLocalNumElements(); i++) {
      if (hi_to_lo_map[i] != lo_invalid) {
        lo_col_data[idx] = cvec_data[i];
        idx++;
      }
    }
  }

  lo_columnMap = Xpetra::MapFactory<LO, GO, NO>::Build(lo_domainMap.lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), lo_col_data(), lo_domainMap.getIndexBase(), lo_domainMap.getComm());
}

/*********************************************************************************************************/
// Generates a list of "representative candidate" hi dofs for each lo dof on the reference element.  This is to be used in global numbering.
// Input:
//  basis                  - The low order basis
//  ReferenceNodeLocations - FC<Scalar> of size (#hidofs, dim) Locations of higher order nodes on the reference element
//  threshold              - tolerance for equivalance testing
// Output:
//  representative_node_candidates - std::vector<std::vector<size_t> > of lists of "representative candidate" hi dofs for each lo dof
template <class Basis, class SCFieldContainer>
void GenerateRepresentativeBasisNodes(const Basis &basis, const SCFieldContainer &ReferenceNodeLocations, const double threshold, std::vector<std::vector<size_t>> &representative_node_candidates) {
  typedef SCFieldContainer FC;
  typedef typename FC::data_type SC;

  // Evaluate the linear basis functions at the Pn nodes
  size_t numFieldsHi = ReferenceNodeLocations.extent(0);
  // size_t dim         = ReferenceNodeLocations.extent(1);
  size_t numFieldsLo = basis.getCardinality();

  FC LoValues("LoValues", numFieldsLo, numFieldsHi);

  basis.getValues(LoValues, ReferenceNodeLocations, Intrepid2::OPERATOR_VALUE);

  Kokkos::fence();  // for kernel in getValues

#if 0
  printf("** LoValues[%d,%d] **\n",(int)numFieldsLo,(int)numFieldsHi);
  for(size_t i=0; i<numFieldsLo; i++) {
    for(size_t j=0; j<numFieldsHi; j++)
      printf("%6.4e ",LoValues(i,j));
    printf("\n");
  }
  printf("**************\n");fflush(stdout);
#endif

  representative_node_candidates.resize(numFieldsLo);
  auto LoValues_host = Kokkos::create_mirror_view(LoValues);
  Kokkos::deep_copy(LoValues_host, LoValues);
  for (size_t i = 0; i < numFieldsLo; i++) {
    // 1st pass: find the max value
    typename Teuchos::ScalarTraits<SC>::magnitudeType vmax = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<SC>::magnitudeType>::zero();
    for (size_t j = 0; j < numFieldsHi; j++)
      vmax = std::max(vmax, Teuchos::ScalarTraits<SC>::magnitude(LoValues_host(i, j)));

    // 2nd pass: Find all values w/i threshhold of target
    for (size_t j = 0; j < numFieldsHi; j++) {
      if (Teuchos::ScalarTraits<SC>::magnitude(vmax - LoValues_host(i, j)) < threshold * vmax)
        representative_node_candidates[i].push_back(j);
    }
  }

  // Sanity check
  for (size_t i = 0; i < numFieldsLo; i++)
    if (!representative_node_candidates[i].size())
      throw std::runtime_error("ERROR: GenerateRepresentativeBasisNodes: No candidates found!");
}

}  // namespace MueLuIntrepid

/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GenerateLinearCoarsening_pn_kirby_to_p1(const LOFieldContainer &hi_elemToNode,
                                                                                                                 const std::vector<bool> &hi_nodeIsOwned,
                                                                                                                 const SCFieldContainer &hi_DofCoords,
                                                                                                                 const std::vector<size_t> &lo_node_in_hi,
                                                                                                                 const Basis &lo_basis,
                                                                                                                 const std::vector<LocalOrdinal> &hi_to_lo_map,
                                                                                                                 const Teuchos::RCP<const Map> &lo_colMap,
                                                                                                                 const Teuchos::RCP<const Map> &lo_domainMap,
                                                                                                                 const Teuchos::RCP<const Map> &hi_map,
                                                                                                                 Teuchos::RCP<Matrix> &P) const {
  typedef SCFieldContainer FC;
  // Evaluate the linear basis functions at the Pn nodes
  size_t numFieldsHi     = hi_elemToNode.extent(1);
  size_t numFieldsLo     = lo_basis.getCardinality();
  LocalOrdinal LOINVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  FC LoValues_at_HiDofs("LoValues_at_HiDofs", numFieldsLo, numFieldsHi);
  lo_basis.getValues(LoValues_at_HiDofs, hi_DofCoords, Intrepid2::OPERATOR_VALUE);
  auto LoValues_at_HiDofs_host = Kokkos::create_mirror_view(LoValues_at_HiDofs);
  Kokkos::deep_copy(LoValues_at_HiDofs_host, LoValues_at_HiDofs);
  Kokkos::fence();  // for kernel in getValues

  typedef typename Teuchos::ScalarTraits<SC>::halfPrecision SClo;
  typedef typename Teuchos::ScalarTraits<SClo>::magnitudeType MT;
  MT effective_zero = Teuchos::ScalarTraits<MT>::eps();

  // Allocate P
  P                   = rcp(new CrsMatrixWrap(hi_map, lo_colMap, numFieldsHi));  // FIXLATER: Need faster fill
  RCP<CrsMatrix> Pcrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Slow-ish fill
  size_t Nelem = hi_elemToNode.extent(0);
  std::vector<bool> touched(hi_map->getLocalNumElements(), false);
  Teuchos::Array<GO> col_gid(1);
  Teuchos::Array<SC> val(1);
  auto hi_elemToNode_host = Kokkos::create_mirror_view(hi_elemToNode);
  Kokkos::deep_copy(hi_elemToNode_host, hi_elemToNode);
  for (size_t i = 0; i < Nelem; i++) {
    for (size_t j = 0; j < numFieldsHi; j++) {
      LO row_lid = hi_elemToNode_host(i, j);
      GO row_gid = hi_map->getGlobalElement(row_lid);
      if (hi_nodeIsOwned[row_lid] && !touched[row_lid]) {
        for (size_t k = 0; k < numFieldsLo; k++) {
          // Get the local id in P1's column map
          LO col_lid = hi_to_lo_map[hi_elemToNode_host(i, lo_node_in_hi[k])];
          if (col_lid == LOINVALID) continue;

          col_gid[0] = {lo_colMap->getGlobalElement(col_lid)};
          val[0]     = LoValues_at_HiDofs_host(k, j);

          // Skip near-zeros
          if (Teuchos::ScalarTraits<SC>::magnitude(val[0]) >= effective_zero)
            P->insertGlobalValues(row_gid, col_gid(), val());
        }
        touched[row_lid] = true;
      }
    }
  }
  P->fillComplete(lo_domainMap, hi_map);
}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GenerateLinearCoarsening_pn_kirby_to_pm(const LOFieldContainer &hi_elemToNode,
                                                                                                                 const std::vector<bool> &hi_nodeIsOwned,
                                                                                                                 const SCFieldContainer &hi_DofCoords,
                                                                                                                 const LOFieldContainer &lo_elemToHiRepresentativeNode,
                                                                                                                 const Basis &lo_basis,
                                                                                                                 const std::vector<LocalOrdinal> &hi_to_lo_map,
                                                                                                                 const Teuchos::RCP<const Map> &lo_colMap,
                                                                                                                 const Teuchos::RCP<const Map> &lo_domainMap,
                                                                                                                 const Teuchos::RCP<const Map> &hi_map,
                                                                                                                 Teuchos::RCP<Matrix> &P) const {
  typedef SCFieldContainer FC;
  // Evaluate the linear basis functions at the Pn nodes
  size_t numFieldsHi = hi_elemToNode.extent(1);
  size_t numFieldsLo = lo_basis.getCardinality();
  FC LoValues_at_HiDofs("LoValues_at_HiDofs", numFieldsLo, numFieldsHi);
  lo_basis.getValues(LoValues_at_HiDofs, hi_DofCoords, Intrepid2::OPERATOR_VALUE);
  auto LoValues_at_HiDofs_host            = Kokkos::create_mirror_view(LoValues_at_HiDofs);
  auto hi_elemToNode_host                 = Kokkos::create_mirror_view(hi_elemToNode);
  auto lo_elemToHiRepresentativeNode_host = Kokkos::create_mirror_view(lo_elemToHiRepresentativeNode);
  Kokkos::deep_copy(LoValues_at_HiDofs_host, LoValues_at_HiDofs);
  Kokkos::deep_copy(hi_elemToNode_host, hi_elemToNode);
  Kokkos::deep_copy(lo_elemToHiRepresentativeNode_host, lo_elemToHiRepresentativeNode);
  Kokkos::fence();  // for kernel in getValues

  typedef typename Teuchos::ScalarTraits<SC>::halfPrecision SClo;
  typedef typename Teuchos::ScalarTraits<SClo>::magnitudeType MT;
  MT effective_zero = Teuchos::ScalarTraits<MT>::eps();

  // Allocate P
  P                   = rcp(new CrsMatrixWrap(hi_map, lo_colMap, numFieldsHi));  // FIXLATER: Need faster fill
  RCP<CrsMatrix> Pcrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Slow-ish fill
  size_t Nelem = hi_elemToNode.extent(0);
  std::vector<bool> touched(hi_map->getLocalNumElements(), false);
  Teuchos::Array<GO> col_gid(1);
  Teuchos::Array<SC> val(1);
  for (size_t i = 0; i < Nelem; i++) {
    for (size_t j = 0; j < numFieldsHi; j++) {
      LO row_lid = hi_elemToNode_host(i, j);
      GO row_gid = hi_map->getGlobalElement(row_lid);
      if (hi_nodeIsOwned[row_lid] && !touched[row_lid]) {
        for (size_t k = 0; k < numFieldsLo; k++) {
          // Get the local id in P1's column map
          LO col_lid = hi_to_lo_map[lo_elemToHiRepresentativeNode_host(i, k)];
          col_gid[0] = {lo_colMap->getGlobalElement(col_lid)};
          val[0]     = LoValues_at_HiDofs_host(k, j);

          // Skip near-zeros
          if (Teuchos::ScalarTraits<SC>::magnitude(val[0]) >= effective_zero)
            P->insertGlobalValues(row_gid, col_gid(), val());
        }
        touched[row_lid] = true;
      }
    }
  }
  P->fillComplete(lo_domainMap, hi_map);
}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("pcoarsen: hi basis");
  SET_VALID_ENTRY("pcoarsen: lo basis");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");

  validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<RCP<const FactoryBase>>("pcoarsen: element to node map", Teuchos::null, "Generating factory of the element to node map");
  return validParamList;
}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level & /* coarseLevel */) const {
  Input(fineLevel, "A");
  Input(fineLevel, "pcoarsen: element to node map");
  Input(fineLevel, "Nullspace");
}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

/*********************************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "P Coarsening", coarseLevel);
  std::string levelIDs     = toString(coarseLevel.GetLevelID());
  const std::string prefix = "MueLu::IntrepidPCoarsenFactory(" + levelIDs + "): ";

  // NOTE: This is hardwired to double on purpose.  See the note below.
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;
  typedef Kokkos::DynRankView<double, typename Node::device_type> FC;

  // Level Get
  RCP<Matrix> A                                                          = Get<RCP<Matrix>>(fineLevel, "A");
  RCP<MultiVector> fineNullspace                                         = Get<RCP<MultiVector>>(fineLevel, "Nullspace");
  Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Acrs = dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> &>(*A);

  if (restrictionMode_) {
    SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
    A = Utilities::Transpose(*A, true);  // build transpose of A explicitly
  }

  // Find the Dirichlet rows in A
  std::vector<LocalOrdinal> A_dirichletRows;
  Utilities::FindDirichletRows(A, A_dirichletRows);

  // Build final prolongator
  RCP<Matrix> finalP;

  // Reuse pattern if available
  RCP<ParameterList> APparams = rcp(new ParameterList);
  if (coarseLevel.IsAvailable("AP reuse data", this)) {
    GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous AP data" << std::endl;

    APparams = coarseLevel.Get<RCP<ParameterList>>("AP reuse data", this);

    if (APparams->isParameter("graph"))
      finalP = APparams->get<RCP<Matrix>>("graph");
  }
  const ParameterList &pL = GetParameterList();

  /*******************/
  // FIXME LATER: Allow these to be manually specified instead of Intrepid
  // Get the Intrepid bases
  // NOTE: To make sure Stokhos works we only instantiate these guys with double.  There's a lot
  // of stuff in the guts of Intrepid2 that doesn't play well with Stokhos as of yet.
  int lo_degree, hi_degree;
  RCP<Basis> hi_basis = MueLuIntrepid::BasisFactory<double, typename Node::device_type::execution_space>(pL.get<std::string>("pcoarsen: hi basis"), hi_degree);
  RCP<Basis> lo_basis = MueLuIntrepid::BasisFactory<double, typename Node::device_type::execution_space>(pL.get<std::string>("pcoarsen: lo basis"), lo_degree);

  // Useful Output
  GetOStream(Statistics1) << "P-Coarsening from basis " << pL.get<std::string>("pcoarsen: hi basis") << " to " << pL.get<std::string>("pcoarsen: lo basis") << std::endl;

  /*******************/
  // Get the higher-order element-to-node map
  const Teuchos::RCP<FCi> Pn_elemToNode = Get<Teuchos::RCP<FCi>>(fineLevel, "pcoarsen: element to node map");

  // Calculate node ownership (the quick and dirty way)
  // NOTE: This exploits two things:
  //  1) domainMap == rowMap
  //  2) Standard [e|t]petra ordering (namely the local unknowns are always numbered first).
  // This routine does not work in general.
  RCP<const Map> rowMap    = A->getRowMap();
  RCP<const Map> colMap    = Acrs.getColMap();
  RCP<const Map> domainMap = A->getDomainMap();
  int NumProc              = rowMap->getComm()->getSize();
  assert(rowMap->isSameAs(*domainMap));
  std::vector<bool> Pn_nodeIsOwned(colMap->getLocalNumElements(), false);
  LO num_owned_rows = 0;
  for (size_t i = 0; i < rowMap->getLocalNumElements(); i++) {
    if (rowMap->getGlobalElement(i) == colMap->getGlobalElement(i)) {
      Pn_nodeIsOwned[i] = true;
      num_owned_rows++;
    }
  }

  // Used in all cases
  FC hi_DofCoords;
  Teuchos::RCP<FCi> P1_elemToNode = rcp(new FCi());

  std::vector<bool> P1_nodeIsOwned;
  int P1_numOwnedNodes;
  std::vector<LO> hi_to_lo_map;

  // Degree-1 variables
  std::vector<size_t> lo_node_in_hi;

  // Degree-n variables
  FCi lo_elemToHiRepresentativeNode;

  // Get Dirichlet unknown information
  RCP<Xpetra::Vector<int, LocalOrdinal, GlobalOrdinal, Node>> hi_isDirichletRow, hi_isDirichletCol;
  Utilities::FindDirichletRowsAndPropagateToCols(A, hi_isDirichletRow, hi_isDirichletCol);

#if 0
    printf("[%d] isDirichletRow = ",A->getRowMap()->getComm()->getRank());
    for(size_t i=0;i<hi_isDirichletRow->getMap()->getLocalNumElements(); i++)
      printf("%d ",hi_isDirichletRow->getData(0)[i]);
    printf("\n");
    printf("[%d] isDirichletCol = ",A->getRowMap()->getComm()->getRank());
    for(size_t i=0;i<hi_isDirichletCol->getMap()->getLocalNumElements(); i++)
      printf("%d ",hi_isDirichletCol->getData(0)[i]);
    printf("\n");
    fflush(stdout);
#endif

  /*******************/
  if (lo_degree == 1) {
    // Get reference coordinates and the lo-to-hi injection list for the reference element
    MueLuIntrepid::IntrepidGetP1NodeInHi(hi_basis, lo_node_in_hi, hi_DofCoords);

    // Generate lower-order element-to-node map
    MueLuIntrepid::BuildLoElemToNode(*Pn_elemToNode, Pn_nodeIsOwned, lo_node_in_hi, hi_isDirichletCol->getData(0), *P1_elemToNode, P1_nodeIsOwned, hi_to_lo_map, P1_numOwnedNodes);
    assert(hi_to_lo_map.size() == colMap->getLocalNumElements());
  } else {
    // Get lo-order candidates
    double threshold = 1e-10;
    std::vector<std::vector<size_t>> candidates;
    Kokkos::resize(hi_DofCoords, hi_basis->getCardinality(), hi_basis->getBaseCellTopology().getDimension());
    hi_basis->getDofCoords(hi_DofCoords);

    MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<Basis, FC>(*lo_basis, hi_DofCoords, threshold, candidates);

    // Generate the representative nodes
    MueLu::MueLuIntrepid::GenerateLoNodeInHiViaGIDs(candidates, *Pn_elemToNode, colMap, lo_elemToHiRepresentativeNode);
    MueLu::MueLuIntrepid::BuildLoElemToNodeViaRepresentatives(*Pn_elemToNode, Pn_nodeIsOwned, lo_elemToHiRepresentativeNode, *P1_elemToNode, P1_nodeIsOwned, hi_to_lo_map, P1_numOwnedNodes);
  }
  MUELU_LEVEL_SET_IF_REQUESTED_OR_KEPT(coarseLevel, "pcoarsen: element to node map", P1_elemToNode);

  /*******************/
  // Generate the P1_domainMap
  // HOW: Since we know how many each proc has, we can use the non-uniform contiguous map constructor to do the work for us
  RCP<const Map> P1_domainMap = MapFactory::Build(rowMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), P1_numOwnedNodes, rowMap->getIndexBase(), rowMap->getComm());
  MUELU_LEVEL_SET_IF_REQUESTED_OR_KEPT(coarseLevel, "CoarseMap", P1_domainMap);

  // Generate the P1_columnMap
  RCP<const Map> P1_colMap;
  if (NumProc == 1)
    P1_colMap = P1_domainMap;
  else
    MueLuIntrepid::GenerateColMapFromImport<LO, GO, NO>(*Acrs.getCrsGraph()->getImporter(), hi_to_lo_map, *P1_domainMap, P1_nodeIsOwned.size(), P1_colMap);

  /*******************/
  // Generate the coarsening
  if (lo_degree == 1)
    GenerateLinearCoarsening_pn_kirby_to_p1(*Pn_elemToNode, Pn_nodeIsOwned, hi_DofCoords, lo_node_in_hi, *lo_basis, hi_to_lo_map, P1_colMap, P1_domainMap, A->getRowMap(), finalP);
  else
    GenerateLinearCoarsening_pn_kirby_to_pm(*Pn_elemToNode, Pn_nodeIsOwned, hi_DofCoords, lo_elemToHiRepresentativeNode, *lo_basis, hi_to_lo_map, P1_colMap, P1_domainMap, A->getRowMap(), finalP);

  /*******************/
  // Zero out the Dirichlet rows in P
  Utilities::ZeroDirichletRows(finalP, A_dirichletRows);

  /*******************/
  // Build the nullspace
  RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P1_domainMap, fineNullspace->getNumVectors());
  finalP->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS);
  Set(coarseLevel, "Nullspace", coarseNullspace);

  // Level Set
  if (!restrictionMode_) {
    // The factory is in prolongation mode
    Set(coarseLevel, "P", finalP);

    APparams->set("graph", finalP);
    MUELU_LEVEL_SET_IF_REQUESTED_OR_KEPT(coarseLevel, "AP reuse data", APparams);

    if (IsPrint(Statistics1)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*finalP, "P", params);
    }
  } else {
    // The factory is in restriction mode
    RCP<Matrix> R;
    {
      SubFactoryMonitor m2(*this, "Transpose P", coarseLevel);
      R = Utilities::Transpose(*finalP, true);
    }

    Set(coarseLevel, "R", R);

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*R, "R", params);
    }
  }

}  // Build()

}  // namespace MueLu

#endif  // MUELU_IPCFACTORY_DEF_HPP
