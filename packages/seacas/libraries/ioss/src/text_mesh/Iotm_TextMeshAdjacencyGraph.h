// Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#pragma once

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <ctype.h>                                   // for toupper
#include <stddef.h>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <map>
#include <set>                                       // for set
#include <sstream>                                   // for operator<<, etc
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector
#include <unordered_map>
#include <sstream>                       // for ostringstream
#include <iostream>
#include <functional>
#include <stdexcept>
#include <numeric>
#include <strings.h>

#include "Iotm_TextMeshFuncs.h"
#include "Iotm_TextMeshDataTypes.h"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace Iotm {
  namespace text_mesh {

    using ErrorHandler = std::function<void(const std::ostringstream &)>;

    template <typename EntityId, typename Topology> struct TextMeshData;

    template <typename EntityId, typename Topology> struct ElementData;

    // Side (or face-face) adjacency is a means to compute neighbors of elements based on
    // node connectivity of the sides bordering neighboring elements. The provided elements
    // are defined in the [ std::vector<ElementData<EntityId, Topology>> elementDataVec ]
    // member of the [ TextMeshData<EntityId, Topology> ] object and since neighbor connectivity
    // based on element id will involve a search on this object, this implementation is based on
    // indexing into this object for performance. Indices are zero-based and sides are one-based
    //
    // The adjacency graph may be created with element filtering based on proc and/or a list of
    // selected element blocks. Any element which is selected by the filter (or all elements)
    // is used as a basis to compute the list of filter nodes which is then used to compute a
    // final list of "local" and "aura" elements which have any connection to these filter nodes.
    // A map is also constructed for the list of these connected elements per filter node
    //
    // For every element in this final list, we can get the list of nodes of per side (based on
    // topology) and using an intersection, find the list of common elements that share all the
    // nodes on that side. These are the neighboring elements and we must determine how to connect
    // the neighboring elements based on solid/shell connectivity rules. The basic rule is that
    // shells act as "blockers" to solid-solid connections
    //
    // Polarity of faces is one determining factor and this is defined based on the relative
    // permutation of the nodes on both faces under consideration
    // +ve/+ve or -ve/-ve => same polarity
    // -ve/+ve or +ve/-ve => opposite polarity
    //
    // {this face / that face}
    // -----------------------
    // Solid/solid connection:
    //    -> always allowed unless there is a shell already connected to either face or both
    //       faces have the same polarity (partially connected solids not allowed)
    // Solid/shell connection:
    //    -> only allowed if both faces are opposite polarity. If the solid face is already
    //    connected
    //       to another solid, break the existing connection to insert the shell connection
    // Shell/solid connection:
    //    -> only allowed if both faces are opposite polarity. If the solid face is already
    //    connected
    //       to another solid, break the existing connection to insert the shell connection
    // Shell/shell connection:
    //    -> only allowed if both faces are opposite polarity AND neither shell face is already
    //       connected to another shell face
    //
    // Stacked shells on a solid element must also obey the rule that all the shells must live on
    // the same processor. This is enforced as it is an error and unsupported in STK
    //
    template <typename EntityId, typename Topology> class SideAdjacencyGraph
    {
    public:
      using IndexType = int64_t;

      static constexpr int       ANY_PROC      = -1;
      static constexpr int       INVALID_SIDE  = -1;
      static constexpr IndexType INVALID_INDEX = -1;

      struct FaceConnection
      {
        FaceConnection()
            : thisSide(INVALID_SIDE), thatElement(INVALID_INDEX), thatSide(INVALID_SIDE){};

        FaceConnection(int thisSide_, IndexType otherElement_, int otherSide_)
            : thisSide(thisSide_), thatElement(otherElement_), thatSide(otherSide_)
        {
        }

        bool operator==(const FaceConnection &rhs) const
        {
          return ((thisSide == rhs.thisSide) && (thatElement == rhs.thatElement) &&
                  (thatSide == rhs.thatSide));
        }

        bool operator!=(const FaceConnection &rhs) const
        {
          return ((thisSide != rhs.thisSide) || (thatElement != rhs.thatElement) ||
                  (thatSide != rhs.thatSide));
        }

        bool operator<(const FaceConnection &rhs) const
        {
          if (thisSide < rhs.thisSide)
            return true;
          else if (thisSide == rhs.thisSide && thatElement < rhs.thatElement)
            return true;
          else if (thisSide == rhs.thisSide && thatElement == rhs.thatElement &&
                   thatSide < rhs.thatSide)
            return true;
          else
            return false;
        }

        int       thisSide;
        IndexType thatElement;
        int       thatSide;
      };

      struct FaceConnections
      {
        FaceConnections() : numSides(0u) {}
        explicit FaceConnections(unsigned numSides_)
            : numSides(numSides_), sideReference(numSides_, 0)
        {
        }

        bool has_any_connection(int thisSide) const
        {
          FaceConnection                              minConnection(thisSide, INVALID_INDEX, 0);
          typename std::set<FaceConnection>::iterator lowerBound =
              std::lower_bound(connections.begin(), connections.end(), minConnection);

          return (lowerBound != connections.end() && lowerBound->thisSide == thisSide);
        }

        bool has_any_connection(int thisSide, IndexType thatElem) const
        {
          FaceConnection                              minConnection(thisSide, thatElem, 0);
          typename std::set<FaceConnection>::iterator lowerBound =
              std::lower_bound(connections.begin(), connections.end(), minConnection);

          return (lowerBound != connections.end() && lowerBound->thisSide == thisSide &&
                  lowerBound->thatElement == thatElem);
        }

        void add(const FaceConnection &connection)
        {
          auto iter = connections.find(connection);
          if (iter == connections.end()) {
            sideReference[connection.thisSide - 1]++;
            connections.insert(connection);
          }
        }

        void remove(const FaceConnection &connection)
        {
          auto iter = connections.find(connection);
          if (iter != connections.end()) {
            sideReference[connection.thisSide - 1]--;
            connections.erase(iter);
          }
        }

        // Number of sides defined by topology for this element entry
        unsigned numSides{0};

        // Number of references to given side based on entries to the connections member: for
        // performance
        std::vector<int> sideReference{};

        // Collection of all face connections for the element. Number might be more than numSides
        // based on solid/shell connections e.e solid connected to multiple stacked shells
        std::set<FaceConnection> connections{};
      };

      virtual size_t get_num_elements() const                                                   = 0;
      virtual int    get_element_proc(const size_t elemIndex) const                             = 0;
      virtual bool   element_has_any_node_on_proc(const size_t elemIndex, int proc) const       = 0;
      virtual const std::string           &get_element_block_name(const size_t elemIndex) const = 0;
      virtual const std::vector<EntityId> &get_element_node_ids(const size_t elemIndex) const   = 0;
      virtual const Topology              &get_element_topology(const size_t elemIndex) const   = 0;
      virtual EntityId                     get_element_id(const size_t elemIndex) const         = 0;

      SideAdjacencyGraph()
      {
        ErrorHandler errorHandler = [](const std::ostringstream &errmsg) {
          default_error_handler(errmsg);
        };
        set_error_handler(errorHandler);
      }

      virtual ~SideAdjacencyGraph() = default;

      void create_graph(int proc = ANY_PROC) { create_graph({}, proc); }

      void create_graph(const std::vector<std::string> &selectedBlocks, int proc = ANY_PROC)
      {
        m_indexGraph.clear();
        std::vector<size_t> localAndAuraElementIndex =
            get_local_and_aura_elements(selectedBlocks, proc);
        std::unordered_map<EntityId, std::set<size_t>> elementsForNode =
            get_elements_for_node_map(localAndAuraElementIndex);
        build_side_connectivity_graph(localAndAuraElementIndex, elementsForNode);
      }

      void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

      void dump(std::ostream &out = std::cout)
      {
        for (const auto &iter : m_indexGraph) {
          if (!iter.second.connections.empty()) {
            out << "Element index: " << iter.first << std::endl;
          }
          for (const auto &entry : iter.second.connections) {
            out << "\tConnected on side: " << entry.thisSide
                << " to element index: " << entry.thatElement << " and side: " << entry.thatSide
                << std::endl;
          }
        }
      }

      void dump(const std::vector<ElementData<EntityId, Topology>> &elemDataVec,
                std::ostream                                       &out = std::cout)
      {
        for (const auto &iter : m_indexGraph) {
          if (!iter.second.connections.empty()) {
            out << "Element: " << elemDataVec[iter.first].identifier << " {"
                << elemDataVec[iter.first].topology << "}" << std::endl;
          }
          for (const auto &entry : iter.second.connections) {
            out << "\tConnected on side: " << entry.thisSide
                << " to element: " << elemDataVec[entry.thatElement].identifier << " {"
                << elemDataVec[entry.thatElement].topology << "}" << " and side: " << entry.thatSide
                << std::endl;
          }
        }
      }

      const FaceConnections &operator[](const IndexType elemIndex) const
      {
        auto it = m_indexGraph.find(elemIndex);

        if (it == m_indexGraph.end()) {
          std::ostringstream errmsg;
          errmsg << "Could not find element index " << elemIndex;
          m_errorHandler(errmsg);
        }

        return it->second;
      }

      size_t size() const { return m_indexGraph.size(); }

      typename std::unordered_map<size_t, FaceConnections>::const_iterator begin() const
      {
        return m_indexGraph.begin();
      }
      typename std::unordered_map<size_t, FaceConnections>::const_iterator end() const
      {
        return m_indexGraph.end();
      }

    protected:
      using OrdinalType     = typename Topology::Ordinal;
      using PermutationType = typename Topology::Permutation;

      struct StringCaseCompLess
      {
        bool operator()(const std::string &lhs, const std::string &rhs)
        {
          return (strcasecmp(lhs.c_str(), rhs.c_str()) < 0);
        };
      };

      // Current element and side under consideration
      struct CurrentAdjacency
      {
        CurrentAdjacency(size_t elementIndex_, int side_) : elementIndex(elementIndex_), side(side_)
        {
        }

        size_t elementIndex;
        int    side;

        // Scratch space that is filled on-the-fly: for performance
        std::vector<int> connectedSides;
      };

      std::set<size_t> get_element_indices_with_common_nodes_on_side(
          const size_t elemIndex, int side,
          const std::unordered_map<EntityId, std::set<size_t>> &elementsForNode)
      {
        std::set<size_t>         neighbors;
        std::vector<OrdinalType> sideNodeIndices =
            get_element_topology(elemIndex).side_topology_node_indices(side);

        const std::vector<EntityId> &elementNodeIds = get_element_node_ids(elemIndex);

        if (!sideNodeIndices.empty()) {
          EntityId firstSideNode = elementNodeIds[sideNodeIndices[0]];
          neighbors              = elementsForNode.at(firstSideNode);
        }

        for (size_t i = 1; i < sideNodeIndices.size(); ++i) {
          OrdinalType             sideNodeIndex          = sideNodeIndices[i];
          EntityId                sideNode               = elementNodeIds[sideNodeIndex];
          const std::set<size_t> &sideNodeElementIndices = elementsForNode.at(sideNode);

          std::set<size_t> intersection;
          std::set_intersection(neighbors.begin(), neighbors.end(), sideNodeElementIndices.begin(),
                                sideNodeElementIndices.end(),
                                std::inserter(intersection, intersection.begin()));

          neighbors = intersection;
        }

        return neighbors;
      }

      std::vector<EntityId> get_side_nodes(const size_t elemIndex, const int side)
      {
        std::vector<OrdinalType> sideNodeIndices =
            get_element_topology(elemIndex).side_topology_node_indices(side);
        std::vector<EntityId>        sideNodes(sideNodeIndices.size());
        const std::vector<EntityId> &elementNodeIds = get_element_node_ids(elemIndex);

        for (size_t i = 0; i < sideNodeIndices.size(); ++i) {
          OrdinalType sideNodeIndex = sideNodeIndices[i];
          sideNodes[i]              = elementNodeIds[sideNodeIndex];
        }

        return sideNodes;
      }

      std::vector<EntityId> get_sorted_side_nodes(const size_t elemIndex, const int side)
      {
        std::vector<EntityId> sideNodes = get_side_nodes(elemIndex, side);
        std::sort(sideNodes.begin(), sideNodes.end());

        return sideNodes;
      }

      // Fill scratch space with side connection info to neighbor element
      void internal_fill_sides_for_connected_element(CurrentAdjacency &adjacency,
                                                     size_t            neighborElementIndex)
      {
        adjacency.connectedSides.clear();

        std::vector<EntityId> sideNodes =
            get_sorted_side_nodes(adjacency.elementIndex, adjacency.side);

        for (int otherSide = 1; otherSide <= get_element_topology(neighborElementIndex).num_sides();
             ++otherSide) {
          std::vector<EntityId> otherSideNodes =
              get_sorted_side_nodes(neighborElementIndex, otherSide);
          if (sideNodes == otherSideNodes) {
            adjacency.connectedSides.push_back(otherSide);
          }
        }
      }

      // Find and validate connected sides to neighbor element
      void fill_sides_for_connected_element(CurrentAdjacency &adjacency,
                                            size_t            neighborElementIndex)
      {
        internal_fill_sides_for_connected_element(adjacency, neighborElementIndex);

        if (adjacency.connectedSides.empty()) {
          std::ostringstream errmsg;
          errmsg << "Neighboring reciprocity check for elements "
                 << get_element_id(adjacency.elementIndex) << " and "
                 << get_element_id(neighborElementIndex) << " failed.";
          m_errorHandler(errmsg);
        }
      }

      // Check to see if two set of nodes are equivalent based on permutation list
      bool equivalent_node_permutation(const std::vector<EntityId>    &controlNodes,
                                       const std::vector<EntityId>    &permutedNodes,
                                       const std::vector<OrdinalType> &permutationOrdinals)
      {
        const size_t numNodes = permutationOrdinals.size();

        if ((numNodes > permutedNodes.size()) || (numNodes > controlNodes.size())) {
          return false;
        }

        bool equivalent = true;
        for (size_t i = 0; equivalent && i < numNodes; ++i) {
          equivalent = controlNodes[permutationOrdinals[i]] == permutedNodes[i];
        }

        return equivalent;
      }

      // Get the permutation that makes one set of nodes equivalent to the other
      std::pair<bool, PermutationType> get_permutation(const Topology              &topology,
                                                       const std::vector<EntityId> &controlNodes,
                                                       const std::vector<EntityId> &permutedNodes,
                                                       PermutationType              numPermutations)
      {
        PermutationType permutation = Topology::InvalidPermutation;
        bool            equivalent  = false;

        if (controlNodes.size() != permutedNodes.size())
          return std::make_pair(equivalent, permutation);

        if (numPermutations > topology.num_permutations()) {
          std::ostringstream errmsg;
          errmsg << "Invalid number of permutations to check: " << numPermutations;
          m_errorHandler(errmsg);
        }

        std::vector<OrdinalType> permutationOrdinals;

        for (PermutationType i = 0; i < numPermutations; ++i) {
          if (topology.fill_permutation_indices(i, permutationOrdinals)) {
            if (equivalent_node_permutation(controlNodes, permutedNodes, permutationOrdinals)) {
              equivalent  = true;
              permutation = i;
              break;
            }
          }
        }

        return std::make_pair(equivalent, permutation);
      }

      std::pair<bool, PermutationType> get_permutation(const Topology              &topology,
                                                       const std::vector<EntityId> &controlNodes,
                                                       const std::vector<EntityId> &permutedNodes)
      {
        return get_permutation(topology, controlNodes, permutedNodes, topology.num_permutations());
      }

      std::pair<bool, PermutationType>
      get_positive_permutation(const Topology &topology, const std::vector<EntityId> &controlNodes,
                               const std::vector<EntityId> &permutedNodes)
      {
        return get_permutation(topology, controlNodes, permutedNodes,
                               topology.num_positive_permutations());
      }

      bool has_same_polarity(const size_t thisElem, const int thisSide, const size_t thatElem,
                             const int thatSide)
      {
        std::vector<EntityId> thisNodes = get_side_nodes(thisElem, thisSide);
        std::vector<EntityId> thatNodes = get_side_nodes(thatElem, thatSide);

        const Topology &sideTopo = get_element_topology(thisElem).side_topology(thisSide);
        std::pair<bool, PermutationType> result =
            get_positive_permutation(sideTopo, thisNodes, thatNodes);

        bool samePolarity = result.first;
        return samePolarity;
      }

      // Check for valid element/side pair
      bool verify_entry(IndexType elemIndex, int side)
      {
        if (INVALID_INDEX == elemIndex)
          return false;

        auto it = m_indexGraph.find(elemIndex);
        if (it == m_indexGraph.end())
          return false;

        int numSides = it->second.numSides;
        if (side < 1 || side > numSides)
          return false;

        return true;
      }

      // Given a current element/side pair, find all other defined reciprocal element connections
      std::vector<FaceConnection> get_reciprocity(CurrentAdjacency &adjacency)
      {
        std::vector<FaceConnection> reciprocity;
        const FaceConnections      &thisEntry = m_indexGraph[adjacency.elementIndex];

        for (const FaceConnection &connectionToThatElement : thisEntry.connections) {
          IndexType thatIndex = connectionToThatElement.thatElement;
          int       thatSide  = connectionToThatElement.thatSide;

          if (!verify_entry(thatIndex, thatSide))
            continue;
          const FaceConnections &thatEntry = m_indexGraph[thatIndex];

          const FaceConnection connectionToThisElement(thatSide, adjacency.elementIndex,
                                                       adjacency.side);
          if (std::binary_search(thatEntry.connections.begin(), thatEntry.connections.end(),
                                 connectionToThisElement)) {
            reciprocity.push_back(connectionToThatElement);
          }
        }

        return reciprocity;
      }

      // Break all current reciprocal connections to this element/side pair
      void break_reciprocal_connections(CurrentAdjacency &adjacency)
      {
        std::vector<FaceConnection> reciprocity = get_reciprocity(adjacency);

        if (!reciprocity.empty()) {
          FaceConnections &thisEntry = m_indexGraph[adjacency.elementIndex];

          for (const FaceConnection &connectionToThatElement : reciprocity) {
            IndexType thatIndex = connectionToThatElement.thatElement;
            int       thatSide  = connectionToThatElement.thatSide;

            FaceConnections &thatEntry = m_indexGraph[thatIndex];
            FaceConnection   connectionToThisElement(thatSide, adjacency.elementIndex,
                                                     adjacency.side);

            thisEntry.remove(connectionToThatElement);
            thatEntry.remove(connectionToThisElement);
          }
        }
      }

      inline bool is_shell_shell_connection(const Topology &thisElemTopology,
                                            const Topology &thatElemTopology)
      {
        return thisElemTopology.is_shell() && thatElemTopology.is_shell();
      }

      inline bool is_shell_shell_connection(const size_t thisElem, const size_t thatElem)
      {
        return is_shell_shell_connection(get_element_topology(thisElem),
                                         get_element_topology(thatElem));
      }

      inline bool is_shell_solid_connection(const Topology &thisElemTopology,
                                            const Topology &thatElemTopology)
      {
        return thisElemTopology.is_shell() && !thatElemTopology.is_shell();
      }

      inline bool is_shell_solid_connection(const size_t thisElem, const size_t thatElem)
      {
        return is_shell_solid_connection(get_element_topology(thisElem),
                                         get_element_topology(thatElem));
      }

      inline bool is_solid_shell_connection(const Topology &thisElemTopology,
                                            const Topology &thatElemTopology)
      {
        return !thisElemTopology.is_shell() && thatElemTopology.is_shell();
      }

      inline bool is_solid_shell_connection(const size_t thisElem, const size_t thatElem)
      {
        return is_solid_shell_connection(get_element_topology(thisElem),
                                         get_element_topology(thatElem));
      }

      inline bool is_solid_solid_connection(const Topology &thisElemTopology,
                                            const Topology &thatElemTopology)
      {
        return !thisElemTopology.is_shell() && !thatElemTopology.is_shell();
      }

      inline bool is_solid_solid_connection(const size_t thisElem, const size_t thatElem)
      {
        return is_solid_solid_connection(get_element_topology(thisElem),
                                         get_element_topology(thatElem));
      }

      using Criterion = std::function<bool(const Topology &topo1, const Topology &topo2)>;

      bool has_connection_type_on_side(size_t thisIndex, int thisSide, Criterion criterion)
      {
        const FaceConnections &thisEntry = m_indexGraph[thisIndex];

        for (const FaceConnection &connection : thisEntry.connections) {
          IndexType thatIndex = connection.thatElement;

          if (connection.thisSide != thisSide || INVALID_INDEX == thatIndex)
            continue;

          if (criterion(get_element_topology(thisIndex), get_element_topology(thatIndex))) {
            return true;
          }
        }

        return false;
      }

      bool has_any_shell_connection_on_side(size_t thisIndex, int thisSide)
      {
        Criterion criterion = [&](const Topology & /* topo1 */, const Topology &topo2) {
          return topo2.is_shell();
        };

        return has_connection_type_on_side(thisIndex, thisSide, criterion);
      }

      bool has_shell_shell_connection_on_side(size_t thisIndex, int thisSide)
      {
        Criterion criterion = [&](const Topology &topo1, const Topology &topo2) {
          return is_shell_shell_connection(topo1, topo2);
        };

        return has_connection_type_on_side(thisIndex, thisSide, criterion);
      }

      bool has_shell_solid_connection_on_side(size_t thisIndex, int thisSide)
      {
        Criterion criterion = [&](const Topology &topo1, const Topology &topo2) {
          return is_shell_solid_connection(topo1, topo2);
        };

        return has_connection_type_on_side(thisIndex, thisSide, criterion);
      }

      bool has_solid_shell_connection_on_side(size_t thisIndex, int thisSide)
      {
        Criterion criterion = [&](const Topology &topo1, const Topology &topo2) {
          return is_solid_shell_connection(topo1, topo2);
        };

        return has_connection_type_on_side(thisIndex, thisSide, criterion);
      }

      bool has_solid_solid_connection_on_side(size_t thisIndex, int thisSide)
      {
        Criterion criterion = [&](const Topology &topo1, const Topology &topo2) {
          return is_solid_solid_connection(topo1, topo2);
        };

        return has_connection_type_on_side(thisIndex, thisSide, criterion);
      }

      void add_connection(CurrentAdjacency &adjacency, size_t connectedElementIndex, int otherSide)
      {
        const size_t thisElem = adjacency.elementIndex;
        const size_t thatElem = connectedElementIndex;

        bool doConnect       = false;
        bool breakConnection = false;

        if (is_shell_solid_connection(thisElem, thatElem)) {
          doConnect = !has_same_polarity(thisElem, adjacency.side, thatElem, otherSide);

          if (has_solid_solid_connection_on_side(connectedElementIndex, otherSide)) {
            breakConnection = doConnect;
          }
        }
        else if (is_solid_shell_connection(thisElem, thatElem)) {
          doConnect = !has_same_polarity(thisElem, adjacency.side, thatElem, otherSide);

          if (has_solid_solid_connection_on_side(adjacency.elementIndex, adjacency.side)) {
            breakConnection = doConnect;
          }
        }
        else if (is_solid_solid_connection(thisElem, thatElem)) {
          doConnect = !has_any_shell_connection_on_side(adjacency.elementIndex, adjacency.side) &&
                      !has_any_shell_connection_on_side(connectedElementIndex, otherSide) &&
                      !has_same_polarity(thisElem, adjacency.side, thatElem, otherSide);
        }

        if (breakConnection) {
          break_reciprocal_connections(adjacency);
        }

        if (doConnect) {
          m_indexGraph[adjacency.elementIndex].add(
              FaceConnection(adjacency.side, connectedElementIndex, otherSide));
          m_indexGraph[connectedElementIndex].add(
              FaceConnection(otherSide, adjacency.elementIndex, adjacency.side));
        }
      }

      void set_side_connectivity(CurrentAdjacency &adjacency, size_t connectedElementIndex)
      {
        fill_sides_for_connected_element(adjacency, connectedElementIndex);

        for (int otherSide : adjacency.connectedSides) {
          add_connection(adjacency, connectedElementIndex, otherSide);
        }
      }

      void enforce_coincident_shell_ownership(IndexType connectedElemIndex1,
                                              IndexType connectedElemIndex2)
      {
        if (get_element_topology(connectedElemIndex1).is_shell() &&
            get_element_topology(connectedElemIndex2).is_shell()) {
          if (get_element_proc(connectedElemIndex1) != get_element_proc(connectedElemIndex2)) {
            std::ostringstream errmsg;
            errmsg << "Invalid proc ownership for co-incident shells "
                   << get_element_id(connectedElemIndex1) << " (proc "
                   << get_element_proc(connectedElemIndex1) << ") and "
                   << get_element_proc(connectedElemIndex2) << " (proc "
                   << get_element_proc(connectedElemIndex2) << ")."
                   << " Co-incident shells must all exist on the same processor";
            m_errorHandler(errmsg);
          }
        }
      }

      void process_side_connectivity(
          CurrentAdjacency                                     &adjacency,
          const std::unordered_map<EntityId, std::set<size_t>> &elementsForNode)
      {
        int side = adjacency.side;

        std::set<size_t> elementIndicesConnectedToSide =
            get_element_indices_with_common_nodes_on_side(adjacency.elementIndex, side,
                                                          elementsForNode);
        for (size_t connectedElementIndex : elementIndicesConnectedToSide) {
          if (connectedElementIndex != adjacency.elementIndex) {
            enforce_coincident_shell_ownership(adjacency.elementIndex, connectedElementIndex);
            set_side_connectivity(adjacency, connectedElementIndex);
          }
        }
      }

      void build_side_connectivity_graph(
          const std::vector<size_t>                            &elementIndices,
          const std::unordered_map<EntityId, std::set<size_t>> &elementsForNode)
      {
        initialize_side_connectivity_graph(elementIndices);

        for (size_t elementIndex : elementIndices) {
          int numSides = get_element_topology(elementIndex).num_sides();
          for (int side = 1; side <= numSides; ++side) {
            if (m_indexGraph[elementIndex].sideReference[side - 1] == 0) {
              CurrentAdjacency adjacency(elementIndex, side);
              process_side_connectivity(adjacency, elementsForNode);
            }
          }
        }
      }

      void initialize_side_connectivity_graph(const std::vector<size_t> &elementIndices)
      {
        for (size_t elementIndex : elementIndices) {
          m_indexGraph[elementIndex] =
              FaceConnections(get_element_topology(elementIndex).num_sides());
        }
      }

      std::unordered_map<EntityId, std::set<size_t>>
      get_elements_for_node_map(const std::vector<size_t> &elementIndices)
      {
        std::unordered_map<EntityId, std::set<size_t>> elementsForNode;
        for (size_t index : elementIndices) {
          for (const EntityId nodeId : get_element_node_ids(index)) {
            elementsForNode[nodeId].insert(index);
          }
        }

        return elementsForNode;
      }

      bool element_is_in_selected_blocks(const size_t                    elemIndex,
                                         const std::vector<std::string> &sortedSelectedBlocks)
      {
        if (sortedSelectedBlocks.empty())
          return true;

        const std::string &partName = get_element_block_name(elemIndex);

        return std::binary_search(sortedSelectedBlocks.begin(), sortedSelectedBlocks.end(),
                                  partName, StringCaseCompLess());
      }

      bool is_selected_element(const size_t                    elemIndex,
                               const std::vector<std::string> &sortedSelectedBlocks, int proc)
      {
        bool isGloballySelected = (ANY_PROC == proc);
        bool isLocallySelected  = (get_element_proc(elemIndex) == proc);
        bool hasLocalNode       = element_has_any_node_on_proc(elemIndex, proc);
        bool isInSelectedBlocks = element_is_in_selected_blocks(elemIndex, sortedSelectedBlocks);

        return isInSelectedBlocks && (isGloballySelected || isLocallySelected || hasLocalNode);
      }

      std::vector<size_t>
      get_local_and_aura_elements(const std::vector<std::string> &selectedBlocks, int proc)
      {
        std::vector<size_t> localAndAuraElementIndex;
        localAndAuraElementIndex.reserve(get_num_elements());

        std::vector<std::string> sortedSelectedBlocks;
        for (const std::string &block : selectedBlocks) {
          sortedSelectedBlocks.push_back(block);
        }
        std::sort(sortedSelectedBlocks.begin(), sortedSelectedBlocks.end(), StringCaseCompLess());

        for (size_t i = 0; i < get_num_elements(); ++i) {
          if (is_selected_element(i, sortedSelectedBlocks, proc)) {
            localAndAuraElementIndex.push_back(i);
          }
        }

        localAndAuraElementIndex.resize(localAndAuraElementIndex.size());
        return localAndAuraElementIndex;
      }

      ErrorHandler                                m_errorHandler;
      std::unordered_map<size_t, FaceConnections> m_indexGraph;
    };

  } // namespace text_mesh
} // namespace Iotm
