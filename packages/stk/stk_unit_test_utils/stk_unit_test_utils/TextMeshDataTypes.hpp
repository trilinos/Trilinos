// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#ifndef TextMeshDataTypes_hpp
#define TextMeshDataTypes_hpp

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

#include "TextMeshFuncs.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace text_mesh {

using ErrorHandler = std::function<void(const std::ostringstream &)>;

template <typename T>
class TopologyMapping
{
 public:
  using Topology = T;

  virtual ~TopologyMapping() {}

  Topology topology(const std::string &textMeshName) const
  {
    auto it = m_nameToTopology.find(textMeshName);
    return (it != m_nameToTopology.end() ? it->second : invalid_topology());
  }

  virtual Topology invalid_topology() const = 0;
  virtual void initialize_topology_map() = 0;

 protected:
  std::unordered_map<std::string, Topology> m_nameToTopology;
};

class PartIdMapping
{
 public:
  PartIdMapping() : m_idsAssigned(false)
  {
    set_error_handler([](const std::ostringstream &errmsg) { default_error_handler(errmsg); });
  }

  void register_part_name(const std::string &name)
  {
    m_partNames.push_back(name);
    handle_block_part(name);
  }

  void register_part_name_with_id(const std::string &name, unsigned id)
  {
    register_part_name(name);
    assign(name, id);
  }

  unsigned get(const std::string &name) const
  {
    if (!m_idsAssigned) assign_ids();
    return get_part_id(name);
  }

  std::string get(unsigned id) const
  {
    if (!m_idsAssigned) assign_ids();
    return get_part_name(id);
  }

  unsigned size() const
  {
    if (!m_idsAssigned) assign_ids();
    return m_ids.size();
  }

  std::vector<std::string> get_part_names_sorted_by_id() const
  {
    if (!m_idsAssigned) assign_ids();

    std::vector<std::string> names;
    names.reserve(m_parts.size());

    for (auto iter : m_parts) {
      names.push_back(iter.second);
    }

    return names;
  }

  bool is_registered(const std::string &name) const { return m_ids.count(name) > 0; }

  const std::vector<std::string> &get_part_names() const { return m_partNames; }

  void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

  const std::string get_group_type() const { return "element block"; }

  void finalize_parse() {
    if (!m_idsAssigned) assign_ids();
  }

 private:
  void handle_block_part(const std::string &name)
  {
    auto result = get_id_from_part_name(name, "BLOCK_");

    if (!result.second) return;

    assign(name, result.first);
  }

  void assign_ids() const
  {
    unsigned nextPartId = 1;
    for (const std::string &name : m_partNames) {
      if (m_ids.find(name) == m_ids.end()) {
        while (is_assigned(nextPartId)) nextPartId++;

        assign(name, nextPartId);
      }
    }

    m_idsAssigned = true;
  }

  void assign(const std::string &name, unsigned id) const
  {
    validate_name_and_id(name, id);
    m_ids[name] = id;
    m_parts[id] = name;
  }

  void validate_name_and_id(const std::string &name, unsigned id) const
  {
    if (is_registered(name)) {
      if (m_ids[name] != id) {
        std::ostringstream errmsg;
        errmsg << "Cannot assign part '" << name << "' two different ids: " << m_ids[name] << " and " << id;
        m_errorHandler(errmsg);
      }
    } else {
      if (is_assigned(id)) {
        std::ostringstream errmsg;
        errmsg << "Part id " << id << " has already been assigned, cannot assign it to part '" << name << "'";
        m_errorHandler(errmsg);
      }
    }
  }

  bool is_assigned(unsigned id) const { return m_parts.count(id) > 0; }

  unsigned get_part_id(const std::string &name) const
  {
    auto it = m_ids.find(name);
    if (it == m_ids.end()) {
      std::ostringstream errmsg;
      errmsg << "PartIdMapping has no ID for invalid part name " << name;
      m_errorHandler(errmsg);
    }
    return it->second;
  }

  std::string get_part_name(unsigned id) const
  {
    auto it = m_parts.find(id);
    if (it == m_parts.end()) {
      std::ostringstream errmsg;
      errmsg << "PartIdMapping has no part name for invalid id " << id;
      m_errorHandler(errmsg);
    }
    return it->second;
  }

  std::vector<std::string> m_partNames{};
  mutable std::unordered_map<std::string, unsigned> m_ids;
  mutable std::map<unsigned, std::string> m_parts;
  mutable bool m_idsAssigned{false};

  ErrorHandler m_errorHandler;
};

template <typename EntityId>
class Coordinates
{
 public:
  Coordinates()
  {
    set_error_handler([](const std::ostringstream &errmsg) { default_error_handler(errmsg); });
  }

  const std::vector<double> &operator[](const EntityId nodeId) const
  {
    auto it = m_nodalCoords.find(nodeId);

    if (it == m_nodalCoords.end()) {
      std::ostringstream errmsg;
      errmsg << "Could not find node id " << nodeId;
      m_errorHandler(errmsg);
    }

    return it->second;
  }

  void set_coordinate_data(
      const unsigned spatialDim, const std::set<EntityId> &nodeIds, const std::vector<double> &coordinates)
  {
    if (!coordinates.empty()) {
      validate_num_coordinates(spatialDim, nodeIds, coordinates);
      fill_coordinate_map(spatialDim, nodeIds, coordinates);
      m_hasCoordinateData = true;
    }
  }

  void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

  bool has_coordinate_data() const { return m_hasCoordinateData; }

 private:
  void validate_num_coordinates(
      const unsigned spatialDim, const std::set<EntityId> &nodeIds, const std::vector<double> &coordinates)
  {
    if (coordinates.size() != nodeIds.size() * spatialDim) {
      std::ostringstream errmsg;
      errmsg << "Number of coordinates: " << coordinates.size() << ", Number of nodes: " << nodeIds.size()
             << ", Spatial dimension: " << spatialDim;
      m_errorHandler(errmsg);
    }
  }

  void fill_coordinate_map(
      const unsigned spatialDim, const std::set<EntityId> &nodeIds, const std::vector<double> &coordinates)
  {
    std::vector<double>::const_iterator coordIter = coordinates.begin();
    for (const EntityId &nodeId : nodeIds) {
      m_nodalCoords[nodeId] = std::vector<double>(coordIter, coordIter + spatialDim);
      coordIter += spatialDim;
    }
  }

  bool m_hasCoordinateData{false};
  std::unordered_map<EntityId, std::vector<double>> m_nodalCoords;
  ErrorHandler m_errorHandler;
};

template <typename EntityId, typename Topology>
struct ElementData {
  int proc;
  EntityId identifier;
  Topology topology;
  std::vector<EntityId> nodeIds{};
  std::string partName = "";

  operator EntityId() const { return identifier; }
};

template <typename EntityId, typename Topology>
struct ElementDataLess {
  bool operator()(const ElementData<EntityId, Topology> &lhs, const ElementData<EntityId, Topology> &rhs)
  {
    return lhs.identifier < rhs.identifier;
  };

  bool operator()(const ElementData<EntityId, Topology> &lhs, const EntityId rhs) { return lhs.identifier < rhs; };

  bool operator()(const EntityId lhs, const ElementData<EntityId, Topology> &rhs) { return lhs < rhs.identifier; };

  bool operator()(const EntityId lhs, const EntityId rhs) { return lhs < rhs; };
};

template <typename EntityId, typename Topology>
class Sidesets;

template <typename EntityId>
class Nodesets;

template <typename EntityId>
class Assemblies;

template <typename EntityId, typename Topology>
struct TextMeshData {
  unsigned spatialDim{0};
  std::vector<ElementData<EntityId, Topology>> elementDataVec{};
  PartIdMapping partIds;
  std::set<EntityId> nodeIds{};
  Coordinates<EntityId> coords;
  Sidesets<EntityId, Topology> sidesets;
  Nodesets<EntityId> nodesets;
  Assemblies<EntityId> assemblies;

  TextMeshData() : spatialDim(0) {}

  void add_element(const ElementData<EntityId, Topology> &elem)
  {
    elementDataVec.push_back(elem);
    for (const EntityId &nodeId : elem.nodeIds) {
      nodeIds.insert(nodeId);
      associate_node_with_proc(nodeId, elem.proc);
    }
  }

  const std::set<EntityId> &nodes_on_proc(int proc) const
  {
    auto it = m_nodesOnProc.find(proc);
    return it != m_nodesOnProc.end() ? it->second : m_emptyNodes;
  }

  unsigned num_nodes_on_proc(int proc) const
  {
    auto it = m_nodesOnProc.find(proc);
    return it != m_nodesOnProc.end() ? it->second.size() : 0;
  }

  const std::set<int> &procs_for_node(const EntityId nodeId) const
  {
    auto it = m_procsForNode.find(nodeId);
    return it != m_procsForNode.end() ? it->second : m_emptyProcs;
  }

 private:
  void associate_node_with_proc(const EntityId nodeId, int proc)
  {
    m_procsForNode[nodeId].insert(proc);
    m_nodesOnProc[proc].insert(nodeId);
  }

  std::unordered_map<EntityId, std::set<int>> m_procsForNode;
  std::unordered_map<int, std::set<EntityId>> m_nodesOnProc;

  std::set<int> m_emptyProcs{};
  std::set<EntityId> m_emptyNodes{};
};


}  // namespace text_mesh

#endif
