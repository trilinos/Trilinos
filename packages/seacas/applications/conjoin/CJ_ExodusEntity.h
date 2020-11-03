// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef SEACAS_ExodusEntity_H
#define SEACAS_ExodusEntity_H

#define NO_NETCDF_2
#include "CJ_ObjectType.h"
#include <copy_string_cpp.h>
#include <cstring>
#include <exodusII.h>
#include <fmt/ostream.h>
#include <iostream>
#include <string>
#include <vector>

namespace Excn {
  using IntVector  = std::vector<int>;
  using DistVector = std::vector<char>;

  template <typename INT> struct Mesh
  {
    Mesh() = default;

    size_t count(ObjectType type) const
    {
      switch (type) {
      case EBLK: return blockCount;
      case NSET: return nodesetCount;
      case SSET: return sidesetCount;
      case NODE: return nodeCount;
      case ELEM: return elementCount;
      case TIME: return timestepCount;
      case DIM: return dimensionality;
      default: return 0;
      }
    }

    std::vector<INT> localNodeToGlobal{};
    std::vector<INT> localElementToGlobal{};

    std::string title{};

    size_t dimensionality{0};
    size_t nodeCount{0};
    size_t elementCount{0};
    size_t blockCount{0};
    size_t nodesetCount{0};
    size_t sidesetCount{0};
    size_t timestepCount{0};

    bool isActive{true};
  };

  struct Block
  {
    Block() { copy_string(elType, ""); }

    Block(const Block &other)
        : truthTable(other.truthTable), attributeNames(other.attributeNames), name_(other.name_),
          id(other.id), elementCount(other.elementCount), nodesPerElement(other.nodesPerElement),
          attributeCount(other.attributeCount), offset_(other.offset_), position_(other.position_)
    {
      copy_string(elType, other.elType);
    }

    ~Block() = default;

    size_t entity_count() const { return elementCount; }

    IntVector                truthTable{};
    std::vector<std::string> attributeNames{};
    std::string              name_{};
    int64_t                  id{0};
    size_t                   elementCount{0};
    size_t                   nodesPerElement{0};
    size_t                   attributeCount{0};
    size_t                   offset_{0};
    size_t                   position_{0};
    char                     elType[MAX_STR_LENGTH + 1]{};

    Block &operator=(const Block &other)
    {
      truthTable      = other.truthTable;
      attributeNames  = other.attributeNames;
      id              = other.id;
      elementCount    = other.elementCount;
      nodesPerElement = other.nodesPerElement;
      attributeCount  = other.attributeCount;
      attributeNames  = other.attributeNames;
      offset_         = other.offset_;
      position_       = other.position_;
      copy_string(elType, other.elType);
      name_ = other.name_;
      return *this;
    }
  };

  template <typename INT> struct NodeSet
  {
    NodeSet() = default;

    IntVector   truthTable{};
    int64_t     id{0};
    size_t      nodeCount{0};
    size_t      dfCount{0};
    size_t      offset_{0};
    size_t      position_{0};
    std::string name_{};

    std::vector<INT> nodeSetNodes{};
    std::vector<INT> nodeOrderMap{};
    DistVector       distFactors{};

    size_t entity_count() const { return nodeCount; }

    void dump() const
    {
      fmt::print("NodeSet {}, Name: '{}', {:n} nodes, {:n} df,\torder = {}\n", id, name_, nodeCount,
                 dfCount, position_);
    }

    void dump_order() const
    {
      dump();
      for (size_t i = 0; i < nodeCount; i++) {
        fmt::print("{}, ", nodeOrderMap[i]);
      }
      fmt::print("\n");
    }
  };

  typedef std::pair<int, int> Side;
  template <typename INT> struct SideSet
  {
    SideSet() = default;

    IntVector   truthTable{};
    int64_t     id{0};
    size_t      sideCount{0};
    size_t      dfCount{0};
    size_t      offset_{0};
    size_t      position_{0};
    std::string name_{};

    std::vector<INT> elems{};
    std::vector<INT> sides{};

    // For conjoin only. Maps the location (of elems, sides, vars) within this sideset into
    // the location in the corresponding global sideset
    std::vector<INT> elemOrderMap{};
    DistVector       distFactors{};

    size_t entity_count() const { return sideCount; }

    void dump() const
    {
      fmt::print("SideSet {}, Name: '{}', {:n} sides, {:n} df\toffset = {}, order = {}\n", id,
                 name_, sideCount, dfCount, offset_, position_);
    }
  };

  struct CommunicationMap
  {
    CommunicationMap() = default;
    CommunicationMap(size_t the_id, size_t count, char the_type)
        : id(the_id), entityCount(count), type(the_type)
    {
    }
    int64_t id{0};
    size_t  entityCount{0};
    char    type{'U'}; // 'n' for node, 'e' for element
  };

  struct CommunicationMetaData
  {
    CommunicationMetaData()                              = default;
    CommunicationMetaData(const CommunicationMetaData &) = delete;

    std::vector<CommunicationMap> nodeMap{};
    std::vector<CommunicationMap> elementMap{};

    size_t processorId{0};
    size_t processorCount{0};
    size_t globalNodes{0};
    size_t globalElements{0};
    size_t nodesInternal{0};
    size_t nodesBorder{0};
    size_t nodesExternal{0};
    size_t elementsInternal{0};
    size_t elementsBorder{0};
  };
} // namespace Excn
#endif /* SEACAS_ExodusEntity_H */
