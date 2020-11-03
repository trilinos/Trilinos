/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SEACAS_ExodusEntity_H
#define SEACAS_ExodusEntity_H

#define NO_NETCDF_2
#include "EP_ObjectType.h"
#include <copy_string_cpp.h>
#include <exodusII.h>
#include <fmt/ostream.h>
#include <iostream>
#include <string>
#include <vector>

namespace Excn {
  using IntVector   = std::vector<int>;
  using Int64Vector = std::vector<int64_t>;
  using DistVector  = std::vector<char>;

  class Mesh
  {
  public:
    Mesh() = default;

    size_t count(ObjectType type) const
    {
      switch (type) {
      case EBLK: return blockCount;
      case NSET: return nodesetCount;
      case SSET: return sidesetCount;
      case NODE: return nodeCount;
      case ELEM: return elementCount;
      default: return 0;
      }
    }

    IntVector truthTable[3];

    std::string title{};
    int         dimensionality{0};
    int64_t     nodeCount{0};
    int64_t     elementCount{0};
    int         blockCount{0};
    int         nodesetCount{0};
    int         sidesetCount{0};
    bool        needNodeMap{true};
    bool        needElementMap{true};
  };

  class Block
  {
  public:
    Block() { copy_string(elType, ""); }

    Block(const Block &other)
        : name_(other.name_), id(other.id), elementCount(other.elementCount),
          nodesPerElement(other.nodesPerElement), attributeCount(other.attributeCount),
          offset_(other.offset_), position_(other.position_)
    {
      copy_string(elType, other.elType);
    }

    ~Block() = default;

    size_t entity_count() const { return elementCount; }

    char                     elType[MAX_STR_LENGTH + 1]{};
    std::string              name_{""};
    std::vector<std::string> attributeNames{};
    int64_t                  id{0};
    int64_t                  elementCount{0};
    int                      nodesPerElement{0};
    int                      attributeCount{0};
    int64_t                  offset_{0};
    int                      position_{0};

    Block &operator=(const Block &other)
    {
      copy_string(elType, other.elType);
      name_           = other.name_;
      id              = other.id;
      elementCount    = other.elementCount;
      nodesPerElement = other.nodesPerElement;
      attributeCount  = other.attributeCount;
      attributeNames  = other.attributeNames;
      offset_         = other.offset_;
      position_       = other.position_;
      return *this;
    }
  };

  template <typename INT> class NodeSet
  {
  public:
    NodeSet() = default;

    ex_entity_id id{0};
    int64_t      nodeCount{0};
    int64_t      dfCount{0};
    int64_t      offset_{0};
    int          position_{-1};
    std::string  name_{""};

    std::vector<INT> nodeSetNodes{};
    std::vector<INT> nodeOrderMap{};
    DistVector       distFactors{};

    size_t entity_count() const { return nodeCount; }

    void dump() const
    {
      fmt::print(stderr, "NodeSet {}, Name: {}, {} nodes, {} df,\torder = {}\n", id, name_,
                 nodeCount, dfCount, position_);
    }

    void dump_order() const
    {
      dump();
      for (int64_t i = 0; i < nodeCount; i++) {
        fmt::print(stderr, "{}, ", nodeOrderMap[i]);
      }
      fmt::print("\n");
    }
  };

  typedef std::pair<int64_t, int64_t> Side;
  template <typename INT> class SideSet
  {
  public:
    SideSet() = default;

    ex_entity_id id{0};
    int64_t      sideCount{0};
    int64_t      dfCount{0};
    int64_t      offset_{-1};
    int          position_{-1};
    std::string  name_{""};

    std::vector<INT> elems{};
    std::vector<INT> sides{};
    DistVector       distFactors{};

    size_t entity_count() const { return sideCount; }

    void dump() const
    {
      fmt::print(stderr, "SideSet {}, Name: {}, {} sides, {} df\toffset = {}, order = {}\n", id,
                 name_, sideCount, dfCount, offset_, position_);
    }
  };

  class CommunicationMap
  {
  public:
    CommunicationMap() = default;
    CommunicationMap(int the_id, int64_t count, char the_type)
        : id(the_id), entityCount(count), type(the_type)
    {
    }
    int64_t id{0};
    int64_t entityCount{0};
    char    type{'U'}; // 'n' for node, 'e' for element
  };

  class CommunicationMetaData
  {
  public:
    CommunicationMetaData()                              = default;
    CommunicationMetaData(const CommunicationMetaData &) = delete;
    CommunicationMetaData &operator=(const CommunicationMetaData &other) = delete;

    std::vector<CommunicationMap> nodeMap{};
    std::vector<CommunicationMap> elementMap{};

    int     processorId{0};
    int     processorCount{0};
    int64_t globalNodes{0};
    int64_t globalElements{0};
    int64_t nodesInternal{0};
    int64_t nodesBorder{0};
    int64_t nodesExternal{0};
    int64_t elementsInternal{0};
    int64_t elementsBorder{0};
  };
} // namespace Excn
#endif /* SEACAS_ExodusEntity_H */
