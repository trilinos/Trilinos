/*
 * Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

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
    size_t count(ObjectType type) const
    {
      switch (type) {
      case Excn::ObjectType::EBLK: return blockCount;
      case Excn::ObjectType::NSET: return nodesetCount;
      case Excn::ObjectType::SSET: return sidesetCount;
      case Excn::ObjectType::NODE: return nodeCount;
      case Excn::ObjectType::ELEM: return elementCount;
      case Excn::ObjectType::EDGE: return edgeCount;
      case Excn::ObjectType::FACE: return faceCount;
      case Excn::ObjectType::ASSM: return assemblyCount;
      case Excn::ObjectType::EDBLK: return edgeBlockCount;
      case Excn::ObjectType::FABLK: return faceBlockCount;
      default: return 0;
      }
    }

    IntVector truthTable[5];

    std::string title{};
    int         dimensionality{0};
    int64_t     nodeCount{0};
    int64_t     elementCount{0};
    int64_t     edgeCount{0};
    int64_t     faceCount{0};
    int         blockCount{0};
    int         nodesetCount{0};
    int         sidesetCount{0};
    int         assemblyCount{0};
    int         edgeBlockCount{0};
    int         faceBlockCount{0};
    bool        needNodeMap{true};
    bool        needElementMap{true};
  };

  class Assembly
  {
  public:
    size_t     entity_count() const { return entityCount; }
    ObjectType entity_type() const { return type_; }

    ex_entity_id         id{0};
    std::string          name_;
    ObjectType           type_{Excn::ObjectType::UNSET};
    int                  entityCount{0};
    std::vector<int64_t> entityList;
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

    size_t entity_count() const { return elementCount; }

    char                     elType[MAX_STR_LENGTH + 1]{};
    std::string              name_;
    std::vector<std::string> attributeNames{};
    ex_entity_id             id{0};
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
    ex_entity_id id{0};
    int64_t      nodeCount{0};
    int64_t      dfCount{0};
    int64_t      offset_{0};
    int          position_{-1};
    std::string  name_;

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

  template <typename INT> class SideSet
  {
  public:
    ex_entity_id id{0};
    int64_t      sideCount{0};
    int64_t      dfCount{0};
    int64_t      offset_{-1};
    int          position_{-1};
    std::string  name_;

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

  template <typename INT> class EdgeBlock
  {
  public:
    EdgeBlock() { copy_string(elType, ""); }

    EdgeBlock(const EdgeBlock &other)
        : name_(other.name_), id(other.id), edgeCount(other.edgeCount),
          nodesPerEdge(other.nodesPerEdge), attributeCount(other.attributeCount),
          offset_(other.offset_), position_(other.position_)
    {
      copy_string(elType, other.elType);
    }

    char                     elType[MAX_STR_LENGTH + 1]{};
    std::string              name_;
    std::vector<std::string> attributeNames{};
    ex_entity_id             id{0};
    int64_t                  edgeCount{0};
    int                      nodesPerEdge{0};
    int                      attributeCount{0};
    int64_t                  offset_{0};
    int                      position_{0};

    size_t entity_count() const { return edgeCount; }

    void dump() const
    {
      fmt::print(stderr, "EdgeBlock {}, Name: {}, {} edges\n", id, name_, edgeCount);
    }

    EdgeBlock &operator=(const EdgeBlock &other)
    {
      copy_string(elType, other.elType);
      name_          = other.name_;
      id             = other.id;
      edgeCount      = other.edgeCount;
      nodesPerEdge   = other.nodesPerEdge;
      attributeCount = other.attributeCount;
      attributeNames = other.attributeNames;
      offset_        = other.offset_;
      position_      = other.position_;
      return *this;
    }
  };

  template <typename INT> class FaceBlock
  {
  public:
    FaceBlock() { copy_string(elType, ""); }

    FaceBlock(const FaceBlock &other)
        : name_(other.name_), id(other.id), faceCount(other.faceCount),
          nodesPerFace(other.nodesPerFace), attributeCount(other.attributeCount),
          offset_(other.offset_), position_(other.position_)
    {
      copy_string(elType, other.elType);
    }

    char                     elType[MAX_STR_LENGTH + 1]{};
    std::string              name_;
    std::vector<std::string> attributeNames{};
    ex_entity_id             id{0};
    int64_t                  faceCount{0};
    int                      nodesPerFace{0};
    int                      attributeCount{0};
    int64_t                  offset_{0};
    int                      position_{0};

    size_t entity_count() const { return faceCount; }

    void dump() const
    {
      fmt::print(stderr, "FaceBlock {}, Name: {}, {} faces\n", id, name_, faceCount);
    }

    FaceBlock &operator=(const FaceBlock &other)
    {
      copy_string(elType, other.elType);
      name_          = other.name_;
      id             = other.id;
      faceCount      = other.faceCount;
      nodesPerFace   = other.nodesPerFace;
      attributeCount = other.attributeCount;
      attributeNames = other.attributeNames;
      offset_        = other.offset_;
      position_      = other.position_;
      return *this;
    }
  };

  class CommunicationMap
  {
  public:
    CommunicationMap(int the_id, int64_t count, char the_type)
        : id(the_id), entityCount(count), type(the_type)
    {
    }
    ex_entity_id id{0};
    int64_t      entityCount{0};
    char         type{'U'}; // 'n' for node, 'e' for element
  };

  class CommunicationMetaData
  {
  public:
    CommunicationMetaData()                                              = default;
    CommunicationMetaData(const CommunicationMetaData &)                 = delete;
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
