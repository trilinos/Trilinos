// Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#define NO_NETCDF_2
#include "CJ_ObjectType.h"
#include <array>
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
    size_t count(ObjectType type) const
    {
      switch (type) {
      case ObjectType::EBLK: return blockCount;
      case ObjectType::NSET: return nodesetCount;
      case ObjectType::SSET: return sidesetCount;
      case ObjectType::NODE: return nodeCount;
      case ObjectType::ELEM: return elementCount;
      case ObjectType::TIME: return timestepCount;
      case ObjectType::DIM: return dimensionality;
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
    size_t entity_count() const { return elementCount; }

    IntVector                truthTable{};
    std::vector<std::string> attributeNames{};
    std::string              name_{};
    ex_entity_id             id{0};
    size_t                   elementCount{0};
    size_t                   nodesPerElement{0};
    size_t                   attributeCount{0};
    size_t                   offset_{0};
    mutable size_t           position_{0};
    std::string              elType{};
  };

  template <typename INT> struct NodeSet
  {
    IntVector    truthTable{};
    ex_entity_id id{0};
    size_t       nodeCount{0};
    size_t       dfCount{0};
    size_t       offset_{0};
    size_t       position_{0};
    std::string  name_{};

    std::vector<INT> nodeSetNodes{};
    std::vector<INT> nodeOrderMap{};
    DistVector       distFactors{};

    size_t entity_count() const { return nodeCount; }

    void dump() const
    {
      fmt::print("NodeSet {}, Name: '{}', {} nodes, {} df,\torder = {}\n", id, name_,
                 fmt::group_digits(nodeCount), fmt::group_digits(dfCount), position_);
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

  template <typename INT> struct SideSet
  {
    IntVector    truthTable{};
    ex_entity_id id{0};
    size_t       sideCount{0};
    size_t       dfCount{0};
    size_t       offset_{0};
    size_t       position_{0};
    std::string  name_{};

    std::vector<INT> elems{};
    std::vector<INT> sides{};

    // For conjoin only. Maps the location (of elems, sides, vars) within this sideset into
    // the location in the corresponding global sideset
    std::vector<INT> elemOrderMap{};
    DistVector       distFactors{};

    size_t entity_count() const { return sideCount; }

    void dump() const
    {
      fmt::print("SideSet {}, Name: '{}', {} sides, {} df\toffset = {}, order = {}\n", id, name_,
                 fmt::group_digits(sideCount), fmt::group_digits(dfCount), offset_, position_);
    }
  };

  struct CommunicationMap
  {
    CommunicationMap(size_t the_id, size_t count, char the_type)
        : id(the_id), entityCount(count), type(the_type)
    {
    }
    ex_entity_id id{0};
    size_t       entityCount{0};
    char         type{'U'}; // 'n' for node, 'e' for element
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
