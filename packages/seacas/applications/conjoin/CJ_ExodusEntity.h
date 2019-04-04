// Copyright(C) 2009-2010-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#ifndef SEACAS_ExodusEntity_H
#define SEACAS_ExodusEntity_H

#define NO_NETCDF_2
#include "CJ_ObjectType.h"
#include <cstring>
#include <exodusII.h>
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

    std::vector<INT> localNodeToGlobal;
    std::vector<INT> localElementToGlobal;

    std::string title;

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
    Block() { std::strcpy(elType, ""); }

    Block(const Block &other)
        : name_(other.name_), id(other.id), elementCount(other.elementCount),
          nodesPerElement(other.nodesPerElement), attributeCount(other.attributeCount),
          offset_(other.offset_), position_(other.position_)
    {
      std::strcpy(elType, other.elType);
    }

    ~Block() = default;

    size_t entity_count() const { return elementCount; }

    IntVector                truthTable;
    std::vector<std::string> attributeNames;
    std::string              name_;
    int64_t                  id{0};
    size_t                   elementCount{0};
    size_t                   nodesPerElement{0};
    size_t                   attributeCount{0};
    size_t                   offset_{0};
    size_t                   position_{0};
    char                     elType[MAX_STR_LENGTH + 1]{};

    Block &operator=(const Block &other)
    {
      id              = other.id;
      elementCount    = other.elementCount;
      nodesPerElement = other.nodesPerElement;
      attributeCount  = other.attributeCount;
      attributeNames  = other.attributeNames;
      offset_         = other.offset_;
      position_       = other.position_;
      std::strcpy(elType, other.elType);
      name_ = other.name_;
      return *this;
    }
  };

  template <typename INT> struct NodeSet
  {
    NodeSet() = default;

    IntVector   truthTable;
    int64_t     id{0};
    size_t      nodeCount{0};
    size_t      dfCount{0};
    size_t      offset_{0};
    size_t      position_{0};
    std::string name_;

    std::vector<INT> nodeSetNodes;
    std::vector<INT> nodeOrderMap;
    DistVector       distFactors;

    size_t entity_count() const { return nodeCount; }

    void dump() const
    {
      std::cerr << "NodeSet " << id << ", Name: " << name_ << ", " << nodeCount << " nodes, "
                << dfCount << " df,\torder = " << position_ << "\n";
    }

    void dump_order() const
    {
      dump();
      for (size_t i = 0; i < nodeCount; i++) {
        std::cerr << nodeOrderMap[i] << ", ";
      }
      std::cerr << "\n";
    }
  };

  typedef std::pair<int, int> Side;
  template <typename INT> struct SideSet
  {
    SideSet() = default;

    IntVector   truthTable;
    int64_t     id{0};
    size_t      sideCount{0};
    size_t      dfCount{0};
    size_t      offset_{0};
    size_t      position_{0};
    std::string name_;

    std::vector<INT> elems;
    std::vector<INT> sides;

    // For conjoin only. Maps the location (of elems, sides, vars) within this sideset into
    // the location in the corresponding global sideset
    std::vector<INT> elemOrderMap;
    DistVector       distFactors;

    size_t entity_count() const { return sideCount; }

    void dump() const
    {
      std::cerr << "SideSet " << id << ", Name: " << name_ << ", " << sideCount << " sides, "
                << dfCount << " df\toffset = " << offset_ << ", order = " << position_ << "\n";
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

    std::vector<CommunicationMap> nodeMap;
    std::vector<CommunicationMap> elementMap;

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
