/*
 * Copyright(C) 2010-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef SEACAS_ExodusEntity_H
#define SEACAS_ExodusEntity_H

#define NO_NETCDF_2
#include "EP_ObjectType.h"
#include <copy_string_cpp.h>
#include <exodusII.h>
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
    Mesh()
        : dimensionality(0), nodeCount(0), elementCount(0), blockCount(0), nodesetCount(0),
          sidesetCount(0), needNodeMap(true), needElementMap(true)
    {
    }

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

    std::string title;
    int         dimensionality;
    int64_t     nodeCount;
    int64_t     elementCount;
    int         blockCount;
    int         nodesetCount;
    int         sidesetCount;
    bool        needNodeMap;
    bool        needElementMap;
  };

  class Block
  {
  public:
    Block()
        : name_(""), id(0), elementCount(0), nodesPerElement(0), attributeCount(0), offset_(0),
          position_(0)
    {
      copy_string(elType, "");
    }

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
    std::string              name_;
    std::vector<std::string> attributeNames;
    int64_t                  id;
    int64_t                  elementCount;
    int                      nodesPerElement;
    int                      attributeCount;
    int64_t                  offset_;
    int                      position_;

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
    NodeSet() : id(0), nodeCount(0), dfCount(0), offset_(0), position_(-1), name_("") {}

    ex_entity_id id;
    int64_t      nodeCount;
    int64_t      dfCount;
    int64_t      offset_;
    int          position_;
    std::string  name_;

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
      for (int64_t i = 0; i < nodeCount; i++) {
        std::cerr << nodeOrderMap[i] << ", ";
      }
      std::cerr << "\n";
    }
  };

  typedef std::pair<int64_t, int64_t> Side;
  template <typename INT> class SideSet
  {
  public:
    SideSet() : id(0), sideCount(0), dfCount(0), offset_(-1), position_(-1), name_("") {}

    ex_entity_id id;
    int64_t      sideCount;
    int64_t      dfCount;
    int64_t      offset_;
    int          position_;
    std::string  name_;

    std::vector<INT> elems;
    std::vector<INT> sides;
    DistVector       distFactors;

    size_t entity_count() const { return sideCount; }

    void dump() const
    {
      std::cerr << "SideSet " << id << ", Name: " << name_ << ", " << sideCount << " sides, "
                << dfCount << " df\toffset = " << offset_ << ", order = " << position_ << "\n";
    }
  };

  class CommunicationMap
  {
  public:
    CommunicationMap() : id(0), entityCount(0), type('U') {}
    CommunicationMap(int the_id, int64_t count, char the_type)
        : id(the_id), entityCount(count), type(the_type)
    {
    }
    int64_t id;
    int64_t entityCount;
    char    type; // 'n' for node, 'e' for element
  };

  class CommunicationMetaData
  {
  public:
    CommunicationMetaData()
        : processorId(0), processorCount(0), globalNodes(0), globalElements(0), nodesInternal(0),
          nodesBorder(0), nodesExternal(0), elementsInternal(0), elementsBorder(0)
    {
    }

    std::vector<CommunicationMap> nodeMap;
    std::vector<CommunicationMap> elementMap;

    int     processorId;
    int     processorCount;
    int64_t globalNodes;
    int64_t globalElements;
    int64_t nodesInternal;
    int64_t nodesBorder;
    int64_t nodesExternal;
    int64_t elementsInternal;
    int64_t elementsBorder;

  private:
    CommunicationMetaData(const CommunicationMetaData &);
  };
} // namespace Excn
#endif /* SEACAS_ExodusEntity_H */
