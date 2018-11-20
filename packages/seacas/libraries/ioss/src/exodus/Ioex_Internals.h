/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
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
 */

#ifndef IOSS_Ioex_Internals_h
#define IOSS_Ioex_Internals_h

#include "Ioss_ParallelUtils.h" // for ParallelUtils
#include <cstdint>              // for int64_t
#include <cstring>              // for strcpy, strncpy
#include <exodusII.h>           // for MAX_LINE_LENGTH, etc
#include <string>               // for string
#include <vector>               // for vector
namespace Ioss {
  class EdgeBlock;
} // namespace Ioss
namespace Ioss {
  class EdgeSet;
} // namespace Ioss
namespace Ioss {
  class ElementBlock;
} // namespace Ioss
namespace Ioss {
  class ElementSet;
} // namespace Ioss
namespace Ioss {
  class FaceBlock;
} // namespace Ioss
namespace Ioss {
  class FaceSet;
} // namespace Ioss
namespace Ioss {
  class NodeBlock;
} // namespace Ioss
namespace Ioss {
  class NodeSet;
} // namespace Ioss
namespace Ioss {
  class SideBlock;
} // namespace Ioss
namespace Ioss {
  class SideSet;
} // namespace Ioss

using entity_id = int64_t;

namespace Ioss {
} // namespace Ioss
/*!
 * This set of classes provides a thin wrapper around the exodusII
 * internals.  It supplants several of the exodusII API calls in
 * order to avoid ncredef calls which totally rewrite the existing
 * database and can be very expensive.  These routines provide all
 * required variable, dimension, and attribute definitions to the
 * underlying netcdf file with only a single ncredef call.
 *
 * To use the application must create an Internals instance
 * and call the Internals::write_meta_data() function.  This
 * function requires several classes as arguments including:
 * <ul>
 * <li> Mesh -- defines mesh global metadata
 * <li> Block -- defines metadata for each block
 * <li> NodeSet -- defines metadata for each nodeset
 * <li> SideSet -- defines metadata for each sideset
 * <li> CommunicationMetaData -- global metadata relating to
 * parallel info.
 * </ul>
 *
 * Calling Internals::write_meta_data(), replaces the
 * following exodusII and nemesis API calls:
 * <ul>
 * <li> ex_put_init(),
 * <li> ex_put_elem_block(),
 * <li> ex_put_node_set_param(),
 * <li> ex_put_side_set_param(),
 * <li> ne_put_init_info(),
 * <li> ne_put_loadbal_param(),
 * <li> ne_put_cmap_params(),
 * </ul>
 */
namespace Ioex {
  struct NodeBlock
  {
    NodeBlock()                       = default;
    NodeBlock(const NodeBlock &other) = default;
    explicit NodeBlock(const Ioss::NodeBlock &other);

    NodeBlock &operator=(const NodeBlock &other);

    ~NodeBlock() = default;

    bool operator==(const NodeBlock &) const;
    bool operator!=(const NodeBlock &other) const { return !(*this == other); }

    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     localOwnedCount{0};
    int64_t     attributeCount{0};
    int64_t     procOffset{0};

  private:
  };

  struct EdgeBlock
  {
    EdgeBlock() { std::strcpy(elType, ""); }

    EdgeBlock(const EdgeBlock &other)
        : name(other.name), id(other.id), entityCount(other.entityCount),
          nodesPerEntity(other.nodesPerEntity), attributeCount(other.attributeCount),
          procOffset(other.procOffset)
    {
      std::strcpy(elType, other.elType);
    }

    explicit EdgeBlock(const Ioss::EdgeBlock &other);

    EdgeBlock &operator=(const EdgeBlock &other);

    ~EdgeBlock() = default;

    bool operator==(const EdgeBlock & /*other*/) const;
    bool operator!=(const EdgeBlock &other) const { return !(*this == other); }

    char        elType[MAX_STR_LENGTH + 1]{};
    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     nodesPerEntity{0};
    int64_t     attributeCount{0};
    int64_t     procOffset{0};

  private:
  };

  struct FaceBlock
  {
    FaceBlock()
        : name(""), id(0), entityCount(0), nodesPerEntity(0), edgesPerEntity(0), attributeCount(0),
          procOffset(0)
    {
      std::strcpy(elType, "");
    }

    FaceBlock(const FaceBlock &other)
        : name(other.name), id(other.id), entityCount(other.entityCount),
          nodesPerEntity(other.nodesPerEntity), edgesPerEntity(other.edgesPerEntity),
          attributeCount(other.attributeCount), procOffset(other.procOffset)
    {
      std::strcpy(elType, other.elType);
    }

    explicit FaceBlock(const Ioss::FaceBlock &other);

    FaceBlock &operator=(const FaceBlock &other);

    ~FaceBlock() = default;

    bool operator==(const FaceBlock & /*other*/) const;
    bool operator!=(const FaceBlock &other) const { return !(*this == other); }

    char        elType[MAX_STR_LENGTH + 1]{};
    std::string name;
    entity_id   id;
    int64_t     entityCount;
    int64_t     nodesPerEntity;
    int64_t     edgesPerEntity;
    int64_t     attributeCount;
    int64_t     procOffset;

  private:
  };

  struct ElemBlock
  {
    ElemBlock()
        : name(""), id(0), entityCount(0), nodesPerEntity(0), edgesPerEntity(0), facesPerEntity(0),
          attributeCount(0), offset_(-1), procOffset(0)
    {
      std::strcpy(elType, "");
    }

    ElemBlock(const ElemBlock &other)
        : name(other.name), id(other.id), entityCount(other.entityCount),
          nodesPerEntity(other.nodesPerEntity), edgesPerEntity(other.edgesPerEntity),
          facesPerEntity(other.facesPerEntity), attributeCount(other.attributeCount),
          offset_(other.offset_), procOffset(other.procOffset)
    {
      std::strcpy(elType, other.elType);
    }

    explicit ElemBlock(const Ioss::ElementBlock &other);

    ElemBlock &operator=(const ElemBlock &other);

    ~ElemBlock() = default;

    bool operator==(const ElemBlock & /*other*/) const;
    bool operator!=(const ElemBlock &other) const { return !(*this == other); }

    char        elType[MAX_STR_LENGTH + 1]{};
    std::string name;
    entity_id   id;
    int64_t     entityCount;
    int64_t     nodesPerEntity;
    int64_t     edgesPerEntity;
    int64_t     facesPerEntity;
    int64_t     attributeCount;
    int64_t     offset_;
    int64_t     procOffset;
  };

  struct NodeSet
  {
    NodeSet()                     = default;
    NodeSet(const NodeSet &other) = default;
    explicit NodeSet(const Ioss::NodeSet &other);
    bool operator==(const NodeSet & /*other*/) const;
    bool operator!=(const NodeSet &other) const { return !(*this == other); }

    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     localOwnedCount{0};
    int64_t     attributeCount{0};
    int64_t     dfCount{0};
    int64_t     procOffset{0};
  };

  struct EdgeSet
  {
    EdgeSet()                     = default;
    EdgeSet(const EdgeSet &other) = default;
    explicit EdgeSet(const Ioss::EdgeSet &other);
    bool operator==(const EdgeSet & /*other*/) const;
    bool operator!=(const EdgeSet &other) const { return !(*this == other); }

    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     attributeCount{0};
    int64_t     dfCount{0};
    int64_t     procOffset{0};
  };

  struct FaceSet
  {
    FaceSet()                     = default;
    FaceSet(const FaceSet &other) = default;
    explicit FaceSet(const Ioss::FaceSet &other);
    bool operator==(const FaceSet & /*other*/) const;
    bool operator!=(const FaceSet &other) const { return !(*this == other); }

    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     attributeCount{0};
    int64_t     dfCount{0};
    int64_t     procOffset{0};
  };

  struct ElemSet
  {
    ElemSet()                     = default;
    ElemSet(const ElemSet &other) = default;
    explicit ElemSet(const Ioss::ElementSet &other);
    bool operator==(const ElemSet & /*other*/) const;
    bool operator!=(const ElemSet &other) const { return !(*this == other); }

    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     attributeCount{0};
    int64_t     dfCount{0};
    int64_t     procOffset{0};
  };

  struct SideSet
  {
    SideSet() = default;
    explicit SideSet(const Ioss::SideBlock &other);
    explicit SideSet(const Ioss::SideSet &other);
    bool operator==(const SideSet & /*other*/) const;
    bool operator!=(const SideSet &other) const { return !(*this == other); }

    std::string name{};
    entity_id   id{0};
    int64_t     entityCount{0};
    int64_t     dfCount{0};
    int64_t     procOffset{0};
    int64_t     dfProcOffset{0};
  };

  struct CommunicationMap
  {
    CommunicationMap() = default;
    CommunicationMap(entity_id the_id, int64_t count, char the_type)
        : id(the_id), entityCount(count), type(the_type)
    {
    }
    bool      operator==(const CommunicationMap & /*other*/) const;
    bool      operator!=(const CommunicationMap &other) const { return !(*this == other); }
    entity_id id{0};
    int64_t   entityCount{0};
    char      type{'U'}; // 'n' for node, 'e' for element
  };

  struct CommunicationMetaData
  {
    CommunicationMetaData()
        : processorId(0), processorCount(0), globalNodes(0), globalElements(0),
          globalElementBlocks(0), globalNodeSets(0), globalSideSets(0), nodesInternal(0),
          nodesBorder(0), nodesExternal(0), elementsInternal(0), elementsBorder(0),
          outputNemesis(false)
    {
    }

    std::vector<CommunicationMap> nodeMap;
    std::vector<CommunicationMap> elementMap;
    int                           processorId;
    int                           processorCount;
    int64_t                       globalNodes;
    int64_t                       globalElements;
    int64_t                       globalElementBlocks;
    int64_t                       globalNodeSets;
    int64_t                       globalSideSets;
    int64_t                       nodesInternal;
    int64_t                       nodesBorder;
    int64_t                       nodesExternal;
    int64_t                       elementsInternal;
    int64_t                       elementsBorder;
    bool                          outputNemesis;

  private:
    CommunicationMetaData(const CommunicationMetaData &);
  };

  class Redefine
  {
  public:
    explicit Redefine(int exoid);
    Redefine(const Redefine &from) = delete;
    Redefine &operator=(const Redefine &from) = delete;
    ~Redefine();

  private:
    int exodusFilePtr;
  };

  class Mesh
  {
  public:
    Mesh() : title(), dimensionality(0), file_per_processor(true) {}

    Mesh(int dim, char *the_title, bool file_pp) : dimensionality(dim), file_per_processor(file_pp)
    {
      std::strncpy(title, the_title, MAX_LINE_LENGTH + 1);
      title[MAX_LINE_LENGTH] = '\0';
    }

    void populate(Ioss::Region *region);

    char title[MAX_LINE_LENGTH + 1]{};
    int  dimensionality;
    bool file_per_processor;

    std::vector<NodeBlock> nodeblocks;
    std::vector<EdgeBlock> edgeblocks;
    std::vector<FaceBlock> faceblocks;
    std::vector<ElemBlock> elemblocks;
    std::vector<NodeSet>   nodesets;
    std::vector<EdgeSet>   edgesets;
    std::vector<FaceSet>   facesets;
    std::vector<ElemSet>   elemsets;
    std::vector<SideSet>   sidesets;
    CommunicationMetaData  comm;
  };

  class Internals
  {
  public:
    Internals(int exoid, int maximum_name_length, const Ioss::ParallelUtils &util);
    Internals(const Internals &from) = delete;
    Internals &operator=(const Internals &from) = delete;

    int initialize_state_file(Mesh &mesh, const ex_var_params &var_params,
                              const std::string &base_file_name);

    int write_meta_data(Mesh &mesh);

    /*!  A restart file may contain an attribute which contains
     *   information about the processor count and current processor id
     *   * when the file was written.  This code checks whether that
     *   information matches the current processor count and id.  If it
     *   * exists, but doesn't match, a warning message is printed.
     *   Eventually, this will be used to determine whether certain
     *   decomposition-related data in the file is valid or has been
     *   invalidated by a join/re-spread to a different number of
     *   processors.
     */

  private:
    void get_global_counts(Mesh &mesh);

    int put_metadata(const Mesh &mesh, const CommunicationMetaData &comm);
    int put_metadata(const std::vector<NodeBlock> &nodeblocks, bool count_only = false);
    int put_metadata(const std::vector<EdgeBlock> &blocks, bool count_only = false);
    int put_metadata(const std::vector<FaceBlock> &blocks, bool count_only = false);
    int put_metadata(const std::vector<ElemBlock> &blocks, bool count_only = false);

    int put_metadata(const std::vector<NodeSet> &nodesets, bool count_only = false);
    int put_metadata(const std::vector<EdgeSet> &edgesets, bool count_only = false);
    int put_metadata(const std::vector<FaceSet> &facesets, bool count_only = false);
    int put_metadata(const std::vector<ElemSet> &elemsets, bool count_only = false);

    int put_metadata(const std::vector<SideSet> &sidesets, bool count_only = false);

    int put_non_define_data(const CommunicationMetaData &comm);
    int put_non_define_data(const std::vector<NodeBlock> &nodeblocks);
    int put_non_define_data(const std::vector<EdgeBlock> &blocks);
    int put_non_define_data(const std::vector<FaceBlock> &blocks);
    int put_non_define_data(const std::vector<ElemBlock> &blocks);

    int put_non_define_data(const std::vector<NodeSet> &nodesets);
    int put_non_define_data(const std::vector<EdgeSet> &edgesets);
    int put_non_define_data(const std::vector<FaceSet> &facesets);
    int put_non_define_data(const std::vector<ElemSet> &elemsets);

    int put_non_define_data(const std::vector<SideSet> &sidesets);

    int max_name_length() const { return maximumNameLength; }

    int                 exodusFilePtr;
    int                 nodeMapVarID[3];
    int                 elementMapVarID[2];
    int                 commIndexVar{0};
    int                 elemCommIndexVar{0};
    int                 maximumNameLength{32};
    Ioss::ParallelUtils parallelUtil;
  };
} // namespace Ioex
#endif /* IOSS_Ioex_Internals_h */
