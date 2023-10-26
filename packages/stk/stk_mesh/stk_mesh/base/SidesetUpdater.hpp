// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_SIDESET_UPDATER_HPP
#define STK_SIDESET_UPDATER_HPP

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequire
#include <stk_util/environment/Env.hpp>
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector
#include "stk_mesh/base/SideSetEntry.hpp"
#include "stk_mesh/base/SideSetHelper.hpp"
#include "stk_mesh/base/SideSetUtil.hpp"
#include <stk_mesh/baseImpl/elementGraph/GraphTypes.hpp>

namespace stk { namespace mesh {

struct part_compare_by_ordinal {
  bool operator() (const Part * const i, const Part * const j) const {
      if(i == j) return false;
      if(nullptr == i) return true;
      if(nullptr == j) return false;
      return (i->mesh_meta_data_ordinal() < j->mesh_meta_data_ordinal());
  }
};


class SidesetUpdater : public ModificationObserver
{
public:

  SidesetUpdater(BulkData &bulk, Selector& active_selector)
    : ModificationObserver(ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY)
    , m_bulkData(bulk)
    , m_metaData(bulk.mesh_meta_data())
    , m_activeSelector(active_selector)
    , m_isActive(true)
    , m_outputStream(nullptr)
    , m_helper(bulk, active_selector)
    {
      set_warn_about_internal_sideset(false);
    }

    void set_active(bool active) { m_isActive = active; }
    bool get_active_flag() const { return m_isActive; }

    void set_output_stream(std::ostream& ostrm) {
      if(m_outputStream != nullptr) return;

      m_outputStream = &ostrm;
      m_helper.set_output_stream(ostrm);
    }

    void clear_output_stream() {
      m_outputStream = nullptr;
    }

    std::ostream* get_output_stream() { return m_outputStream; }

    void set_debug_file(const std::string& fileNamePrefix)
    {
      if(m_debugOutputFileStream.is_open()) {
        m_debugOutputFileStream.close();
      }

      std::ostringstream oss;
      oss << fileNamePrefix << "." << m_bulkData.parallel_rank();
      m_debugOutputFileStream = std::ofstream(oss.str(), std::ofstream::out);
      STK_ThrowRequireMsg(m_debugOutputFileStream.fail() == false, "Failed to open debug file: " << oss.str());
      set_output_stream(m_debugOutputFileStream);
    }

    void set_warn_about_internal_sideset(bool flag)
    {
      m_warnAboutInternalSideset = flag;
      m_helper.set_warn_about_internal_sideset(flag);
    }

    const impl::ParallelPartInfo& get_parallel_part_info() const { return m_helper.get_parallel_part_info(); }

protected:
    BulkData &m_bulkData;
    MetaData& m_metaData;
    Selector m_activeSelector;
    bool m_isActive = true;
    bool m_warnAboutInternalSideset = false;

    std::ostream* m_outputStream = nullptr;

    SideSetHelper m_helper;

    std::ofstream m_debugOutputFileStream;
};

class ReconstructionSidesetUpdater : public SidesetUpdater
{
public:
  ReconstructionSidesetUpdater(BulkData &bulk, Selector activeSelector)
    : SidesetUpdater(bulk, activeSelector)
    , m_internalSidesetWarningHasBeenIssued(false)
    {
    }

    void entity_added(Entity entity) override;

    void entity_deleted(Entity entity) override;

    void entity_parts_added(Entity entity, const OrdinalVector& parts) override;

    void entity_parts_removed(Entity entity, const OrdinalVector& parts) override;

    void modification_begin_notification() override;

    void finished_modification_end_notification() override;

    void started_modification_end_notification() override;

    void fill_values_to_reduce(std::vector<size_t> &valuesToReduce) override;

    void set_reduced_values(const std::vector<size_t> &reducedValues) override;

    void elements_about_to_move_procs_notification(const EntityProcVec &elemProcPairsToMove) override;

    void elements_moved_procs_notification(const EntityProcVec &elemProcPairsToMove) override;

private:
    std::set<const Part*, part_compare_by_ordinal> m_stkSideSets;
    bool m_internalSidesetWarningHasBeenIssued;
    std::set<const Part*> m_sidesetPartsWithDeletedEntries;

    void tag_sideset(Entity entity, const Part& part);
    void insert_parts(Entity entity, const ConstPartVector& parts);
    void insert_parts(Entity entity, const OrdinalVector& parts);
    void fill_sidesets_element_belongs_to(Entity elem);
    void reconstruct_noninternal_sidesets(const std::vector<size_t> &reducedValues);
    void reconstruct_sidesets();
    void update_sidesets_without_surface_block_mapping();
};

struct RelationUpdate
{
  RelationUpdate(const Entity face_, const Entity element_, const ConnectivityOrdinal ordinal_)
  : face(face_), element(element_), ordinal(ordinal_), removed(false) { }
  RelationUpdate() : face(), element(), ordinal(INVALID_CONNECTIVITY_ORDINAL), removed(false) {}

  Entity face;
  Entity element;
  ConnectivityOrdinal ordinal;
  bool removed;

  inline bool operator==( const RelationUpdate & rhs ) const
  {
    return face == rhs.face && element == rhs.element && ordinal == rhs.ordinal;
  }

  inline bool operator!=( const RelationUpdate & rhs ) const
  {
    return face != rhs.face || element != rhs.element || ordinal != rhs.ordinal;
  }

  inline bool operator<( const RelationUpdate & rhs ) const
  {
    if(face < rhs.face) {
      return true;
    }
    else if(face == rhs.face) {
      if(element < rhs.element) {
        return true;
      }
      else if (element == rhs.element && ordinal < rhs.ordinal) {
        return true;
      }
      else {
        return false;
      }
    }

    return false;
  }
};

struct BlockToSurfaceMapping
{
  BlockToSurfaceMapping()
    : block(nullptr) { }

  BlockToSurfaceMapping(const Part* block_)
    : block(block_) { }

  BlockToSurfaceMapping(const Part* block_, const ConstPartVector& surfaces_)
    : block(block_), surfaces(surfaces_) { }

  const Part* block;
  ConstPartVector surfaces;

  inline bool operator==( const BlockToSurfaceMapping & rhs ) const
  {
    return block->mesh_meta_data_ordinal() == rhs.block->mesh_meta_data_ordinal();
  }

  inline bool operator!=( const BlockToSurfaceMapping & rhs ) const
  {
    return block->mesh_meta_data_ordinal() != rhs.block->mesh_meta_data_ordinal();
  }

  inline bool operator<( const BlockToSurfaceMapping & rhs ) const
  {
    if(block->mesh_meta_data_ordinal() < rhs.block->mesh_meta_data_ordinal()) {
      return true;
    }

    return false;
  }
};

class IncrementalSidesetUpdater : public SidesetUpdater
{
public:
    IncrementalSidesetUpdater(BulkData &bulk, Selector& activeSelector)
    : SidesetUpdater(bulk, activeSelector)
    , m_relationUpdatesAreSorted(true)
    , m_relationUpdatesRemoved(false)
    , m_elemOrSideChangedRankedParts(false)
    {
    }

    void modification_begin_notification() override;

    void started_modification_end_notification() override;

    void finished_modification_end_notification() override;

    void entity_deleted(Entity entity) override;

    void entity_parts_added(Entity entity, const OrdinalVector& parts) override;

    void entity_parts_removed(Entity entity, const OrdinalVector& parts) override;

    void relation_destroyed(Entity from, Entity to, ConnectivityOrdinal ordinal) override;

    void relation_declared(Entity from, Entity to, ConnectivityOrdinal ordinal) override;

private:
    std::vector<SideSet*> m_sidesets;
    std::map<const Part*, Selector> m_cachedSidesetSelectors;
    std::vector<RelationUpdate> m_relationUpdates;
    bool m_relationUpdatesAreSorted;
    bool m_relationUpdatesRemoved;
    ConstPartVector m_surfacesTouchingBlocks;
    std::vector<BlockToSurfaceMapping> m_cachedBlockToSurfaceMapping;
    bool m_elemOrSideChangedRankedParts;

    PartChangeAccumulatorVector m_accumulatedElementPartChanges;

    const Selector& get_cached_sideset_selector(const Part* part);
    const ConstPartVector& get_cached_surfaces_touching_block(const Part& part);

    void update_sideset_vector();
    void update_surfaces_touching_blocks_vector();

    void add_sidesets_from_part(const Part& part, std::vector<SideSet*>& sidesets);
    void add_sidesets_from_parts(const PartVector& parts, std::vector<SideSet*>& sidesets);
    void add_sidesets_from_parts(const ConstPartVector& parts, std::vector<SideSet*>& sidesets);

    void fill_sidesets_from_part_ordinals(const OrdinalVector& parts, std::vector<SideSet*>& sidesets);
    void fill_sidesets_from_parts(const PartVector& parts, std::vector<SideSet*>& sidesets);

    void add_sidesets_and_selectors_from_part(const Part& part, bool checkSubsets, SideSetSelectorVector& sidesets);
    void fill_sidesets_and_selectors_from_part_ordinals(const OrdinalVector& parts, SideSetSelectorVector& sidesets);
    void fill_sidesets_and_selectors_from_parts(const PartVector& parts, SideSetSelectorVector& sidesets);
    void fill_sidesets_and_selectors(SideSetSelectorVector& sidesets);
    void fill_sidesets_and_selectors_from_blocks_touching_surfaces(const OrdinalVector& parts,
                                                                   SideSetSelectorVector& sidesetsAndSelectors);

    void fill_surfaces_touching_block(const Part* part, ConstPartVector& surfacesTouchingBlock);

    void remove_relation(Entity from, Entity to, ConnectivityOrdinal ordinal);
    void add_relation(Entity from, Entity to, ConnectivityOrdinal ordinal);

    SideSet* get_or_create_connectivity_based_sideset(const Part& part, bool checkSubsets);

    void resolve_relation_updates();

    void resolve_faces_affected_by_connected_element_part_change();

    void accumulate_element_block_part_changes(Entity entity, const OrdinalVector& parts);
    void fill_info_for_element_affected_by_block_part_change(Entity entity,
                                                             SideSetSelectorVector& sidesetsAndSelectors,
                                                             std::vector<std::pair<Entity, unsigned>>& sideAndPartOrdinalVec);
};

}} // end stk mesh namespaces

#endif
