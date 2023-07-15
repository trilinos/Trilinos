/*
 * Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "ioexnl_export.h"

#include <Ioss_CoordinateFrame.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Utils.h>

#include <cassert>
#include <exodusII.h>
#include <set>
#include <string>
#include <vector>

#define EXU_USE_HOPSCOTCH
#if defined EXU_USE_HOPSCOTCH
#include <hopscotch_map.h>
#elif defined EXU_USE_ROBIN
#include <robin_map.h>
#endif

// Contains code that is common between the file-per-processor and
// parallel exodus and base exodus classes.

namespace Ioss {
  class GroupingEntity;
  using CoordinateFrameContainer = std::vector<CoordinateFrame>;
} // namespace Ioss

namespace Ioexnl {
  using EntityIdSet = std::set<std::pair<int64_t, int64_t>>;
  using SideSetSet  = std::set<std::string>;
  using SideSetMap  = std::map<std::string, const std::string, std::less<>>;

  using NameTopoKey = std::pair<std::string, const Ioss::ElementTopology *>;
  struct IOEXNL_EXPORT NameTopoKeyCompare
  {
    bool operator()(const NameTopoKey &lhs, const NameTopoKey &rhs) const
    {
      assert(lhs.second != nullptr);
      assert(rhs.second != nullptr);
      return lhs.first < rhs.first ||
             (!(rhs.first < lhs.first) && lhs.second->name() < rhs.second->name());
    }
  };

  struct IOEXNL_EXPORT NameTopoKeyHash
  {
    size_t operator()(const NameTopoKey &name_topo) const
    {
      return std::hash<std::string>{}(name_topo.first) +
             std::hash<size_t>{}((size_t)name_topo.second);
    }
  };

#if defined EXU_USE_HOPSCOTCH
  using TopologyMap = tsl::hopscotch_map<NameTopoKey, int, NameTopoKeyHash>;
#elif defined EXU_USE_ROBIN
  using TopologyMap = tsl::robin_map<NameTopoKey, int, NameTopoKeyHash>;
#else
  // This is the original method that was used in IOSS prior to using hopscotch or robin map.
  using TopologyMap = std::map<NameTopoKey, int, NameTopoKeyCompare>;
#endif

  IOEXNL_EXPORT const char *Version();
  IOEXNL_EXPORT bool        check_processor_info(const std::string &filename, int exodusFilePtr,
                                                 int processor_count, int processor_id);

  IOEXNL_EXPORT Ioss::EntityType map_exodus_type(ex_entity_type type);
  IOEXNL_EXPORT ex_entity_type   map_exodus_type(Ioss::EntityType type);

  IOEXNL_EXPORT void update_last_time_attribute(int exodusFilePtr, double value);
  IOEXNL_EXPORT bool read_last_time_attribute(int exodusFilePtr, double *value);

  IOEXNL_EXPORT bool    type_match(const std::string &type, const char *substring);
  IOEXNL_EXPORT int64_t extract_id(const std::string &name_id);
  IOEXNL_EXPORT bool    set_id(const Ioss::GroupingEntity *entity, Ioexnl::EntityIdSet *idset);
  IOEXNL_EXPORT int64_t get_id(const Ioss::GroupingEntity *entity, Ioexnl::EntityIdSet *idset);
  IOEXNL_EXPORT void    decode_surface_name(Ioexnl::SideSetMap &fs_map, Ioexnl::SideSetSet &fs_set,
                                            const std::string &name);
  IOEXNL_EXPORT void    fix_bad_name(char *name);

  IOEXNL_EXPORT void exodus_error(int exoid, int lineno, const char *function,
                                  const char *filename);
  IOEXNL_EXPORT void exodus_error(int exoid, int lineno, const char *function, const char *filename,
                                  const std::string &extra);

  IOEXNL_EXPORT int add_map_fields(int exoid, Ioss::ElementBlock *block, int64_t my_element_count,
                                   size_t name_length);

  IOEXNL_EXPORT void add_coordinate_frames(int exoid, Ioss::Region *region);
  IOEXNL_EXPORT void write_coordinate_frames(int                                   exoid,
                                             const Ioss::CoordinateFrameContainer &frames);

  IOEXNL_EXPORT bool find_displacement_field(Ioss::NameList             &fields,
                                             const Ioss::GroupingEntity *block, int ndim,
                                             std::string *disp_name);

  IOEXNL_EXPORT std::string get_entity_name(int exoid, ex_entity_type type, int64_t id,
                                            const std::string &basename, int length,
                                            bool &db_has_name);

  IOEXNL_EXPORT void filter_element_list(Ioss::Region *region, Ioss::Int64Vector &elements,
                                         Ioss::Int64Vector &sides, bool remove_omitted_elements);

  IOEXNL_EXPORT bool filter_node_list(Ioss::Int64Vector                &nodes,
                                      const std::vector<unsigned char> &node_connectivity_status);

  template <typename T>
  void filter_node_list(T *data, std::vector<T> &dbvals,
                        const std::vector<int64_t> &active_node_index)
  {
    for (size_t i = 0; i < active_node_index.size(); i++) {
      data[i] = dbvals[active_node_index[i]];
    }
  }

  IOEXNL_EXPORT void filter_element_list(Ioss::Region *region, Ioss::Int64Vector &elements,
                                         Ioss::Int64Vector &sides, bool remove_omitted_elements);

  IOEXNL_EXPORT void separate_surface_element_sides(Ioss::Int64Vector &element,
                                                    Ioss::Int64Vector &sides, Ioss::Region *region,
                                                    Ioexnl::TopologyMap   &topo_map,
                                                    Ioexnl::TopologyMap   &side_map,
                                                    Ioss::SurfaceSplitType split_type,
                                                    const std::string     &surface_name);

  IOEXNL_EXPORT void         write_reduction_attributes(int exoid, const Ioss::GroupingEntity *ge);
  template <typename T> void write_reduction_attributes(int exoid, const std::vector<T *> &entities)
  {
    // For the entity, write all "reduction attributes"
    for (const auto &ge : entities) {
      write_reduction_attributes(exoid, ge);
    }
  }
} // namespace Ioexnl
