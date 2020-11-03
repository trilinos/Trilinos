// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioad_Constants_h
#define IOSS_Ioad_Constants_h

#include <map>
#include <set>
#include <string>
#include <vector>

namespace Ioad {
  // Constant variables
  const std::string     Schema_version_string = "IOSS_adios_version";
  const std::string     Sideblock_separator   = "::";
  const std::string     Name_separator        = "/";
  const std::string     Role_meta             = "role";
  const std::string     Var_type_meta         = "var_type";
  const std::string     Topology_meta         = "topology";
  const std::string     property_meta         = "property_";
  const std::string     Parent_topology_meta  = "parent_topology";
  const std::string     Time_scale_factor     = "time_scale_factor";
  const std::string     Time_meta             = "time";
  const std::string     Processor_id_meta     = "processor_id";
  const std::string     Processor_number_meta = "processor_number";
  const std::string     globals_entity_type   = "globals";
  const std::string     globals_entity_name   = "";
  const std::string     region_name           = "no_name";
  const std::string     original_name         = "original_name";
  constexpr const char *sideblock_names       = "sideblock_names";

  const std::string                                  coordinate_frame_name = "CoordinateFrame";
  const std::map<std::string, std::set<std::string>> Use_transformed_storage_map = {
      {"ElementBlock", {"connectivity_edge", "connectivity_face"}},
      {"FaceBlock", {"connectivity_edge"}}};
  const std::map<std::string, std::set<std::string>> Ignore_fields = {
      {"NodeBlock",
       {"connectivity", "connectivity_raw", "node_connectivity_status", "implicit_ids",
        "mesh_model_coordinates_x", "mesh_model_coordinates_y", "mesh_model_coordinates_z"}},
      {"ElementBlock", {"implicit_ids"}},
      {"FaceBlock", {"connectivity_raw"}},
      {"EdgeBlock", {"connectivity_raw"}},
      {"CommSet", {"ids"}},
      {"SideSet", {"ids"}},
      {"SideBlock", {"side_ids", "ids", "connectivity", "connectivity_raw"}}};
  const std::vector<std::string> Ignore_properties = {{
      "name", // Name is already known as it is how it is encoded in the output file.
      "_base_stk_part_name", "db_name", // Not necessary
      "streaming_status", "streaming",  // Properties added during processing. Should not be saved.
      "entity_count" // Set in GroupingEntity constructor and can be different across mpi
                     // processes.
  }};

} // namespace Ioad

#endif
