// Copyright(C) 1999-2010 National Technology & Engineering Solutions
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
//
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
      "entity_count" // Set in GroupingEntity constructor and can be different accross mpi
                     // processes.
  }};

} // namespace Ioad

#endif
