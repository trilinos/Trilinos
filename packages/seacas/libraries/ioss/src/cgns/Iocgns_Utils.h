// Copyright(C) 1999-2017 National Technology & Engineering Solutions
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

#ifndef IOSS_IOCGNS_UTILS_H
#define IOSS_IOCGNS_UTILS_H

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Region.h>
#include <Ioss_SideBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Utils.h>
#include <cgnslib.h>
#include <ostream>
#include <string>

// Used in Iocgns_DatabaseIO.C and Iocgns_ParallelDatabase.C
#define CGCHECK(funcall)                                                                           \
  if ((funcall) != CG_OK) {                                                                        \
    Iocgns::Utils::cgns_error(cgnsFilePtr, __FILE__, __func__, __LINE__, myProcessor);             \
  }

#define CGCHECKNP(funcall)                                                                         \
  if ((funcall) != CG_OK) {                                                                        \
    Iocgns::Utils::cgns_error(cgnsFilePtr, __FILE__, __func__, __LINE__, -1);                      \
  }

// Used in Iocgns_Decomposition.C
#define CGCHECK2(funcall)                                                                          \
  if ((funcall) != CG_OK) {                                                                        \
    Iocgns::Utils::cgns_error(filePtr, __FILE__, __func__, __LINE__, m_decomposition.m_processor); \
  }

namespace Iocgns {
  class Utils
  {
  public:
    Utils()  = default;
    ~Utils() = default;

    static const size_t CG_CELL_CENTER_FIELD_ID = 1ul << 33;
    static const size_t CG_VERTEX_FIELD_ID      = 1ul << 34;

    static size_t index(const Ioss::Field &field);

    static void cgns_error(int cgnsid, const char *file, const char *function, int lineno,
                           int processor);

    static void set_field_index(const Ioss::Field &field, size_t index, CG_GridLocation_t location);
    static bool is_cell_field(const Ioss::Field &field);

    template <typename INT>
    static void map_cgns_face_to_ioss(const Ioss::ElementTopology *parent_topo, size_t num_to_get,
                                      INT *idata)
    {
      // The {topo}_map[] arrays map from CGNS face# to IOSS face#.
      // See http://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
      // NOTE: '0' for first entry is to account for 1-based face numbering.

      switch (parent_topo->shape()) {
      case Ioss::ElementShape::HEX:
        static int hex_map[] = {0, 5, 1, 2, 3, 4, 6};
        for (size_t i = 0; i < num_to_get; i++) {
          idata[2 * i + 1] = hex_map[idata[2 * i + 1]];
        }
        break;

      case Ioss::ElementShape::TET:
        static int tet_map[] = {0, 4, 1, 2, 3};
        for (size_t i = 0; i < num_to_get; i++) {
          idata[2 * i + 1] = tet_map[idata[2 * i + 1]];
        }
        break;

      case Ioss::ElementShape::PYRAMID:
        static int pyr_map[] = {0, 5, 1, 2, 3, 4};
        for (size_t i = 0; i < num_to_get; i++) {
          idata[2 * i + 1] = pyr_map[idata[2 * i + 1]];
        }
        break;

      case Ioss::ElementShape::WEDGE:
#if 0
	  static int wed_map[] = {0, 1, 2, 3, 4, 5}; // Same
	  // Not needed -- maps 1 to 1
	  for (size_t i=0; i < num_to_get; i++) {
	    idata[2*i+1] = wed_map[idata[2*i+1]];
	  }
#endif
        break;
      default:;
      }
    }

    static void map_ioss_face_to_cgns(const Ioss::ElementTopology *parent_topo, size_t num_to_get,
                                      std::vector<cgsize_t> &data)
    {
      // The {topo}_map[] arrays map from CGNS face# to IOSS face#.
      // See http://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
      // NOTE: '0' for first entry is to account for 1-based face numbering.

      switch (parent_topo->shape()) {
      case Ioss::ElementShape::HEX:
        static int hex_map[] = {0, 2, 3, 4, 5, 1, 6};
        for (size_t i = 0; i < num_to_get; i++) {
          data[num_to_get * 2 + i] = hex_map[data[num_to_get * 2 + i]];
        }
        break;

      case Ioss::ElementShape::TET:
        static int tet_map[] = {0, 2, 3, 4, 1};
        for (size_t i = 0; i < num_to_get; i++) {
          data[num_to_get * 2 + i] = tet_map[data[num_to_get * 2 + i]];
        }
        break;

      case Ioss::ElementShape::PYRAMID:
        static int pyr_map[] = {0, 2, 3, 4, 5, 1};
        for (size_t i = 0; i < num_to_get; i++) {
          data[num_to_get * 2 + i] = pyr_map[data[num_to_get * 2 + i]];
        }
        break;

      case Ioss::ElementShape::WEDGE:
#if 0
	  static int wed_map[] = {0, 1, 2, 3, 4, 5}; // Same
	  // Not needed -- maps 1 to 1
	  for (size_t i=0; i < num_to_get; i++) {
	    data[num_to_get * 2 + i] = wed_map[data[num_to_get * 2 + i]];
	  }
#endif
        break;
      default:;
      }
    }

    static void          write_flow_solution_metadata(int file_ptr, Ioss::Region *region, int state,
                                                      int *vertex_solution_index,
                                                      int *cell_center_solution_index);
    static int           find_solution_index(int cgnsFilePtr, int base, int zone, int step,
                                             CG_GridLocation_t location);
    static CG_ZoneType_t check_zone_type(int cgnsFilePtr);
    static void          common_write_meta_data(int file_ptr, const Ioss::Region &region,
                                                std::vector<size_t> &zone_offset);
    static size_t        resolve_nodes(Ioss::Region &region, int my_processor, bool is_parallel);
    static void          resolve_shared_nodes(Ioss::Region &region, int my_processor);

    static CG_ElementType_t map_topology_to_cgns(const std::string &name);
    static std::string      map_cgns_to_topology_type(CG_ElementType_t type);
    static void             add_sidesets(int cgnsFilePtr, Ioss::DatabaseIO *db);
    static void add_structured_boundary_conditions(int cgnsFilePtr, Ioss::StructuredBlock *block);
    static void finalize_database(int cgnsFilePtr, const std::vector<double> &timesteps,
                                  Ioss::Region *region, int myProcessor);
    static int get_step_times(int cgnsFilePtr, std::vector<double> &timesteps, Ioss::Region *region,
                              double timeScaleFactor, int myProcessor);
    static void add_transient_variables(int cgnsFilePtr, const std::vector<double> &timesteps,
                                        Ioss::Region *region, int myProcessor);
  };
} // namespace Iocgns

#endif
