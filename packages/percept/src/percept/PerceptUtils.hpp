// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_PerceptUtils_hpp
#define percept_PerceptUtils_hpp

#include <string>
#include <vector>
#include <mpi.h>

namespace stk {
struct topology;
namespace mesh {
  struct Entity;
  class FieldBase;
}
namespace diag {
  class Timer;
}
}

struct CellTopologyData;

namespace percept {
  
  std::vector<std::string> entity_rank_names_and_family_tree();

  void computeCentroid(stk::mesh::Entity entity, double centroid[3], const stk::mesh::FieldBase & coord_field);

  double volume(stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data);

  stk::diag::Timer& rootTimerStructured();

  void printTimersTableStructured();

  double MegaByte(const size_t x);

  void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg, size_t& hwm_sum);
  void get_memory_now_across_processors(            MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg, size_t& hwm_sum);
  std::string print_memory_high_water_mark(MPI_Comm comm);
  std::string print_memory_now(MPI_Comm comm);
  std::string print_memory_both(MPI_Comm comm);

  enum ExecutionStage {INITIALIZE=0, READ_MESH, REFINE_MESH, WRITE_MESH, FINALIZE, INVALID_STAGE};

  class MemoryInfo
  {
    public:
      void fill_memory_info_across_processors(MPI_Comm comm);

      void write_memory_log_header(const std::string& logfile_name);
      void write_memory_log_data  (const std::string& logfile_name, ExecutionStage stage);

    private:
      enum MemorySlot {SUM=0, MIN, MAX, AVG, INVALID_SLOT};

      size_t now[4];
      //size_t hwm[4]; // currently not output

      std::string header[INVALID_SLOT] = {"SUM", "MIN", "MAX", "AVG"};
      std::string stages[INVALID_STAGE] = {"INITIALIZE", "READ_MESH", "REFINE_MESH", "WRITE_MESH", "FINALIZE"};

      const int col_width = 14;
  };

  void write_memory_logfile(MPI_Comm comm, ExecutionStage stage, const std::string& memory_logfile_name);

  void convert_stk_topology_to_ioss_name(
          const stk::topology stk_topo,
          std::string& ioss_topo);
}

#endif
