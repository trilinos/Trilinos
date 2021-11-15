// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/PerceptUtils.hpp>
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <Shards_CellTopologyData.h>

namespace percept {

std::vector<std::string> entity_rank_names_and_family_tree()
{
  std::vector<std::string> rank_names = stk::mesh::entity_rank_names();
  rank_names.push_back("FAMILY_TREE");
  return rank_names;
}

void computeCentroid(stk::mesh::Entity elem, double centroid[3], const stk::mesh::FieldBase & coord_field)
{
  int spaceDim = coord_field.mesh_meta_data().spatial_dimension();
  
  const stk::mesh::BulkData& mesh = coord_field.get_mesh();
  
  const stk::mesh::Entity * nodes = mesh.begin(elem, stk::topology::NODE_RANK);
  const unsigned        num_nodes = mesh.num_connectivity(elem, stk::topology::NODE_RANK);

  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  double nnode = num_nodes;

  for (unsigned ii=0; ii < num_nodes; ++ii)
    {
      double * node_coord_data_0 = (double*)stk::mesh::field_data(coord_field, nodes[ii]);
      for (int iSpaceDimOrd = 0; iSpaceDimOrd < spaceDim; iSpaceDimOrd++)
        {
          centroid[iSpaceDimOrd] += node_coord_data_0[iSpaceDimOrd] / nnode;
        }
    }
}

double volume(stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data)
{
  VolumeUtil jacA;
  shards::CellTopology cell_topo(cell_topo_data);
  double volScale = jacA.getJacobianToVolumeScale(cell_topo);
  
  double jacobian = 0.0;
  jacA(jacobian, element, coord_field, cell_topo_data);
  double cellVol = jacobian*volScale;
  return cellVol;
}

stk::diag::Timer& rootTimerStructured() { 
  static stk::diag::TimerSet s_timerSet(sierra::Diag::TIMER_ALL);
  static stk::diag::Timer s_timer = stk::diag::createRootTimer("Structured", s_timerSet);
  return s_timer;
}

void printTimersTableStructured() {
  std::ostringstream str;
  rootTimerStructured().stop();
  stk::diag::printTimersTable(str, rootTimerStructured(), stk::diag::METRICS_ALL, false);
  
  if (0==stk::parallel_machine_rank(MPI_COMM_WORLD))
    std::cout << str.str() << std::endl;
}

double MegaByte(const size_t x) {return ((double)x/1024.0/1024.0);}

void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg, size_t& hwm_sum)
{
  size_t now=0, hwm=0;
  stk::get_memory_usage(now, hwm);
  stk::all_reduce_max(comm, &hwm, &hwm_max, 1);
  stk::all_reduce_min(comm, &hwm, &hwm_min, 1);
  stk::all_reduce_sum(comm, &hwm, &hwm_sum, 1);
  int num_procs = stk::parallel_machine_size(comm);
  hwm_avg = hwm_sum / size_t(num_procs);
}

void get_memory_now_across_processors(MPI_Comm comm, size_t& now_max, size_t& now_min, size_t& now_avg, size_t& now_sum)
{
  size_t now=0, hwm=0;
  stk::get_memory_usage(now, hwm);
  stk::all_reduce_max(comm, &now, &now_max, 1);
  stk::all_reduce_min(comm, &now, &now_min, 1);
  stk::all_reduce_sum(comm, &now, &now_sum, 1);
  int num_procs = stk::parallel_machine_size(comm);
  now_avg = now_sum / size_t(num_procs);
}

std::string
print_memory_high_water_mark(MPI_Comm comm)
{
  size_t hwm_max=0, hwm_min=0, hwm_avg=0, hwm_sum=0;
  get_memory_high_water_mark_across_processors(comm, hwm_max, hwm_min, hwm_avg, hwm_sum);
  return stk::human_bytes(hwm_sum) + " [hwm_tot] " + stk::human_bytes(hwm_max) + " [hwm_max]";
}

std::string
print_memory_now(MPI_Comm comm)
{
  size_t now_max=0, now_min=0, now_avg=0, now_sum=0;
  get_memory_now_across_processors(comm, now_max, now_min, now_avg, now_sum);
  return stk::human_bytes(now_sum) + " [now_tot] " + stk::human_bytes(now_max) + " [now_max]";
}

std::string
print_memory_both(MPI_Comm comm)
{
  return print_memory_now(comm)+ " " + print_memory_high_water_mark(comm);
}

void MemoryInfo::
fill_memory_info_across_processors(MPI_Comm comm)
{
    size_t local_now=0, local_hwm=0;
    stk::get_memory_usage(local_now, local_hwm);

    // FIXME - make this a single reduction op
    stk::all_reduce_max(comm, &local_now, &now[MAX], 1);
    stk::all_reduce_min(comm, &local_now, &now[MIN], 1);
    stk::all_reduce_sum(comm, &local_now, &now[SUM], 1);
    now[AVG] = now[SUM] / size_t(stk::parallel_machine_size(comm));

    //stk::all_reduce_max(comm, &local_hwm, &hwm[MAX], 1);
    //stk::all_reduce_min(comm, &local_hwm, &hwm[MIN], 1);
    //stk::all_reduce_sum(comm, &local_hwm, &hwm[SUM], 1);
    //hwm[AVG] = hwm[SUM] / size_t(stk::parallel_machine_size(comm));
}

void MemoryInfo::
write_memory_log_header(const std::string& logfile_name)
{
    std::ofstream logfile(logfile_name.c_str());
    logfile << std::setw(12) << "#Stages \\ MB";
    for (unsigned int i=0; i<INVALID_SLOT; i++)
        logfile << std::setw(col_width) << header[i];
    //for (unsigned int i=0; i<INVALID_SLOT; i++)
    //    logfile << std::setw(col_width) << header[i];
    logfile << std::endl;
}

void MemoryInfo::
write_memory_log_data  (const std::string& logfile_name, ExecutionStage stage)
{
    static const double scale = 1.0 / (1024. * 1024.); // Memory in MB.

    std::ofstream logfile(logfile_name.c_str(), std::fstream::app);
    logfile << std::setw(12) << stages[stage] << " ";
    for (unsigned int i=0; i<INVALID_SLOT; i++)
        logfile << std::fixed << std::setprecision(2) << std::setw(col_width) << now[i] * scale;

    //for (unsigned int i=0; i<INVALID_SLOT; i++)
    //    logfile << std::fixed << std::setprecision(2) << std::setw(col_width) << hwm[i] * scale;
    logfile << std::endl;
}

void write_memory_logfile(MPI_Comm comm, ExecutionStage stage, const std::string& memory_logfile_name)
{
    static MemoryInfo meminfo;

    if (memory_logfile_name.size()==0) return;

    if (stage==INITIALIZE)
        meminfo.write_memory_log_header(memory_logfile_name);

    meminfo.fill_memory_info_across_processors(comm);

    if (stk::parallel_machine_rank(comm) == 0)
        meminfo.write_memory_log_data(memory_logfile_name, stage);
}

void convert_stk_topology_to_ioss_name(
        const stk::topology stk_topo,
        std::string& ioss_topo)
{
    if      (stk_topo == stk::topology::NODE)
      ioss_topo = "sphere";
    else if (stk_topo == stk::topology::LINE_2)
      ioss_topo = "edge2";
    else if (stk_topo == stk::topology::LINE_3)
      ioss_topo = "edge3";
    // 2D
    else if (stk_topo == stk::topology::TRI_3)
      ioss_topo = "tri3";
    else if (stk_topo == stk::topology::TRI_6)
      ioss_topo = "tri6";
    else if (stk_topo == stk::topology::QUAD_4)
      ioss_topo = "quad4";
    else if (stk_topo == stk::topology::QUAD_8)
      ioss_topo = "quad8";
    else if (stk_topo == stk::topology::QUAD_9)
      ioss_topo = "quad9";
    // 3D
    else if (stk_topo == stk::topology::HEX_8)
      ioss_topo = "hex8";
    else if (stk_topo == stk::topology::HEX_20)
      ioss_topo = "hex20";
    else if (stk_topo == stk::topology::HEX_27)
      ioss_topo = "hex27";
    else if (stk_topo == stk::topology::TET_4)
      ioss_topo = "tetra4";
    else if (stk_topo == stk::topology::TET_10)
      ioss_topo = "tetra10";
    else if (stk_topo == stk::topology::WEDGE_6)
      ioss_topo = "wedge6";
    else if (stk_topo == stk::topology::WEDGE_15)
      ioss_topo = "wedge15";
    else if (stk_topo == stk::topology::WEDGE_18)
      ioss_topo = "wedge18";
    else if (stk_topo == stk::topology::PYRAMID_5)
      ioss_topo = "pyramid5";
    else if (stk_topo == stk::topology::SHELL_QUADRILATERAL_4)
      ioss_topo = "shell4";
    else if (stk_topo == stk::topology::SHELL_QUADRILATERAL_8)
      ioss_topo = "shell8";
    else if (stk_topo == stk::topology::SHELL_QUADRILATERAL_9)
      ioss_topo = "shell9";
    else if (stk_topo == stk::topology::SHELL_TRI_3)
      ioss_topo = "trishell3";
    else if (stk_topo == stk::topology::SHELL_TRI_6)
      ioss_topo = "trishell6";
    else {
      std::stringstream errMsg;
      errMsg << "Percept convert_stk_topology_to_ioss_name not implemented for topology "
              << stk_topo << std::endl;
      throw std::runtime_error(errMsg.str().c_str());
    }
}

}
