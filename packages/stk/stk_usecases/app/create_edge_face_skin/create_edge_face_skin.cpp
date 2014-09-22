/*------------------------------------------------------------------------*/
/*  Copyright (c) 2013, Sandia Corporation.
/*  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*  the U.S. Governement retains certain rights in this software.
/*  
/*  Redistribution and use in source and binary forms, with or without
/*  modification, are permitted provided that the following conditions are
/*  met:
/*  
/*      * Redistributions of source code must retain the above copyright
/*        notice, this list of conditions and the following disclaimer.
/*  
/*      * Redistributions in binary form must reproduce the above
/*        copyright notice, this list of conditions and the following
/*        disclaimer in the documentation and/or other materials provided
/*        with the distribution.
/*  
/*      * Neither the name of Sandia Corporation nor the names of its
/*        contributors may be used to endorse or promote products derived
/*        from this software without specific prior written permission.
/*  
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*  
/*------------------------------------------------------------------------*/

#include <string>
#include <iostream>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/environment/perf_util.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>

#include <stk_util/util/human_bytes.hpp>
#include <stk_util/util/ParameterList.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/environment/ProgramOptions.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace {
  void provide_entity_count(stk::mesh::BulkData &bulk, int proc);

  // Do the actual reading of the mesh database and
  // creation and population of the MetaData and BulkData.
  void mesh_read_write(const std::string &type,
		       const std::string &working_directory,
		       const std::string &filename,
		       stk::io::StkMeshIoBroker &mesh_data,
		       bool create_edges,
		       bool create_faces,
		       bool create_skin)
  {
    std::string file = working_directory;
    file += filename;

    size_t input_index = mesh_data.add_mesh_database(file, type, stk::io::READ_MESH);
    mesh_data.set_active_mesh(input_index);

    double start_time = stk::cpu_time();
    mesh_data.create_input_mesh();

    // Unused if !create_skin.
    stk::mesh::Part & skin_part = mesh_data.meta_data().declare_part("skin_part");

    mesh_data.populate_bulk_data();
    double mesh_create_time = stk::cpu_time() - start_time;

    start_time = stk::cpu_time();
    if (create_skin) {
      stk::mesh::PartVector add_parts(1,&skin_part);
      stk::mesh::skin_mesh(mesh_data.bulk_data(), add_parts);
    }
    double create_skin_time = stk::cpu_time() - start_time;

    start_time = stk::cpu_time();
    if (create_faces) {
      stk::mesh::create_faces(mesh_data.bulk_data());
    }
    double create_faces_time = stk::cpu_time() - start_time;

    start_time = stk::cpu_time();
    if (create_edges) {
      stk::mesh::create_edges(mesh_data.bulk_data());
    }
    double create_edges_time = stk::cpu_time() - start_time;

    double total_time = mesh_create_time + create_skin_time + create_faces_time + create_edges_time;

    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(mesh_data.bulk_data(), mesh_counts);

    int proc = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);

    double all_face_time = 0.0;
    MPI_Allreduce(&create_faces_time, &all_face_time, 1, MPI_DOUBLE, MPI_SUM, mesh_data.bulk_data().parallel());

    double all_edge_time = 0.0;
    MPI_Allreduce(&create_edges_time, &all_edge_time, 1, MPI_DOUBLE, MPI_SUM, mesh_data.bulk_data().parallel());

    provide_entity_count(mesh_data.bulk_data(), proc);
    if (proc == 0) {
      std::cout<<"num nodes: "<<mesh_counts[stk::topology::NODE_RANK]<<std::endl;
      std::cout<<"num edges: "<<mesh_counts[stk::topology::EDGE_RANK]<<std::endl;
      std::cout<<"num faces: "<<mesh_counts[stk::topology::FACE_RANK]<<std::endl;
      std::cout<<"num elems: "<<mesh_counts[stk::topology::ELEM_RANK]<<"\n"<<std::endl;
      if (create_edges) {
	std::cout<< "num edges/second: " << mesh_counts[stk::topology::EDGE_RANK]/all_edge_time << "\t" << all_edge_time << std::endl;
      }
      if (create_faces) {
	std::cout<< "num faces/second: " << mesh_counts[stk::topology::FACE_RANK]/all_face_time << "\t" << all_face_time << std::endl;
      }
    }

    size_t now=0, hwm=0;
    stk::get_memory_usage(now,hwm);
    size_t global_hwm=0;
    MPI_Allreduce(&hwm, &global_hwm, 1, MPI_LONG_LONG, MPI_MAX, mesh_data.bulk_data().parallel());

    if (proc == 0) {
      const double timers[] = {mesh_create_time, create_edges_time, create_faces_time, create_skin_time, total_time};
      const char* timer_names[] = {"Create mesh", "Create edges", "Create faces", "Create skin", "Total time"};
      const int NUM_TIMERS = sizeof(timers)/sizeof(timers[0]);
      
      stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);
      std::cout<<"Global HWM: "<<stk::human_bytes(global_hwm)<<std::endl;
    }

    stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
  }

  void driver(stk::ParallelMachine  comm,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type,
	      const std::string &decomp_method,
	      bool create_edges,
	      bool create_faces,
	      bool create_skin,
	      int id_integer_size)
  {
    stk::io::StkMeshIoBroker mesh_data(comm);
    mesh_data.property_add(Ioss::Property("INTEGER_SIZE_API", id_integer_size));

    if (!decomp_method.empty()) {
      mesh_data.property_add(Ioss::Property("DECOMPOSITION_METHOD", decomp_method));
    }

    mesh_read_write(type, working_directory, filename, mesh_data, create_edges, create_faces, create_skin);
  }

  void provide_entity_count(stk::mesh::BulkData &bulk, int proc) {
    std::vector<size_t> counts;
    std::vector<size_t> minCounts;
    std::vector<size_t> maxCounts;
    stk::mesh::comm_mesh_counts(bulk, counts, minCounts, maxCounts);

    if (proc == 0) {
      std::cout << "===========================" << std::endl;
      std::cout << "nodes,    " << std::setw(10) << counts[0] << ", min/max: " << minCounts[0] << "/" << maxCounts[0] << std::endl
		<< "edges,    " << std::setw(10) << counts[1] << ", min/max: " << minCounts[1] << "/" << maxCounts[1] << std::endl
		<< "faces,    " << std::setw(10) << counts[2] << ", min/max: " << minCounts[2] << "/" << maxCounts[2] << std::endl
		<< "elements, " << std::setw(10) << counts[3] << ", min/max: " << minCounts[3] << "/" << maxCounts[3] << std::endl
		<< "===========================" << std::endl;
    }
  }
}

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  std::string working_directory = "";
  std::string decomp_method = "";
  std::string mesh = "";
  std::string type = "exodusii";
  int id_integer_size = 4;
  bool create_faces = false;
  bool create_edges = false;
  bool create_skin  = false;
  //----------------------------------
  // Process the broadcast command line arguments
  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
  stk::BroadcastArg b_arg(comm, argc, argv);
  bopt::options_description desc("options");

  desc.add_options()
    ("help,h", "produce help message")
    ("directory,d",   bopt::value<std::string>(&working_directory),
     "working directory with trailing '/'" )
    ("decomposition,D", bopt::value<std::string>(&decomp_method),
     "decomposition method.  One of: linear, rcb, rib, hsfc, block, cyclic, random, kway, geom_kway, metis_sfc" )
    ("faces", "create all faces" )
    ("edges", "create all edges" )
    ("skin",  "create all boundary faces (skin of the model)" )
    ("id_integer_size", bopt::value<int>(&id_integer_size), "use 4 or 8-byte integers for ids" )
    ("mesh",          bopt::value<std::string>(&mesh),
     "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Can also specify a filename." );

  stk::get_options_description().add(desc);

  bopt::variables_map &vm = stk::get_variables_map();
  try {
    bopt::store(bopt::parse_command_line(b_arg.m_argc, b_arg.m_argv, desc), vm);
    bopt::notify(vm);
  }
  catch (std::exception & /* x */) {
    std::exit(1);
  }

  if (mesh.empty()) {
    std::cerr << "\nERROR: The --mesh option is required\n";
    std::cerr << "\nApplication " << desc << "\n";
    std::exit(EXIT_FAILURE);
  }

  if (vm.count("help")) {
    std::cout << "Usage: " << argv[0] << " " << desc << "\n";
    std::exit(EXIT_SUCCESS);
  }
  
  if (vm.count("faces")) {create_faces = true;}
  if (vm.count("edges")) {create_edges = true;}
  if (vm.count("skin"))  {create_skin  = true;}

  type = "exodusii";
  if (strncasecmp("gen:", mesh.c_str(), 4) == 0) {
    mesh = mesh.substr(4, mesh.size());
    type = "generated";
  }
  if (strncasecmp("dof:", mesh.c_str(), 4) == 0) {
    mesh = mesh.substr(4, mesh.size());
    type = "dof";
  }

  driver(comm, working_directory, mesh, type, decomp_method,
	 create_edges, create_faces, create_skin, id_integer_size);

  stk::parallel_machine_finalize();
  return 0;
}

