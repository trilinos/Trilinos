// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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


#include <Ioss_SubSystem.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Transform.h>
#include <Ioss_FileInfo.h>
#include <Ioss_Utils.h>
#include <Ioss_ParallelUtils.h>
#include <assert.h>
#include <init/Ionit_Initializer.h>
#include <stddef.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <sys/times.h>

#include "shell_interface.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifndef NO_XDMF_SUPPORT
#include <xdmf/Ioxf_Initializer.h>
#endif

#define OUTPUT if (rank == 0) std::cerr

// ========================================================================

namespace {

  double timer()
  {
#ifdef HAVE_MPI
    return MPI_Wtime();
#else
    static double ticks_per_second = 0.0;
    struct tms time_buf;

    if (ticks_per_second == 0.0) {
      ticks_per_second = double(sysconf(_SC_CLK_TCK));
    }

    clock_t ctime = times(&time_buf);
    double time = double(ctime) / ticks_per_second;
    return time;
#endif
  }

  size_t MAX(size_t a, size_t b)
  {
    return b ^ ((a ^ b) & -(a>b));
  }

  template<typename T>  struct remove_pointer     { typedef T type; };
  template<typename T>  struct remove_pointer<T*> { typedef T type; };
  
  // Data space shared by most field input/output routines...
  std::vector<char> data;
  size_t max_field_size = 0;
  int rank = 0;
  int64_t data_read  = 0;
  int64_t data_write = 0;
  double time_read = 0.0;
  double time_write = 0.0;

  void show_usage(const std::string &prog);
  void show_step(int istep, double time);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_edgeblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_faceblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_edgesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_facesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_elemsets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_coordinate_frames(Ioss::Region &region, Ioss::Region &output_region, bool debug);

  template <typename T>
  void transfer_fields(const std::vector<T*>& entities,
		       Ioss::Region &output_region,
		       Ioss::Field::RoleType role,
		       const IOShell::Interface &interface);

  void transfer_fields(Ioss::GroupingEntity *ige,
		       Ioss::GroupingEntity *oge,
		       Ioss::Field::RoleType role,
		       const std::string &prefix = "");

  template <typename T>
  void transfer_field_data(const std::vector<T*>& entities,
			   Ioss::Region &output_region,
			   Ioss::Field::RoleType role,
			   const IOShell::Interface &interface);

  void transfer_field_data(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge,
			   Ioss::Field::RoleType role,
			   const std::string &prefix = "");

  void transfer_properties(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge);

  void transfer_qa_info(Ioss::Region &in, Ioss::Region &out);
  
  void transform_fields(Ioss::GroupingEntity *ige,
			Ioss::GroupingEntity *oge,
			Ioss::Field::RoleType role);

  void transform_field_data(Ioss::GroupingEntity *ige,
			    Ioss::GroupingEntity *oge,
			    Ioss::Field::RoleType role);
  void transfer_field_data_internal(Ioss::GroupingEntity *ige,
				    Ioss::GroupingEntity *oge,
				    const std::string &field_name);

  void file_copy(IOShell::Interface &interface);

  template <typename INT>
  void set_owned_node_count(Ioss::Region &region, int processor, INT dummy);
}
// ========================================================================

namespace {
  std::string codename;
  std::string version = "4.6";
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  IOShell::Interface interface;
  bool success = interface.parse_options(argc, argv);
  if (!success) {
    exit(EXIT_FAILURE);
  }
  
  Ioss::Init::Initializer io;
#ifndef NO_XDMF_SUPPORT
  Ioxf::Initializer ioxf;
#endif

  std::string in_file   = interface.inputFile[0];
  std::string out_file  = interface.outputFile;

  OUTPUT << "Input:    '" << in_file  << "', Type: " << interface.inFiletype  << '\n';
  OUTPUT << "Output:   '" << out_file << "', Type: " << interface.outFiletype << '\n';
  OUTPUT << '\n';
  
  double begin = timer();
  file_copy(interface);
  double end = timer();

#ifdef HAVE_MPI
  // Get total data read/written over all processors, and the max time..
  Ioss::ParallelUtils parallel(MPI_COMM_WORLD);

  // Combine these some time...
  time_read  = parallel.global_minmax(time_read,  Ioss::ParallelUtils::DO_MAX);
  time_write = parallel.global_minmax(time_write, Ioss::ParallelUtils::DO_MAX);
  data_read  = parallel.global_minmax(data_read,  Ioss::ParallelUtils::DO_SUM);
  data_write = parallel.global_minmax(data_write, Ioss::ParallelUtils::DO_SUM);
#endif

  if (interface.statistics) {
    OUTPUT << "\n\tElapsed time = " << end - begin << " seconds."
	   << " (read = " << time_read << ", write = " << time_write << ")\n"
	   << "\tTotal data read  = " << data_read  << " bytes; "
	   << data_read  / time_read  <<" bytes/second.\n"
	   << "\tTotal data write = " << data_write << " bytes; "
	   << data_write / time_write << " bytes/second.\n";
  }

  OUTPUT << "\n" << codename << " execution successful.\n";
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
  void file_copy(IOShell::Interface &interface)
  {
    Ioss::PropertyManager properties;
    if (interface.ints_64_bit) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB",  8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interface.reals_32_bit) {
      properties.add(Ioss::Property("REAL_SIZE_DB",  4));
    }

    if (interface.compression_level > 0 || interface.shuffle) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interface.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", interface.shuffle));
    }
      
    if (interface.compose_output != "none") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
      if (interface.compose_output != "default")
	properties.add(Ioss::Property("PARALLEL_IO_MODE", interface.compose_output));
    }

    if (interface.netcdf4) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }
    
    if (interface.inputFile.size() > 1) {
      properties.add(Ioss::Property("ENABLE_FILE_GROUPS", true));
    }

    if (interface.debug)
      properties.add(Ioss::Property("LOGGING", 1));

    Ioss::PropertyManager properties_in;
    if (!interface.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interface.decomp_method));
    }
    //========================================================================
    // OUTPUT ...
    //========================================================================
    for (size_t input_db = 0; input_db < interface.inputFile.size(); input_db++) {

      Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(interface.outFiletype, interface.outputFile, 
						      Ioss::WRITE_RESTART, (MPI_Comm)MPI_COMM_WORLD, properties);
      if (dbo == NULL || !dbo->ok(true)) {
	std::exit(EXIT_FAILURE);
      }

      // NOTE: 'output_region' owns 'dbo' pointer at this time
      Ioss::Region output_region(dbo, "region_2");
      // Set the qa information...
      output_region.property_add(Ioss::Property(std::string("code_name"), codename));
      output_region.property_add(Ioss::Property(std::string("code_version"), version));

      std::string inpfile = interface.inputFile[input_db];
      
      Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(interface.inFiletype, inpfile,
						      Ioss::READ_MODEL, (MPI_Comm)MPI_COMM_WORLD, properties);
      if (dbi == NULL || !dbi->ok(true)) {
	std::exit(EXIT_FAILURE);
      }

      dbi->set_surface_split_type(Ioss::int_to_surface_split(interface.surface_split_type));
      dbi->set_field_separator(interface.fieldSuffixSeparator);
      if (interface.ints_64_bit)
	dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    
      if (!interface.groupName.empty()) {
	bool success = dbi->open_group(interface.groupName);
	if (!success) {
	  OUTPUT << "ERROR: Unable to open group '" << interface.groupName
		 << "' in file '" << inpfile << "\n";
	  return;
	}
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region region(dbi, "region_1");
    
      if (interface.inputFile.size() > 1) {
	properties.add(Ioss::Property("APPEND_OUTPUT",Ioss::DB_APPEND_GROUP)); 

	if (input_db > 0) {
	  // Putting each file into its own output group...
	  // The name of the group will be the basename portion of the filename...
	  Ioss::FileInfo file(interface.inputFile[input_db]);
	  dbo->create_subgroup(file.tailname());
	}
      }
    
      if (interface.debug) OUTPUT << "DEFINING MODEL ... \n";
      if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
	OUTPUT << "ERROR: Could not put output region into define model state\n";
	std::exit(EXIT_FAILURE);
      }

      // Get all properties of input database...
      transfer_properties(&region, &output_region);
      transfer_qa_info(region, output_region);

      transfer_nodeblock(region, output_region, interface.debug);

#ifdef HAVE_MPI
      // This also assumes that the node order and count is the same for input
      // and output regions... (This is checked during nodeset output)
      if (output_region.get_database()->needs_shared_node_information()) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (interface.ints_64_bit)
	  set_owned_node_count(region, rank, (int64_t)0);
	else
	  set_owned_node_count(region, rank, (int)0);
      }
#endif

      transfer_edgeblocks(region, output_region, interface.debug);
      transfer_faceblocks(region, output_region, interface.debug);
      transfer_elementblocks(region, output_region, interface.debug);

      transfer_nodesets(region, output_region, interface.debug);
      transfer_edgesets(region, output_region, interface.debug);
      transfer_facesets(region, output_region, interface.debug);
      transfer_elemsets(region, output_region, interface.debug);

      transfer_sidesets(region, output_region, interface.debug);
      transfer_commsets(region, output_region, interface.debug);

      transfer_coordinate_frames(region, output_region, interface.debug);

      if (interface.debug) OUTPUT << "END STATE_DEFINE_MODEL... " << '\n';
      output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

      OUTPUT << "Maximum Field size = " << max_field_size << " bytes.\n";
      data.resize(max_field_size);
      OUTPUT << "Resize finished...\n";
    
      if (interface.debug) OUTPUT << "TRANSFERRING MESH FIELD DATA ... " << '\n';
      // Model defined, now fill in the model data...
      output_region.begin_mode(Ioss::STATE_MODEL);

      // Transfer MESH field_data from input to output...
      transfer_field_data(region.get_node_blocks(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_node_blocks(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_edge_blocks(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_edge_blocks(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_face_blocks(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_face_blocks(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_element_blocks(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_element_blocks(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_nodesets(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_nodesets(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_edgesets(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_edgesets(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_facesets(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_facesets(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_elementsets(), output_region, Ioss::Field::MESH,      interface);
      transfer_field_data(region.get_elementsets(), output_region, Ioss::Field::ATTRIBUTE, interface);

      transfer_field_data(region.get_commsets(), output_region, Ioss::Field::MESH,          interface);
      transfer_field_data(region.get_commsets(), output_region, Ioss::Field::ATTRIBUTE,     interface);
      transfer_field_data(region.get_commsets(), output_region, Ioss::Field::COMMUNICATION, interface);

      // Side Sets
      {
	Ioss::SideSetContainer fss = region.get_sidesets();
	Ioss::SideSetContainer::const_iterator I = fss.begin();
	while (I != fss.end()) {
	  std::string name     = (*I)->name();
	  if (interface.debug) OUTPUT << name << ", ";

	  // Find matching output sideset
	  Ioss::SideSet *ofs = output_region.get_sideset(name);

	  if (ofs != NULL) {
	    transfer_field_data(*I, ofs, Ioss::Field::MESH);
	    transfer_field_data(*I, ofs, Ioss::Field::ATTRIBUTE);

	    Ioss::SideBlockContainer fbs = (*I)->get_side_blocks();
	    Ioss::SideBlockContainer::const_iterator J = fbs.begin();
	    while (J != fbs.end()) {

	      // Find matching output sideblock
	      std::string fbname = (*J)->name();
	      if (interface.debug) OUTPUT << fbname << ", ";
	      Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

	      if (ofb != NULL) {
		transfer_field_data(*J, ofb, Ioss::Field::MESH);
		transfer_field_data(*J, ofb, Ioss::Field::ATTRIBUTE);
	      }
	      ++J;
	    }
	  }
	  ++I;
	}
	if (interface.debug) OUTPUT << '\n';
      }
      if (interface.debug) OUTPUT << "END STATE_MODEL... " << '\n';
      output_region.end_mode(Ioss::STATE_MODEL);

      if (interface.debug) OUTPUT << "DEFINING TRANSIENT FIELDS ... " << '\n';
      if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
	if (!interface.debug) {
	  OUTPUT << "\n Number of time steps on database     =" << std::setw(12)
		 << region.get_property("state_count").get_int() << "\n\n";
	}

	output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

	// For each 'TRANSIENT' field in the node blocks and element
	// blocks, transfer to the output node and element blocks.
	transfer_fields(&region, &output_region, Ioss::Field::TRANSIENT);

	transfer_fields(region.get_node_blocks(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_fields(region.get_edge_blocks(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_fields(region.get_face_blocks(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_fields(region.get_element_blocks(), output_region, Ioss::Field::TRANSIENT, interface);

	transfer_fields(region.get_nodesets(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_fields(region.get_edgesets(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_fields(region.get_facesets(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_fields(region.get_elementsets(), output_region, Ioss::Field::TRANSIENT, interface);

	// Side Sets
	{
	  Ioss::SideSetContainer fss = region.get_sidesets();
	  Ioss::SideSetContainer::const_iterator I = fss.begin();
	  while (I != fss.end()) {
	    std::string name     = (*I)->name();
	    if (interface.debug) OUTPUT << name << ", ";

	    // Find matching output sideset
	    Ioss::SideSet *ofs = output_region.get_sideset(name);

	    if (ofs != NULL) {
	      transfer_fields(*I, ofs, Ioss::Field::TRANSIENT);

	      Ioss::SideBlockContainer fbs = (*I)->get_side_blocks();
	      Ioss::SideBlockContainer::const_iterator J = fbs.begin();
	      while (J != fbs.end()) {

		// Find matching output sideblock
		std::string fbname = (*J)->name();
		if (interface.debug) OUTPUT << fbname << ", ";
		Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

		if (ofb != NULL) {
		  transfer_fields(*J, ofb, Ioss::Field::TRANSIENT);
		}
		++J;
	      }
	    }
	    ++I;
	  }
	  if (interface.debug) OUTPUT << '\n';
	}
	if (interface.debug) OUTPUT << "END STATE_DEFINE_TRANSIENT... " << '\n';
	output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
      }

      if (interface.debug) OUTPUT << "TRANSFERRING TRANSIENT FIELDS ... " << '\n';
      output_region.begin_mode(Ioss::STATE_TRANSIENT);
      // Get the timesteps from the input database.  Step through them
      // and transfer fields to output database...

      int step_count = region.get_property("state_count").get_int();

      for (int istep = 1; istep <= step_count; istep++) {
	double time = region.get_state_time(istep);
	if (time < interface.minimum_time)
	  continue;
	if (interface.maximum_time != 0.0 && time > interface.maximum_time)
	  break;
      
	int ostep = output_region.add_state(time);
	show_step(istep, time);

	output_region.begin_state(ostep);
	region.begin_state(istep);

	transfer_field_data(&region, &output_region, Ioss::Field::TRANSIENT);

	transfer_field_data(region.get_node_blocks(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_field_data(region.get_edge_blocks(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_field_data(region.get_face_blocks(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_field_data(region.get_element_blocks(), output_region, Ioss::Field::TRANSIENT, interface);

	transfer_field_data(region.get_nodesets(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_field_data(region.get_edgesets(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_field_data(region.get_facesets(), output_region, Ioss::Field::TRANSIENT, interface);
	transfer_field_data(region.get_elementsets(), output_region, Ioss::Field::TRANSIENT, interface);

	// Side Sets
	{
	  Ioss::SideSetContainer fss = region.get_sidesets();
	  Ioss::SideSetContainer::const_iterator I = fss.begin();
	  while (I != fss.end()) {
	    std::string name     = (*I)->name();
	    if (interface.debug) OUTPUT << name << ", ";
	  
	    // Find matching output sideset
	    Ioss::SideSet *ofs = output_region.get_sideset(name);
	  
	    if (ofs != NULL) {
	      transfer_field_data(*I, ofs, Ioss::Field::TRANSIENT);
	    
	      Ioss::SideBlockContainer fbs = (*I)->get_side_blocks();
	      Ioss::SideBlockContainer::const_iterator J = fbs.begin();
	      while (J != fbs.end()) {
	      
		// Find matching output sideblock
		std::string fbname = (*J)->name();
		if (interface.debug) OUTPUT << fbname << ", ";
		Ioss::SideBlock *ofb = ofs->get_side_block(fbname);
	      
		if (ofb != NULL) {
		  transfer_field_data(*J, ofb, Ioss::Field::TRANSIENT);
		}
		++J;
	      }
	    }
	    ++I;
	  }
	}
	region.end_state(istep);
	output_region.end_state(ostep);
      }
      if (interface.debug) OUTPUT << "END STATE_TRANSIENT... " << '\n';
      output_region.end_mode(Ioss::STATE_TRANSIENT);
    }
  }
  
  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeBlockContainer    nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i = nbs.begin();
    size_t id = 1;
    while (i != nbs.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
      size_t num_nodes = (*i)->get_property("entity_count").get_int();
      size_t degree    = (*i)->get_property("component_degree").get_int();
      if (!debug) {
	OUTPUT << " Number of coordinates per node       =" << std::setw(12) << degree << "\n";
	OUTPUT << " Number of nodes                      =" << std::setw(12) << num_nodes << "\n";
      }

      Ioss::NodeBlock *nb = new Ioss::NodeBlock(output_region.get_database(), name, num_nodes, degree);
      output_region.add(nb);

      
      transfer_properties(*i, nb);

      if (output_region.get_database()->needs_shared_node_information()) {
	// If the "owning_processor" field exists on the input
	// nodeblock, transfer it and the "ids" field to the output
	// nodeblock at this time since it is used to determine
	// per-processor sizes of nodeblocks and nodesets.
	if ((*i)->field_exists("owning_processor")) {
	  size_t isize = (*i)->get_field("ids").get_size();
	  data.resize(isize);
	  double t1 = timer();
	  (*i)->get_field_data("ids", &data[0], isize);
	  double t2 = timer();
	  nb->put_field_data("ids", &data[0], isize);
	  time_write += timer()-t2;
	  time_read += t2-t1;
	  
	  data_read += isize;
	  data_write += isize;

	  isize = (*i)->get_field("owning_processor").get_size();
	  data.resize(isize);
	  t1 = timer();
	  (*i)->get_field_data("owning_processor", &data[0], isize);
	  t2 = timer();
	  nb->put_field_data("owning_processor", &data[0], isize);
	  time_write += timer()-t2;
	  time_read += t2-t1;
	  data_read += isize;
	  data_write += isize;
	}
      }

      transfer_fields(*i, nb, Ioss::Field::MESH);
      transfer_fields(*i, nb, Ioss::Field::ATTRIBUTE);
      ++i;
      ++id;
    }
    if (debug) OUTPUT << '\n';
  }

  template <typename T>
  void transfer_fields(const std::vector<T*>& entities, Ioss::Region &output_region, Ioss::Field::RoleType role,
		       const IOShell::Interface &interface)
  {
    typename std::vector<T*>::const_iterator i = entities.begin();
    while (i != entities.end()) {
      std::string name      = (*i)->name();
      if (interface.debug) OUTPUT << name << ", ";
      
      // Find the corresponding output node_block...
      Ioss::GroupingEntity *oeb = output_region.get_entity(name, (*i)->type());
      if (oeb != NULL) {
	transfer_fields(*i, oeb, role);
	if (interface.do_transform_fields)
	  transform_fields(*i, oeb, role);
      }
      ++i;
    }
    if (interface.debug) OUTPUT << '\n';
  }

  template <typename T>
  void transfer_field_data(const std::vector<T*>& entities, Ioss::Region &output_region, Ioss::Field::RoleType role,
			   const IOShell::Interface &interface)
  {
    typename std::vector<T*>::const_iterator i = entities.begin();
    while (i != entities.end()) {
      T *entity = *i; ++i;
      std::string name = entity->name();

      // Find the corresponding output block...
      Ioss::GroupingEntity *output = output_region.get_entity(name, entity->type());
      if (output != NULL) {
	transfer_field_data(entity, output, role);
	if (interface.do_transform_fields)
	  transform_field_data(entity, output, role);
      }
    }
  }

  template <typename T>
  void transfer_blocks(std::vector<T*>& blocks, Ioss::Region &output_region, bool debug)
  {
    if (!blocks.empty()) {
      typename std::vector<T*>::const_iterator i = blocks.begin();
      size_t total_entities = 0;
      while (i != blocks.end()) {
	std::string name      = (*i)->name();
	if (debug) OUTPUT << name << ", ";
	std::string type      = (*i)->get_property("topology_type").get_string();
	size_t    count  = (*i)->get_property("entity_count").get_int();
	total_entities += count;
	
	T* block = new T(output_region.get_database(), name, type, count);
	output_region.add(block);
	transfer_properties(*i, block);
	transfer_fields(*i, block, Ioss::Field::MESH);
	transfer_fields(*i, block, Ioss::Field::ATTRIBUTE);
	++i;
      }
      if (!debug) {
	OUTPUT << " Number of " << std::setw(14) << (*blocks.begin())->type_string() << "s            ="
	       << std::setw(12) << blocks.size() << "\t"
	       << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      } else {
	OUTPUT << '\n';
      }
    }
  }

  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    transfer_blocks(ebs, output_region, debug);
  }

  void transfer_edgeblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::EdgeBlockContainer ebs = region.get_edge_blocks();
    transfer_blocks(ebs, output_region, debug);
  }

  void transfer_faceblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::FaceBlockContainer ebs = region.get_face_blocks();
    transfer_blocks(ebs, output_region, debug);
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::SideSetContainer      fss = region.get_sidesets();
    Ioss::SideSetContainer::const_iterator i = fss.begin();
    size_t total_sides = 0;
    while (i != fss.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
      Ioss::SideSet *surf = new Ioss::SideSet(output_region.get_database(), name);

      Ioss::SideBlockContainer fbs = (*i)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator j = fbs.begin();
      while (j != fbs.end()) {
	std::string fbname    = (*j)->name();
	if (debug) OUTPUT << fbname << ", ";
	std::string fbtype    = (*j)->get_property("topology_type").get_string();
	std::string partype   = (*j)->get_property("parent_topology_type").get_string();
	size_t    num_side  = (*j)->get_property("entity_count").get_int();
	total_sides += num_side;

	Ioss::SideBlock *block = new Ioss::SideBlock(output_region.get_database(), fbname, fbtype,
						     partype, num_side);
	surf->add(block);
	transfer_properties(*j, block);
	transfer_fields(*j, block, Ioss::Field::MESH);
	transfer_fields(*j, block, Ioss::Field::ATTRIBUTE);
	++j;
      }
      transfer_properties(*i, surf);
      transfer_fields(*i, surf, Ioss::Field::MESH);
      transfer_fields(*i, surf, Ioss::Field::ATTRIBUTE);
      output_region.add(surf);
      ++i;
    }
    if (!debug) {
      OUTPUT << " Number of        SideSets            =" << std::setw(12) << fss.size() << "\t"
	     << "Number of element sides =" << std::setw(12) << total_sides << "\n";
    } else {
      OUTPUT << '\n';
    }
  }

  template <typename T>
  void transfer_sets(std::vector<T*> &sets, Ioss::Region &output_region, bool debug)
  {
    if (!sets.empty()) {
      typename std::vector<T*>::const_iterator i = sets.begin();
      size_t total_entities = 0;
      while (i != sets.end()) {
	std::string name      = (*i)->name();
	if (debug) OUTPUT << name << ", ";
	size_t    count     = (*i)->get_property("entity_count").get_int();
	total_entities += count;
	T* set = new T(output_region.get_database(), name, count);
	output_region.add(set);
	transfer_properties(*i, set);
	transfer_fields(*i, set, Ioss::Field::MESH);
	transfer_fields(*i, set, Ioss::Field::ATTRIBUTE);
	++i;
      }

      if (!debug) {
	OUTPUT << " Number of " << std::setw(14) << (*sets.begin())->type_string() << "s            ="
	       << std::setw(12) << sets.size() << "\t"
	       << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      } else {
	OUTPUT << '\n';
      }

      
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeSetContainer      nss = region.get_nodesets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_edgesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::EdgeSetContainer      nss = region.get_edgesets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_facesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::FaceSetContainer      nss = region.get_facesets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_elemsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::ElementSetContainer      nss = region.get_elementsets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::CommSetContainer      css = region.get_commsets();
    Ioss::CommSetContainer::const_iterator i = css.begin();
    while (i != css.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
      std::string type      = (*i)->get_property("entity_type").get_string();
      size_t    count     = (*i)->get_property("entity_count").get_int();
      Ioss::CommSet *cs = new Ioss::CommSet(output_region.get_database(), name, type, count);
      output_region.add(cs);
      transfer_properties(*i, cs);
      transfer_fields(*i, cs, Ioss::Field::MESH);
      transfer_fields(*i, cs, Ioss::Field::ATTRIBUTE);
      transfer_fields(*i, cs, Ioss::Field::COMMUNICATION);
      ++i;
    }
    if (debug) OUTPUT << '\n';
  }

  void transfer_coordinate_frames(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::CoordinateFrameContainer      cf = region.get_coordinate_frames();
    Ioss::CoordinateFrameContainer::const_iterator i = cf.begin();
    while (i != cf.end()) {
      output_region.add(*i);
      ++i;
    }
    if (debug) OUTPUT << '\n';
  }

  void transfer_fields(Ioss::GroupingEntity *ige,
		       Ioss::GroupingEntity *oge,
		       Ioss::Field::RoleType role,
		       const std::string &prefix)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output
    // database...  If a prefix is specified, only transfer fields
    // whose names begin with the prefix
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      Ioss::Field field = ige->get_field(field_name);
      max_field_size = MAX(max_field_size, field.get_size());
      if (field_name != "ids" && !oge->field_exists(field_name) &&
	  (prefix.length() == 0 || std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
	// If the field does not already exist, add it to the output node block
	oge->field_add(field);
      }
    }
  }

  void transform_fields(Ioss::GroupingEntity *ige,
			Ioss::GroupingEntity *oge,
			Ioss::Field::RoleType role)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output database...
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      std::string out_field_name = field_name + "_mag";
      if (!oge->field_exists(out_field_name)) {
	// If the field does not already exist, add it to the output node block
	Ioss::Field field = ige->get_field(field_name);
	Ioss::Field tr_field(out_field_name, field.get_type(), field.raw_storage(),
			     field.get_role(), field.raw_count());

	Ioss::Transform* transform = Iotr::Factory::create("vector magnitude");
	assert(transform != NULL);
	tr_field.add_transform(transform);

	Ioss::Transform* max_transform = Iotr::Factory::create("absolute_maximum");
	assert(max_transform != NULL);
	tr_field.add_transform(max_transform);

	oge->field_add(tr_field);
      }
    }
  }

  void transform_field_data(Ioss::GroupingEntity *ige,
			    Ioss::GroupingEntity *oge,
			    Ioss::Field::RoleType role)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    ige->field_describe(role, &state_fields);
    // Iterate through mesh description fields and transfer to
    // output database...
    Ioss::NameList::const_iterator IF;

    for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
      std::string field_name = *IF;
      std::string out_field_name = field_name + "_mag";

      assert(oge->field_exists(out_field_name));

      size_t isize = ige->get_field(field_name).get_size();
      size_t osize = oge->get_field(out_field_name).get_size();
      assert (isize == osize);

      assert(data.size() >= isize);
      data_read += isize;
      data_write += isize;
      double t1 = timer();
      ige->get_field_data(field_name,     &data[0], isize);
      double t2 = timer();
      oge->put_field_data(out_field_name, &data[0], osize);
      time_write += timer() - t2;
      time_read += t2-t1;
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge,
			   Ioss::Field::RoleType role,
			   const std::string &prefix)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    Ioss::NameList::const_iterator IF;
    ige->field_describe(role, &state_fields);

    // Complication here is that if the 'role' is 'Ioss::Field::MESH',
    // then the 'ids' field must be transferred first...
    if (role == Ioss::Field::MESH) {
      for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
	std::string field_name = *IF;
	assert(oge->field_exists(field_name));
	if (field_name == "ids") {
	  transfer_field_data_internal(ige, oge, field_name);
	  break;
	}
      }
    }

    for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
      std::string field_name = *IF;
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (field_name == "connectivity" && ige->type() != Ioss::ELEMENTBLOCK)
	continue;


      if (field_name != "ids" &&
	  (prefix.length() == 0 || std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
	assert(oge->field_exists(field_name));
	transfer_field_data_internal(ige, oge, field_name);
      }
    }
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige,
				    Ioss::GroupingEntity *oge,
				    const std::string &field_name)
  {

    size_t isize = ige->get_field(field_name).get_size();
    if (isize != oge->get_field(field_name).get_size()) {
      std::cerr << "Field: " << field_name << "\tIsize = " << isize << "\tOsize = " << oge->get_field(field_name).get_size() << "\n";
      assert(isize == oge->get_field(field_name).get_size());
    }

    if (field_name == "mesh_model_coordinates_x") return;
    if (field_name == "mesh_model_coordinates_y") return;
    if (field_name == "mesh_model_coordinates_z") return;
    if (field_name == "connectivity_raw") return;
    if (field_name == "element_side_raw") return;
    if (field_name == "ids_raw") return;
    if (field_name == "node_connectivity_status") return;
    if (field_name == "owning_processor") return;
    if (field_name == "entity_processor_raw") return;

    if (data.size() < isize) {
      std::cerr << "Field: " << field_name << "\tIsize = " << isize << "\tdata size = " << data.size() << "\n";
      data.resize(isize);
    }

    assert(data.size() >= isize);
    data_read += isize;
    data_write += isize;
    double t1 = timer();
    ige->get_field_data(field_name, &data[0], isize);
    double t2 = timer();
    oge->put_field_data(field_name, &data[0], isize);
    time_write += timer() - t2;
    time_read += t2-t1;
  }

  void transfer_qa_info(Ioss::Region &in, Ioss::Region &out)
  {
    out.add_information_records(in.get_information_records());

    const std::vector<std::string> &qa = in.get_qa_records();
    for (size_t i=0; i < qa.size(); i+=4) {
      out.add_qa_record(qa[i+0], qa[i+1], qa[i+2], qa[i+3]);
    }
  }

  void transfer_properties(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge)
  {
    Ioss::NameList names;
    ige->property_describe(&names);

    // Iterate through properties and transfer to output database...
    Ioss::NameList::const_iterator I;
    for (I = names.begin(); I != names.end(); ++I) {
      //       if (*I == "original_block_order" && oge->property_exists(*I))
      //        oge->property_erase("original_block_order");
      if (!oge->property_exists(*I))
	oge->property_add(ige->get_property(*I));
    }
  }

  void show_step(int istep, double time)
  {
    OUTPUT.setf(std::ios::scientific);
    OUTPUT.setf(std::ios::showpoint);
    OUTPUT << "     Time step " << std::setw(5) << istep
	   << " at time " << std::setprecision(5) << time << '\n';
  }

  template <typename INT>
  void set_owned_node_count(Ioss::Region &region, int my_processor, INT /*dummy*/)
  {
    Ioss::NodeBlock *nb = region.get_node_block("nodeblock_1");
    if (nb->field_exists("owning_processor")) {
      std::vector<INT> my_data;
      double t1 = timer();
      nb->get_field_data("owning_processor", my_data);
      time_read += timer()-t1;
      data_read += my_data.size();

      INT owned = std::count(my_data.begin(), my_data.end(), my_processor);
      nb->property_add(Ioss::Property("locally_owned_count", owned));

      // Set locally_owned_count property on all nodesets...
      Ioss::NodeSetContainer nss = region.get_nodesets();
      for (size_t i=0; i < nss.size(); i++) {
	Ioss::NodeSet *ns = nss[i];
	std::vector<INT> ids;
	t1 = timer();
	ns->get_field_data("ids_raw", ids);
	time_read += timer()-t1;
	data_read += ids.size();
	owned = 0;
	for (size_t n=0; n < ids.size(); n++) {
	  INT id = ids[n];
	  if (my_data[id-1] == my_processor)
	    owned++;
	}
	ns->property_add(Ioss::Property("locally_owned_count", owned));
      }
    }
  }
}
