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

#include <Ioss_CodeTypes.h>

#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <cstring>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifndef NO_XDMF_SUPPORT
#include <xdmf/Ioxf_Initializer.h>
#endif

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    enum FinalPressureType {FINAL, INITIAL, ZERO, OFFSET};
    bool debug;
    bool do_normals;
    bool reverse_normals;
    bool add_sset;
    bool convert_gage;
    FinalPressureType final_pressure;
    double delta_time;
    double maximum_time;
    double minimum_time;
    double offset_time;
    double offset_pressure;
    std::string working_directory;
  };

  void show_usage(const std::string &prog, bool add_sset);
  void show_step(int istep, double time);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void add_sideset(Ioss::Region &ss_region, Ioss::Region &region, Ioss::Region &output_region,
		   Globals &globals);
  void add_sideset_mesh_fields(Ioss::Region &ss_region, Ioss::Region &region, Ioss::Region &output_region,
			       Globals &globals);
  void add_sideset_transient_fields(Ioss::Region &ss_region, Ioss::Region &region, Ioss::Region &output_region,
				    Globals &globals);

  void transfer_sideset_field_data(Ioss::Region &ss_region, Ioss::Region &region,
				   Ioss::Region &output_region, Globals &globals);

  void transfer_fields(Ioss::GroupingEntity *ige,
		       Ioss::GroupingEntity *oge,
		       Ioss::Field::RoleType role,
		       const std::string &prefix = "");

  void transfer_field_data(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge,
			   Ioss::Field::RoleType role,
			   const std::string &prefix = "",
			   bool transfer_connectivity = true);

  void transfer_properties(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge);

  void output_normals(Ioss::Region &region, Ioss::Region &output_region,
		      bool reverse_normals);
  void calculate_normals(std::vector<double> &node_normal,
			 std::vector<double> &face_normal,
			 int num_elem, int num_node_per_elem,
			 std::vector<double> &coord,
			 std::vector<int>    &conn,
			 bool reverse_normals);
  void transfer_field_data_internal(Ioss::GroupingEntity *ige,
				    Ioss::GroupingEntity *oge,
				    const std::string &field_name);

  void file_copy(const std::string& inpfile, const std::string& input_type,
		 const std::string& outfile, const std::string& output_type,
		 const std::string& ss_file, Globals& globals);
}
// ========================================================================

namespace {
  std::string codename;
  std::string version = "$Revision$";
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  
  std::string in_type = "exodusII";
  std::string out_type = "exodusII";
  std::string ss_type = "exodusII";

  Globals globals;

  globals.do_normals = false;
  globals.reverse_normals = false;
  globals.add_sset   = true; // This is the CTH Pressure map executable
  globals.convert_gage = false;
  globals.final_pressure = Globals::FINAL;
  globals.delta_time      = 0.0;
  globals.maximum_time    = 0.0;
  globals.minimum_time    = 0.0;
  globals.offset_time     = 0.0;
  globals.offset_pressure = 0.0;
  
  codename = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos)
    codename = codename.substr(ind+1, codename.size());

  Ioss::Init::Initializer io;
#ifndef NO_XDMF_SUPPORT
  Ioxf::Initializer ioxf;
#endif

  globals.debug = false;

  std::string input_file;
  // Skip past any options...
  int i=1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("-directory", argv[i]) == 0 ||
	std::strcmp("--directory", argv[i]) == 0 ||	
	std::strcmp("-d", argv[i]) == 0) {
      i++;
      globals.working_directory = argv[i++];
    }
    else if (std::strcmp("-i", argv[i]) == 0) {
      i++;
      input_file = argv[i++];
    }
    else if (std::strcmp("-reverse", argv[i]) == 0) {
      globals.reverse_normals = true;
      i++;
    }
    else if (std::strcmp("--Reverse_Normals", argv[i]) == 0) {
      globals.reverse_normals = true;
      i++;
    }
    else if (std::strcmp("--Add_Surface_Fields", argv[i]) == 0) {
      globals.add_sset = true;
      i++;
    }
    else if (std::strcmp("--Calculate_Normals", argv[i]) == 0) {
      globals.do_normals = true;
      globals.add_sset = false;
      i++;
    }
    else if (std::strcmp("--Convert_Gage", argv[i]) == 0) {
      globals.convert_gage = true;
      i++;
    }
    else if (std::strcmp("--Final_Pressure", argv[i]) == 0) {
      i++;
      if (std::strcmp("ZERO", argv[i]) == 0)
	globals.final_pressure = Globals::ZERO;
      else if (std::strcmp("INITIAL", argv[i]) == 0)
	globals.final_pressure = Globals::INITIAL;
      else if (std::strcmp("FINAL", argv[i]) == 0)
	globals.final_pressure = Globals::FINAL;
      else if (std::strcmp("OFFSET", argv[i]) == 0)
	globals.final_pressure = Globals::OFFSET;
      i++;
    }
    else if (std::strcmp("--Final_Time_Delta", argv[i]) == 0) {
      i++;
      globals.delta_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("--Maximum_Time", argv[i]) == 0) {
      i++;
      globals.maximum_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("--Minimum_Time", argv[i]) == 0) {
      i++;
      globals.minimum_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("--Offset_Time", argv[i]) == 0) {
      i++;
      globals.offset_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("--Offset_Pressure", argv[i]) == 0) {
      i++;
      globals.offset_pressure = std::strtod(argv[i++], NULL);
    }

    // Found an option.  See if it has an argument...
    else if (i+1 < argc && argv[i+1][0] == '-') {
      // No argument, another option
      i++;
    } else {
      // Skip the argument...
      i += 2;
    }
  }

  std::string in_file;
  std::string sset_file;
  std::string out_file;

  // If a single argument is specified at end of command line, it is
  // input file which contains the names and types of the files to be
  // converted.
  // If two or more arguments are specified at end of command line, they are
  // the input and output files to be converted.
  // The file types are assumed to be as 'hardwired' above...

  std::string cwd = globals.working_directory;

  if (input_file.empty()) {
    if (argc - i == 3) {
      in_file   = Ioss::Utils::local_filename(argv[i++], in_type,  cwd);
      sset_file = Ioss::Utils::local_filename(argv[i++], ss_type,  cwd);
      out_file  = Ioss::Utils::local_filename(argv[i++], out_type, cwd);
    }
    else if (argc - i == 2) {
      in_file   = Ioss::Utils::local_filename(argv[i++], in_type,  cwd);
      out_file  = Ioss::Utils::local_filename(argv[i++], out_type, cwd);
    }
    else if (argc - i == 1) {
      input_file = Ioss::Utils::local_filename(argv[i++], "text", cwd);
    }
  }

  if (!input_file.empty()) {
    std::ifstream input(input_file.c_str());
    if (!input) {
      std::cerr << "Error opening file '" << input_file << "'.\n";
      show_usage(codename, globals.add_sset);
      return (EXIT_FAILURE);
    }
    // First line should be input_file_name input_file_type
    // If (add_sset), second line is sset file.
    // Last line should be output_file_name output_file_type
    std::string tmp;
    input >> tmp >> in_type;
    in_file = Ioss::Utils::local_filename(tmp, in_type, cwd);

    if (globals.add_sset) {
      input >> tmp >> ss_type;
      sset_file = Ioss::Utils::local_filename(tmp, ss_type, cwd);
    }

    input >> tmp >> out_type;
    out_file = Ioss::Utils::local_filename(tmp, out_type, cwd);
  }
  if (in_file.empty() || out_file.empty() || sset_file.empty()) {
    show_usage(codename, globals.add_sset);
    return(EXIT_FAILURE);
  }

  std::cerr << "Input:    '" << in_file  << "', Type: " << in_type  << '\n';
  if (globals.add_sset) {
    std::cerr << "Pressure: '" << sset_file  << "', Type: " << in_type  << '\n';
  }
  std::cerr << "Output:   '" << out_file << "', Type: " << out_type << '\n';
  std::cerr << '\n';

  if (globals.add_sset || globals.do_normals)
    globals.debug = false;

  if (!(globals.add_sset || globals.do_normals)) {
    std::cerr << "\n" << codename
	      << "ERROR: Either add sideset or do normals must be selected.\n";
    return EXIT_FAILURE;
  }

  file_copy(in_file, in_type, out_file, out_type, sset_file, globals);

  std::cerr << "\n" << codename << " execution successful.\n";
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog, bool add_sset)
  {
    if (add_sset)
      std::cerr << "\nUSAGE: " << prog << " in_file sset_file out_file\n";
    else
      std::cerr << "\nUSAGE: " << prog << " in_file out_file\n";

    std::cerr << "...or: " << prog << " command_file\n";
    std::cerr << "       version: " << version << "\n";
    Ioss::NameList db_types;
    Ioss::IOFactory::describe(&db_types);
    std::cerr << "\nSupports database types:\n\t";
    for (Ioss::NameList::const_iterator IF = db_types.begin(); IF != db_types.end(); ++IF) {
      std::cerr << *IF << "  ";
    }
    std::cerr << "\n\n";
  }

  void file_copy(const std::string& inpfile, const std::string& input_type,
		 const std::string& outfile, const std::string& output_type,
		 const std::string& ss_file, Globals& globals)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
						    (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == NULL || !dbi->ok()) {
      std::cerr << "ERROR: Could not open database '" << inpfile
			      << "' of type '" << input_type << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    //========================================================================
    // OUTPUT ...
    //========================================================================
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(output_type, outfile, Ioss::WRITE_RESTART,
						    (MPI_Comm)MPI_COMM_WORLD);
    if (dbo == NULL || !dbo->ok()) {
      std::cerr << "ERROR: Could not create output database '" << outfile
			      << "' of type '" << output_type << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");

    //========================================================================
    // Optional sideset file...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbs = NULL;
    Ioss::Region     *ss_region = NULL;
    if (globals.add_sset) {
      // Open the file containing the data which will be the new sideset...
      dbs = Ioss::IOFactory::create(input_type, ss_file, Ioss::READ_RESTART,
				    (MPI_Comm)MPI_COMM_WORLD);
      if (dbs == NULL || !dbs->ok()) {
	std::cerr << "ERROR: Could not open database '" << ss_file
				<< "' of type '" << input_type << "'\n";
	std::exit(EXIT_FAILURE);
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      ss_region = new Ioss::Region(dbs, "sideset_region");
    }


    // Set the qa information...
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    if (globals.debug) std::cerr << "DEFINING MODEL ... \n";
    if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
      std::cerr << "ERROR: Could not put output region into define model state\n";
      std::exit(EXIT_FAILURE);
    }

    // Get all properties of input database...
    transfer_properties(&region, &output_region);
    transfer_nodeblock(region, output_region, globals.debug);
    transfer_elementblock(region, output_region, globals.debug);
    transfer_nodesets(region, output_region, globals.debug);
    if (!globals.do_normals) {
      transfer_sidesets(region, output_region, globals.debug);
    }
    transfer_commsets(region, output_region, globals.debug);

    // If adding the sideset record from an external file, do it here...
    if (globals.add_sset) {
      add_sideset(*ss_region, region, output_region, globals);
    }

    if (globals.debug) std::cerr << "END STATE_DEFINE_MODEL... " << '\n';
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    if (globals.debug) std::cerr << "TRANSFERRING MESH FIELD DATA ... " << '\n';
    // Model defined, now fill in the model data...
    output_region.begin_mode(Ioss::STATE_MODEL);

    // Transfer MESH field_data from input to output...
    {
      Ioss::NodeBlock *nb = region.get_node_blocks()[0];
      assert(nb != NULL);

      if (nb != NULL) {
	std::string name = nb->name();
	if (globals.debug) std::cerr << name << ", ";

	// Find the corresponding output element_block...
	Ioss::NodeBlock *onb = output_region.get_node_block(name);
	assert(onb != NULL);
	transfer_field_data(nb, onb, Ioss::Field::MESH);
	transfer_field_data(nb, onb, Ioss::Field::ATTRIBUTE);
      }
      if (globals.debug) std::cerr << '\n';
    }

    // Now do the same for element blocks...
    {
      Ioss::ElementBlockContainer ebs = region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator I = ebs.begin();

      while (I != ebs.end()) {
	std::string name = (*I)->name();
	if (globals.debug) std::cerr << name << ", ";

	// Find the corresponding output element_block...
	Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	assert(oeb != NULL);

	if (oeb != NULL) {
	  transfer_field_data(*I, oeb, Ioss::Field::MESH);
	  transfer_field_data(*I, oeb, Ioss::Field::ATTRIBUTE);
	}
	++I;
      }
      if (globals.debug) std::cerr << '\n';
    }

    // Node Sets
    {
      Ioss::NodeSetContainer nss = region.get_nodesets();
      Ioss::NodeSetContainer::const_iterator I = nss.begin();
      while (I != nss.end()) {
	std::string name     = (*I)->name();
	if (globals.debug) std::cerr << name << ", ";

	// Find matching output nodeset
	Ioss::NodeSet *ons = output_region.get_nodeset(name);
	if (ons != NULL) {
	  transfer_field_data(*I, ons, Ioss::Field::MESH);
	  transfer_field_data(*I, ons, Ioss::Field::ATTRIBUTE);
	}
	++I;
      }
      if (globals.debug) std::cerr << '\n';
    }

    // Comm Sets
    {
      Ioss::CommSetContainer css = region.get_commsets();
      Ioss::CommSetContainer::const_iterator I = css.begin();
      while (I != css.end()) {
	std::string name  = (*I)->name();
	if (globals.debug) std::cerr << name << ", ";

	// Find matching output commset
	Ioss::CommSet *ocs = output_region.get_commset(name);
	if (ocs != NULL) {
	  transfer_field_data(*I, ocs, Ioss::Field::MESH);
	  transfer_field_data(*I, ocs, Ioss::Field::ATTRIBUTE);
	  transfer_field_data(*I, ocs, Ioss::Field::COMMUNICATION);
	}
	++I;
      }
      if (globals.debug) std::cerr << '\n';
    }

    // Side Sets
    if (!globals.do_normals) {
      Ioss::SideSetContainer fss = region.get_sidesets();
      Ioss::SideSetContainer::const_iterator I = fss.begin();
      while (I != fss.end()) {
	std::string name     = (*I)->name();
	if (globals.debug) std::cerr << name << ", ";

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
	    if (globals.debug) std::cerr << fbname << ", ";
	    Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

	    if (ofb != NULL) {
	      transfer_field_data(*J, ofb, Ioss::Field::MESH, "", false);
	      transfer_field_data(*J, ofb, Ioss::Field::ATTRIBUTE, "", false);
	    }
	    ++J;
	  }
	}
	++I;
      }
      if (globals.debug) std::cerr << '\n';
    }

    if (globals.add_sset) {
      add_sideset_mesh_fields(*ss_region, region, output_region, globals);
    }

    if (globals.debug) std::cerr << "END STATE_MODEL... " << '\n';
    output_region.end_mode(Ioss::STATE_MODEL);

    if (globals.do_normals) {
      output_normals(region, output_region, globals.reverse_normals);
    } else if (globals.add_sset) {
      if (ss_region->property_exists("state_count") &&
	  ss_region->get_property("state_count").get_int() > 0) {
	output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
	add_sideset_transient_fields(*ss_region, region, output_region, globals);
	output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
      }
    } else {
      std::cerr << "Internal Error\n";      
      std::abort();
    }

    if (globals.debug) std::cerr << "TRANSFERRING TRANSIENT FIELDS ... " << '\n';
    output_region.begin_mode(Ioss::STATE_TRANSIENT);
    // Get the timesteps from the input database.  Step through them
    // and transfer fields to output database...

    if (globals.do_normals) {
      // Do nothing, normals were already output above...

    } else if (globals.add_sset) {
      transfer_sideset_field_data(*ss_region, region, output_region, globals);
    } else {
      std::cerr << "Internal Error\n";      
      std::abort();
    }
    if (globals.debug) std::cerr << "END STATE_TRANSIENT... " << '\n';
    output_region.end_mode(Ioss::STATE_TRANSIENT);
    delete ss_region;
  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeBlockContainer    nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i = nbs.begin();
    int id = 1;
    while (i != nbs.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr << name << ", ";
      int    num_nodes = (*i)->get_property("entity_count").get_int();
      int    degree    = (*i)->get_property("component_degree").get_int();
      if (!debug) {
	std::cerr << " Number of coordinates per node       =" << std::setw(9) << degree << "\n";
	std::cerr << " Number of nodes                      =" << std::setw(9) << num_nodes << "\n";
      }

      Ioss::NodeBlock *nb = new Ioss::NodeBlock(output_region.get_database(), name, num_nodes, degree);
      output_region.add(nb);

      transfer_properties(*i, nb);
      transfer_fields(*i, nb, Ioss::Field::MESH);
      transfer_fields(*i, nb, Ioss::Field::ATTRIBUTE);
      ++i;
      ++id;
    }
    if (debug) std::cerr << '\n';
  }

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    int total_elements = 0;
    while (i != ebs.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr << name << ", ";
      std::string type      = (*i)->get_property("topology_type").get_string();
      int    num_elem  = (*i)->get_property("entity_count").get_int();
      int    num_attrib= (*i)->get_property("attribute_count").get_int();
      total_elements += num_elem;

      Ioss::ElementBlock *eb = new Ioss::ElementBlock(output_region.get_database(), name, type,
						      num_elem, num_attrib);
      output_region.add(eb);
      transfer_properties(*i, eb);
      transfer_fields(*i, eb, Ioss::Field::MESH);
      transfer_fields(*i, eb, Ioss::Field::ATTRIBUTE);
      ++i;
    }
    if (!debug) {
      std::cerr << " Number of elements                   =" << std::setw(9) << total_elements << "\n";
      std::cerr << " Number of element blocks             =" << std::setw(9) << ebs.size() << "\n\n";
    } else {
      std::cerr << '\n';
    }
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::SideSetContainer      fss = region.get_sidesets();
    Ioss::SideSetContainer::const_iterator i = fss.begin();
    int total_sides = 0;
    while (i != fss.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr << name << ", ";
      Ioss::SideSet *surf = new Ioss::SideSet(output_region.get_database(), name);

      Ioss::SideBlockContainer fbs = (*i)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator j = fbs.begin();
      while (j != fbs.end()) {
	std::string fbname    = (*j)->name();
	if (debug) std::cerr << fbname << ", ";
	std::string fbtype    = (*j)->get_property("topology_type").get_string();
	std::string partype   = (*j)->get_property("parent_topology_type").get_string();
	int    num_side  = (*j)->get_property("entity_count").get_int();
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
      std::cerr << " Number of element side sets          =" << std::setw(9) << fss.size() << "\n";
      std::cerr << "     Number of element sides          =" << std::setw(9) << total_sides << "\n";
    } else {
      std::cerr << '\n';
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeSetContainer      nss = region.get_nodesets();
    Ioss::NodeSetContainer::const_iterator i = nss.begin();
    int total_nodes = 0;
    while (i != nss.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr << name << ", ";
      int    count     = (*i)->get_property("entity_count").get_int();
      total_nodes += count;
      Ioss::NodeSet *ns = new Ioss::NodeSet(output_region.get_database(), name, count);
      output_region.add(ns);
      transfer_properties(*i, ns);
      transfer_fields(*i, ns, Ioss::Field::MESH);
      transfer_fields(*i, ns, Ioss::Field::ATTRIBUTE);
      ++i;
    }
    if (!debug) {
      std::cerr << " Number of nodal point sets           =" << std::setw(9) << nss.size() << "\n";
      std::cerr << "     Length of node list              =" << std::setw(9) << total_nodes << "\n";
    } else {
      std::cerr << '\n';
    }
  }

  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::CommSetContainer      css = region.get_commsets();
    Ioss::CommSetContainer::const_iterator i = css.begin();
    while (i != css.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr << name << ", ";
      std::string type      = (*i)->get_property("entity_type").get_string();
      int    count     = (*i)->get_property("entity_count").get_int();
      Ioss::CommSet *cs = new Ioss::CommSet(output_region.get_database(), name, type, count);
      output_region.add(cs);
      transfer_properties(*i, cs);
      transfer_fields(*i, cs, Ioss::Field::MESH);
      transfer_fields(*i, cs, Ioss::Field::ATTRIBUTE);
      transfer_fields(*i, cs, Ioss::Field::COMMUNICATION);
      ++i;
    }
    if (debug) std::cerr << '\n';
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
      if (field_name != "ids" && !oge->field_exists(field_name) &&
	  (prefix.length() == 0 || std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
	// If the field does not already exist, add it to the output node block
	Ioss::Field field = ige->get_field(field_name);
	oge->field_add(field);
      }
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge,
			   Ioss::Field::RoleType role,
			   const std::string &prefix,
			   bool transfer_connectivity)
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
      if (!transfer_connectivity && field_name == "connectivity")
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
    assert (isize == oge->get_field(field_name).get_size());

    data.resize(isize);

    ige->get_field_data(field_name, &data[0], isize);
    oge->put_field_data(field_name, &data[0], isize);
  }

  void transfer_properties(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge)
  {
    Ioss::NameList names;
    ige->property_describe(&names);

    // Iterate through properties and transfer to output database...
    Ioss::NameList::const_iterator I;
    for (I = names.begin(); I != names.end(); ++I) {
      if (!oge->property_exists(*I))
	oge->property_add(ige->get_property(*I));
    }
  }

  void add_sideset(Ioss::Region &ss_region, Ioss::Region &/* region */, Ioss::Region &output_region, Globals &globals)
  {
    Ioss::SideSet *pressures = new Ioss::SideSet(output_region.get_database(), "cth_pressures");
    // Each element block in the sset file will be a surside in the mesh file...
    {
      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {
	std::string name      = (*i)->name();
	name = "ss" + name;

	if (globals.debug) std::cerr << name << ", ";
	int    num_elem  = (*i)->get_property("entity_count").get_int();
	std::string type      = (*i)->get_property("topology_type").get_string();

	// Should be able to get this from the input mesh element blocks...
	std::string partype   = "unknown";
	Ioss::SideBlock *fb = new Ioss::SideBlock(output_region.get_database(), name, type,
						  partype, num_elem);
	pressures->add(fb);
	++i;
      }
      output_region.add(pressures);
    }
  }

  void add_sideset_mesh_fields(Ioss::Region &ss_region, Ioss::Region &/* region */,
			       Ioss::Region &output_region, Globals &globals)
  {
    // Each element block in the sset file will be a surface in the mesh file...
    {
      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {
	std::string name      = (*i)->name();
	name = "ss"+name;
	Ioss::SideBlock *fb = output_region.get_sideblock(name);
	if (fb == NULL) {
	  std::cerr << "INTERNAL ERROR: Could not find sideblock named '" << name << "'\n";
	  std::exit(EXIT_FAILURE);
	}

	if (globals.debug) std::cerr << name << ", ";

	// The "skin" field should exist on each element block.  It
	// contains the map back to the element/local face in the
	// original (region/output_region) mesh.  It is in the exact
	// same format as a sidesets "element_side" field...
	//
	// See if the field exists...
	if ((*i)->field_exists("skin")) {
	  Ioss::Field skin = (*i)->get_field("skin");

	  int isize = skin.get_size();
	  int osize = fb->get_field("element_side").get_size();
	  assert (isize == osize);

	  data.resize(isize);

	  (*i)->get_field_data("skin", &data[0], isize);
	  fb->put_field_data("element_side", &data[0], osize);
	  ++i;
	} else {
	  std::cerr << "ERROR: Field 'skin' does not exist on element block '"
				  << name << "'.\n";
	  std::exit(EXIT_FAILURE);
	}
      }
    }
  }

  void transfer_sideset_field_data(Ioss::Region &ss_region, Ioss::Region &/* region */,
				   Ioss::Region &output_region, Globals &globals)
  {
    int step_count = ss_region.get_property("state_count").get_int();
    int istep = 1; // Need this to exist after end of step_count loop.

    int initial_state = 1;
    // Find first state >= globals.minimum_time;
    for ( ; initial_state <= step_count; initial_state++) {
      double time = ss_region.get_state_time(initial_state);
      if (time >= globals.minimum_time) {
	globals.minimum_time = time;
	break;
      }
    }

    const std::string cth_pressure = "cth_pressure";

    // If the user specified the "Convert_Gage" option, then we need
    // to save the field values at time 'globals.minimum_time' so they
    // can be subtracted from all field data.
    std::map<std::string, double*> time_zero_field_data;
    if (globals.convert_gage) {
      ss_region.begin_state(1);

      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {

	// The gage pressure conversion is currently only applied to the field "cth_pressure"
	std::string name = (*i)->name();
	if ((*i)->field_exists(cth_pressure)) {
	  int isize = (*i)->get_field(cth_pressure).get_size();
	  void *zdata = new char[isize];
	  time_zero_field_data[name] = (double*)zdata;
	  (*i)->get_field_data(cth_pressure, zdata, isize);
	} else {
	  time_zero_field_data[name] = (double*)NULL;
	}
	++i;
      }
      ss_region.end_state(1);
    }

    // If globals.offset_time is non-zero, then the output times will be
    // offset by the specified amount.  We need to add a time=0.0 step
    // with all pressures set to globals.offset_pressure.
    if (globals.offset_time > 0.0) {
      int ostep = output_region.add_state(0.0);
      output_region.begin_state(ostep);
      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {
	std::string eb_name = (*i)->name();
	std::string name = "ss" + eb_name;
	Ioss::SideBlock *fb = output_region.get_sideblock(name);
	if (fb == NULL) {
	  std::cerr << "INTERNAL ERROR: Could not find sideblock named '" << name << "'\n";
	  std::exit(EXIT_FAILURE);
	}
	
	Ioss::NameList state_fields;
	Ioss::NameList::const_iterator IF;
	(*i)->field_describe(Ioss::Field::TRANSIENT, &state_fields);
	
	for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
	  std::string field_name = *IF;
	  // NOTE: Only dealing with the "cth_" fields here.
	  // If there are other fields, we probably have an invalid
	  // output database...
	  if (std::strncmp("cth_", field_name.c_str(), 4) == 0) {
	    int isize = (*i)->get_field(field_name).get_size();
	    int count = (*i)->get_field(field_name).raw_count();
	    data.resize(isize);
	    double *rdata = (double*)&data[0];
	    for (int ii=0; ii < count; ii++) {
	      rdata[ii] = globals.offset_pressure;
	    }
	    fb->put_field_data(field_name, &data[0], isize);
	    break;
	  }
	}
	++i;
      }
      output_region.end_state(ostep);
    }
      
    for (istep = initial_state; istep <= step_count; istep++) {
      double time = ss_region.get_state_time(istep);
      if (time < globals.minimum_time)
	continue;

      if (globals.maximum_time > 0.0 && time > globals.maximum_time)
	break;

      int ostep = output_region.add_state(time-globals.minimum_time+globals.offset_time);
      show_step(istep, time);

      output_region.begin_state(ostep);
      ss_region.begin_state(istep);
      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {
	std::string eb_name = (*i)->name();
	std::string name = "ss" + eb_name;
	Ioss::SideBlock *fb = output_region.get_sideblock(name);
	if (fb == NULL) {
	  std::cerr << "INTERNAL ERROR: Could not find side block named '" << name << "'\n";
	  std::exit(EXIT_FAILURE);
	}

	Ioss::NameList state_fields;
	Ioss::NameList::const_iterator IF;
	(*i)->field_describe(Ioss::Field::TRANSIENT, &state_fields);

	for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
	  std::string field_name = *IF;
	  if (globals.convert_gage && field_name == cth_pressure) {
	    // Subtract the time zero pressures (stored in
	    // time_zero_field_data) from each time step
	    double* zdata = time_zero_field_data[eb_name];
	    assert(zdata != NULL);

	    int isize = (*i)->get_field(field_name).get_size();
	    int count = (*i)->get_field(field_name).raw_count();
	    data.resize(isize);
	    double *rdata = (double*)&data[0];

	    (*i)->get_field_data(field_name, &data[0], isize);
	    for (int ii=0; ii < count; ii++) {
	      rdata[ii] -= zdata[ii];
	    }

	    if (globals.offset_pressure != 0.0) {
	      for (int ii=0; ii < count; ii++) {
		rdata[ii] += globals.offset_pressure;
	      }
	    }
	    fb->put_field_data(field_name, &data[0], isize);
	  }
	  else if (globals.offset_pressure != 0.0 && field_name == cth_pressure) {
	    int isize = (*i)->get_field(field_name).get_size();
	    int count = (*i)->get_field(field_name).raw_count();
	    data.resize(isize);
	    double *rdata = (double*)&data[0];

	    (*i)->get_field_data(field_name, &data[0], isize);
	    for (int ii=0; ii < count; ii++) {
	      rdata[ii] += globals.offset_pressure;
	    }
	    fb->put_field_data(field_name, &data[0], isize);
	  }
	  else if (std::strncmp("cth_", field_name.c_str(), 4) == 0) {
	    assert(fb->field_exists(field_name));
	    transfer_field_data_internal(*i, fb, field_name);
	  }
	}
	++i;
      }
      ss_region.end_state(istep);
      output_region.end_state(ostep);
    }

    step_count = istep-1;

    // See if special handling of the last step is specified...
    if (globals.final_pressure == Globals::INITIAL ||
	globals.final_pressure == Globals::ZERO ||
	globals.final_pressure == Globals::OFFSET) {
      double time = ss_region.get_state_time(step_count);
      double delta_time;
      if (globals.delta_time > 0.0)
	delta_time = globals.delta_time;
      else
	delta_time = time - ss_region.get_state_time(step_count-1);
      time += delta_time;
      int ostep = output_region.add_state(time-globals.minimum_time+globals.offset_time);

      output_region.begin_state(ostep);
      ss_region.begin_state(initial_state);

      // Repeat the data from the first step at the last step.  This is to bring the model back to an
      // "equilibrium" state in case the CTH analysis was not run out to an equilibrium
      // state. If ZERO was specified, then it simply zeros out the pressure field at the last step.
      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {
	std::string eb_name = (*i)->name();
	std::string name = "ss" + eb_name;
	Ioss::SideBlock *fb = output_region.get_sideblock(name);
	if (fb == NULL) {
          std::ostringstream msg ;
	  msg << " INTERNAL_ERROR: Could not find sideblock '" << name << "'\n";
	  throw std::runtime_error( msg.str() );
	}

	Ioss::NameList state_fields;
	Ioss::NameList::const_iterator IF;
	(*i)->field_describe(Ioss::Field::TRANSIENT, &state_fields);

	for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
	  std::string field_name = *IF;
	  if (std::strncmp("cth_", field_name.c_str(), 4) == 0) {
	    if (field_name == cth_pressure &&
		(globals.final_pressure == Globals::ZERO ||
		 globals.final_pressure == Globals::OFFSET ||
		 globals.convert_gage)) {

	      // If convert_gage is true and we are outputting the
	      // cth_pressure variable, then for both the INITIAL and
	      // ZERO cases, we output a zero-filled field...

	      double value = 0.0;
	      if (globals.final_pressure == Globals::OFFSET)
		value = globals.offset_pressure;
	      
	      int isize = (*i)->get_field(field_name).get_size();
	      int count = (*i)->get_field(field_name).raw_count();
	      data.resize(isize);
	      double *rdata = (double*)&data[0];
	      for (int ii=0; ii < count; ii++) {
		rdata[ii] = value;
	      }
	      fb->put_field_data(field_name, &data[0], isize);
	    }
	    else {
	      // Replicate first field as last field...
	      assert(fb->field_exists(field_name));
	      transfer_field_data_internal(*i, fb, field_name);
	    }
	  }
	}
	++i;
      }
      output_region.end_state(ostep);
      ss_region.end_state(initial_state);
    }
    if (globals.convert_gage) {
      std::map<std::string, double*>::iterator i = time_zero_field_data.begin();
      std::map<std::string, double*>::iterator ie = time_zero_field_data.end();
      while (i != ie) {
	delete [] (char *)(*i).second;
	++i;
      }
    }
  }

  void add_sideset_transient_fields(Ioss::Region &ss_region, Ioss::Region &/* region */,
				    Ioss::Region &output_region, Globals &globals)
  {
    // Each element block in the sset file will be a surface in the mesh file...
    {
      Ioss::ElementBlockContainer ebs = ss_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
      while (i != ebs.end()) {
	std::string name      = (*i)->name();
	name = "ss"+name;
	Ioss::SideBlock *fb = output_region.get_sideblock(name);
	assert(fb != NULL);

	if (globals.debug) std::cerr << name << ", ";

	// Each element variable in the sset file which begins with
	// "cth_" will be a sideset variable in the outptut file...
	transfer_fields((*i), fb, Ioss::Field::TRANSIENT, "cth_");
	++i;
      }
    }
  }

#ifndef MGCVECTOR3_H
#define MGCVECTOR3_H

  // Magic Software, Inc.
  // http://www.magic-software.com
  // Copyright (c) 2000, All Rights Reserved
  //
  // Source code from Magic Software is supplied under the terms of a license
  // agreement and may not be copied or disclosed except in accordance with the
  // terms of that agreement.  The various license agreements may be found at
  // the Magic Software web site.  This file is subject to the license
  //
  // FREE SOURCE CODE
  // http://www.magic-software.com/License/free.pdf

  class Vector3
  {
  public:
    // construction
    Vector3 ();
    Vector3 (double fX, double fY, double fZ);
    explicit Vector3 (double Coordinate[3]);
    Vector3 (const Vector3& from);

    double x, y, z;

    // assignment and comparison
    Vector3& operator= (const Vector3& from);
    bool operator== (const Vector3& from) const;
    bool operator!= (const Vector3& from) const;
    void set(double fX, double fY, double fZ);
    void set(double Coordinate[3]);
    Vector3& reverse();

    // arithmetic operations
    Vector3 operator- () const;

    // arithmetic updates
    Vector3& operator+= (const Vector3& from);
    Vector3& operator-= (const Vector3& from);
    Vector3& operator*= (double scalar);
    Vector3& operator/= (double scalar);

    // vector operations
    double length () const;
    double squared_length () const;
    double dot (const Vector3& from) const;
    double normalize (double tolerance = 1e-06);
    Vector3 cross (const Vector3& from) const;
    static Vector3 plane_normal(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);
  };

  const Vector3 operator* (double scalar, const Vector3& vec);
  const Vector3 operator+ (const Vector3& vec, const Vector3& vec2);
  const Vector3 operator- (const Vector3& vec, const Vector3& vec2);
  const Vector3 operator* (double scalar, const Vector3& vec);
  const Vector3 operator/ (const Vector3& vec, double scalar);

  //----------------------------------------------------------------------------
  inline Vector3 Vector3::cross (const Vector3& from) const
  {
    return Vector3(y*from.z - z*from.y,
		   z*from.x - x*from.z,
		   x*from.y - y*from.x);
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator+= (const Vector3& from)
  {
    x += from.x;
    y += from.y;
    z += from.z;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator-= (const Vector3& from)
  {
    x -= from.x;
    y -= from.y;
    z -= from.z;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator*= (double scalar)
  {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }
#endif

  void output_normals(Ioss::Region &region, Ioss::Region &output_region, bool reverse_normals)
  {
    // Define output fields...
    {
      const Ioss::VariableType *v3d = Ioss::VariableType::factory("vector_3d");
      output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
      Ioss::NodeBlock *nb = (*output_region.get_node_blocks().begin());
      int num_nodes  = nb->get_property("entity_count").get_int();
      Ioss::Field node_normal("node_normal", Ioss::Field::REAL,
			      v3d, Ioss::Field::TRANSIENT, num_nodes);
      nb->field_add(node_normal);

      // Iterate over the element blocks and calculate both node normals and face normals...
      Ioss::ElementBlockContainer ebs = output_region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator ib = ebs.begin();
      while (ib != ebs.end()) {
	Ioss::ElementBlock *eb = *ib; ++ib;
	int num_elem  = eb->get_property("entity_count").get_int();
	Ioss::Field face_normal("face_normal", Ioss::Field::REAL,
				v3d, Ioss::Field::TRANSIENT, num_elem);
	eb->field_add(face_normal);
      }
      output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }

    output_region.begin_mode(Ioss::STATE_TRANSIENT);
    int ostep = output_region.add_state(0.0);
    output_region.begin_state(ostep);

    Ioss::NodeBlock *nb  = (*region.get_node_blocks().begin());
    Ioss::NodeBlock *nbo = (*output_region.get_node_blocks().begin());

    // Get the nodal coordinates...
    int num_nodes  = nb->get_property("entity_count").get_int();
    int coord_size = nb->get_field("mesh_model_coordinates").get_size();

    std::vector<double> coord(3*num_nodes);
    assert(3*num_nodes*sizeof(double) == (size_t)coord_size);

    nb->get_field_data("mesh_model_coordinates", &coord[0], coord_size);

    // Also get an array for the average nodal normal vector...
    std::vector<double> node_normal(3*num_nodes);
    std::fill(node_normal.begin(), node_normal.end(), 0.0);

    // Iterate over the element blocks and calculate both node normals and face normals...
    std::vector<int>    conn;
    std::vector<double> face_normal;
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();

    Ioss::ElementBlockContainer::const_iterator ib  = ebs.begin();
    while (ib != ebs.end()) {
      Ioss::ElementBlock *eb = *ib; ++ib;
      std::string name = (*eb).name();

      Ioss::ElementBlock *ebo = output_region.get_element_block(name);
      if (ebo == NULL) {
	std::cerr << "INTERNAL ERROR: Could not find element block named '" << name << "'\n";
	std::exit(EXIT_FAILURE);
      }

      int num_elem  = eb->get_property("entity_count").get_int();
      int num_node_per_elem = eb->topology()->number_nodes();

      // Get the connectivity array...
      conn.resize(num_elem * num_node_per_elem);
      eb->get_field_data("connectivity", conn);

      // Connectivity is in global id space; change to local...
      for (int i=0; i < num_elem * num_node_per_elem; i++) {
	int local = region.node_global_to_local(conn[i]);
	conn[i] = local-1;
      }
      // Get an array for storing the face normals for this block...
      face_normal.resize(3*num_elem);
      std::fill(face_normal.begin(), face_normal.end(), 0.0);

      calculate_normals(node_normal, face_normal, num_elem, num_node_per_elem, coord, conn,
			reverse_normals);
      ebo->put_field_data("face_normal", face_normal);
    }

    int nsize = node_normal.size();
    for (int i=0; i < nsize; i+=3) {
      Vector3 a(node_normal[i+0], node_normal[i+1], node_normal[i+2]);
      a.normalize();

      node_normal[i+0] = a.x;
      node_normal[i+1] = a.y;
      node_normal[i+2] = a.z;
    }
    nbo->put_field_data("node_normal", node_normal);
    output_region.end_state(ostep);
    output_region.end_mode(Ioss::STATE_TRANSIENT);
  }

  void calculate_normals(std::vector<double> &node_normal,
			 std::vector<double> &face_normal,
			 int num_elem, int num_node_per_elem,
			 std::vector<double> &coord,
			 std::vector<int>    &conn,
			 bool reverse_normals)
  {
    // Iterate the connectivity array and calculate normals...
    // The elements should all be shells with the outward normal
    // defined by the node ordering 1.2.3.
    // The face normal will be the average of the face-corner node
    // normals for a 4-node face.  For a triangular face, all
    // normals are the same...

    for (int iel = 0; iel < num_elem; iel++) {
      int ioff = iel * num_node_per_elem;

      // Triangular faces...
      if (num_node_per_elem == 3) {
	Vector3 local[3];
	for (int i=0; i < 3; i++) {
	  int node = conn[ioff+i];
	  local[i].set(coord[node*3+0], coord[node*3+1], coord[node*3+2]);
	}

	Vector3 plnorm = Vector3::plane_normal(local[0], local[1], local[2]);
	plnorm.normalize();
	if (reverse_normals)
	  plnorm.reverse();

	face_normal[iel*3+0] = plnorm.x;
	face_normal[iel*3+1] = plnorm.y;
	face_normal[iel*3+2] = plnorm.z;

	for (int i=0; i < 3; i++) {
	  int node = conn[ioff+i];
	  node_normal[node*3+0] += plnorm.x;
	  node_normal[node*3+1] += plnorm.y;
	  node_normal[node*3+2] += plnorm.z;
	}
      } else {

	// Quadrilateral faces...
	assert(num_node_per_elem == 4);
	Vector3 local[4];
	for (int i=0; i < 4; i++) {
	  int node = conn[ioff+i];
	  local[i].set(coord[node*3+0], coord[node*3+1], coord[node*3+2]);
	}

	for (int i=0; i < 4; i++) {
	  // at node 0 -- vector from 3-0 X 0-1
	  // at node 1 -- vector from 0-1 X 1-2
	  // at node 2 -- vector from 1-2 X 2-3
	  // at node 3 -- vector from 2-3 X 3-0
	  int nb = (i+3)%4;
	  int na = (i+1)%4;

	  Vector3 a = Vector3::plane_normal(local[nb], local[i], local[na]);
	  a.normalize();
	  if (reverse_normals)
	    a.reverse();

	  face_normal[iel*3+0] += a.x;
	  face_normal[iel*3+1] += a.y;
	  face_normal[iel*3+2] += a.z;

	  int node = conn[ioff+i];
	  node_normal[node*3+0] += a.x;
	  node_normal[node*3+1] += a.y;
	  node_normal[node*3+2] += a.z;
	}
      }
      Vector3 a(face_normal[iel*3+0], face_normal[iel*3+1], face_normal[iel*3+2]);
      a.normalize();
      face_normal[iel*3+0] = a.x;
      face_normal[iel*3+1] = a.y;
      face_normal[iel*3+2] = a.z;
    }
  }

  void show_step(int istep, double time)
  {
    std::cerr.setf(std::ios::scientific);
    std::cerr.setf(std::ios::showpoint);
    std::cerr << "     Time step " << std::setw(5) << istep
			    << " at time " << std::setprecision(5) << time << '\n';
  }
}

// Magic Software, Inc.
// http://www.magic-software.com
// Copyright (c) 2000, All Rights Reserved
//
// Source code from Magic Software is supplied under the terms of a license
// agreement and may not be copied or disclosed except in accordance with the
// terms of that agreement.  The various license agreements may be found at
// the Magic Software web site.  This file is subject to the license
//
// FREE SOURCE CODE
// http://www.magic-software.com/License/free.pdf


//----------------------------------------------------------------------------
Vector3::Vector3 ()
  : x(0.0), y(0.0), z(0.0)
{
}

//----------------------------------------------------------------------------
Vector3::Vector3 (double fX, double fY, double fZ)
  : x(fX), y(fY), z(fZ)
{}

//----------------------------------------------------------------------------
Vector3::Vector3 (double Coordinate[3])
  : x(Coordinate[0]), y(Coordinate[1]), z(Coordinate[2])
{}

//----------------------------------------------------------------------------
Vector3::Vector3 (const Vector3& from)
  : x(from.x), y(from.y), z(from.z)
{}

void Vector3::set(double fX, double fY, double fZ)
{
    x = fX;
    y = fY;
    z = fZ;
}

void Vector3::set(double Coordinate[3])
{
    x = Coordinate[0];
    y = Coordinate[1];
    z = Coordinate[2];
}

Vector3& Vector3::operator= (const Vector3& from)
{
    x = from.x;
    y = from.y;
    z = from.z;
    return *this;
}

Vector3& Vector3::reverse()
{
  x = -x;
  y = -y;
  z = -z;
  return *this;
}


bool Vector3::operator== (const Vector3& from) const
{
    return ( x == from.x && y == from.y && z == from.z );
}

bool Vector3::operator!= (const Vector3& from) const
{
    return ( x != from.x || y != from.y || z != from.z );
}

const Vector3 operator+ (const Vector3& lhs, const Vector3& rhs)
{
  Vector3 Sum(lhs);
  return Sum += rhs;
}

const Vector3 operator- (const Vector3& lhs, const Vector3& rhs)
{
  Vector3 Diff(lhs);
  return Diff -= rhs;
}

const Vector3 operator* (const Vector3& lhs, double scalar)
{
  Vector3 Prod(lhs);
  return Prod *= scalar;
}

Vector3 Vector3::operator- () const
{
    Vector3 Neg;
    return Neg *= -1.0;
}

const Vector3 operator* (double scalar, const Vector3& from)
{
  Vector3 Prod(from);
  return Prod *= scalar;
}

double Vector3::squared_length () const
{
    return x*x + y*y + z*z;
}

double Vector3::dot (const Vector3& from) const
{
    return x*from.x + y*from.y + z*from.z;
}

Vector3 operator/ (const Vector3& lhs, double scalar)
{
  Vector3 Quot(lhs);
  return Quot /= scalar;
}

Vector3& Vector3::operator/= (double scalar)
{
  if ( scalar != 0.0 ) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
  } else {
    x = HUGE_VAL;
    y = HUGE_VAL;
    z = HUGE_VAL;
  }

  return *this;
}

double Vector3::length () const
{
  return sqrt(x*x + y*y + z*z);
}

double Vector3::normalize (double tolerance)
{
  double mylength = length();

  if ( mylength > tolerance ) {
    x /= mylength;
    y /= mylength;
    z /= mylength;
  } else {
    mylength = 0.0;
  }

  return mylength;
}

Vector3 Vector3::plane_normal(const Vector3 &v1,
			      const Vector3 &v2,
			      const Vector3 &v3)
{
  Vector3 v32 = v3;
  v32 -= v2;
  Vector3 v12 = v1;
  v12 -= v2;
  return v32.cross(v12);
}
