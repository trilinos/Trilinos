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
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Transform.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifndef NO_XDMF_SUPPORT
#include <xdmf/Ioxf_Initializer.h>
#endif

#define OUTPUT std::cerr

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    bool debug;
    bool do_transform_fields;
    double maximum_time;
    double minimum_time;
    int  surface_split_type;
    std::string working_directory;
  };

  void show_usage(const std::string &prog);
  void show_step(int istep, double time);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug);

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

  void transform_fields(Ioss::GroupingEntity *ige,
			Ioss::GroupingEntity *oge,
			Ioss::Field::RoleType role);

  void transform_field_data(Ioss::GroupingEntity *ige,
			    Ioss::GroupingEntity *oge,
			    Ioss::Field::RoleType role);
  void transfer_field_data_internal(Ioss::GroupingEntity *ige,
				    Ioss::GroupingEntity *oge,
				    const std::string &field_name);

  void file_copy(const std::string& inpfile, const std::string& input_type,
		 const std::string& outfile, const std::string& output_type,
		 Globals& globals);
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
#endif
  
  std::string in_type = "exodusII";
  std::string out_type = "exodusII";

  Globals globals;

  globals.do_transform_fields = false;
  globals.maximum_time = 0.0;
  globals.minimum_time = 0.0;
  globals.surface_split_type = 1;
  
  codename = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos)
    codename = codename.substr(ind+1, codename.size());

  // Check the program name to see if of the form 'exosaf' or 'safexo'
  // and if it is, set the in_type and out_type accordingly...
  if (std::strncmp(codename.c_str(), "exoxdmf", 7) == 0) {
    codename = "exoxdmf";
    in_type  = "exodusII";
    out_type = "xdmf";
  }
  else {
    codename = "io_shell";
  }

  Ioss::Init::Initializer io;
#ifndef NO_XDMF_SUPPORT
  Ioxf::Initializer ioxf;
#endif

  std::string input_file;
  globals.debug = false;

  // Skip past any options...
  int i=1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("-directory", argv[i]) == 0 ||
	std::strcmp("-d", argv[i]) == 0) {
      i++;
      globals.working_directory = argv[i++];
    }
    else if (std::strcmp("--in_type", argv[i]) == 0) {
      i++;
      in_type = argv[i++];
    }
    else if (std::strcmp("--out_type", argv[i]) == 0) {
      i++;
      out_type = argv[i++];
    }
    else if (std::strcmp("-i", argv[i]) == 0) {
      i++;
      input_file = argv[i++];
    }
    else if (std::strcmp("--debug", argv[i]) == 0) {
      i++;
      globals.debug = true;
    }
    else if (std::strcmp("--Maximum_Time", argv[i]) == 0) {
      i++;
      globals.maximum_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("--Minimum_Time", argv[i]) == 0) {
      i++;
      globals.minimum_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("--Surface_Split_Scheme", argv[i]) == 0) {
      i++;
      char *split_scheme = argv[i++];
      if (std::strcmp(split_scheme, "TOPOLOGY") == 0)
	globals.surface_split_type = 1;
      else if (std::strcmp(split_scheme, "ELEMENT_BLOCK") == 0)
	globals.surface_split_type = 2;
      else if (std::strcmp(split_scheme, "NO_SPLIT") == 0)
	globals.surface_split_type = 3;
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
  std::string out_file;

  // If a single argument is specified at end of command line, it is
  // input file which contains the names and types of the files to be
  // converted.
  // If two or more arguments are specified at end of command line, they are
  // the input and output files to be converted.
  // The file types are assumed to be as 'hardwired' above...

  std::string cwd = globals.working_directory;

  if (input_file.empty()) {
    if (argc - i == 2) {
      in_file   = Ioss::Utils::local_filename(argv[i++], in_type, cwd);
      out_file  = Ioss::Utils::local_filename(argv[i++], out_type, cwd);
    }
    else if (argc - i == 1) {
      input_file = Ioss::Utils::local_filename(argv[i++], "text", cwd);
    }
  }

  if (!input_file.empty()) {
    std::ifstream input(input_file.c_str());
    if (!input) {
      OUTPUT << "Error opening file '" << input_file << "'.\n";
      show_usage(codename);
      return (EXIT_FAILURE);
    }

    // First line should be input_file_name input_file_type
    // Last line should be output_file_name output_file_type
    std::string tmp;
    input >> tmp >> in_type;
    in_file = Ioss::Utils::local_filename(tmp, in_type, cwd);

    input >> tmp >> out_type;
    out_file = Ioss::Utils::local_filename(tmp, out_type, cwd);
  }

  if (in_file.empty() || out_file.empty()) {
    show_usage(codename);
    return(EXIT_FAILURE);
  }

  OUTPUT << "Input:    '" << in_file  << "', Type: " << in_type  << '\n';
  OUTPUT << "Output:   '" << out_file << "', Type: " << out_type << '\n';
  OUTPUT << '\n';

  file_copy(in_file, in_type, out_file, out_type, globals);

  OUTPUT << "\n" << codename << " execution successful.\n";
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog)
  {
    OUTPUT << "\nUSAGE: " << prog << " in_file out_file\n";
    OUTPUT << "...or: " << prog << " command_file\n";
    OUTPUT << "       version: " << version << "\n";
    OUTPUT << "Options:\n";
    OUTPUT << "\t-directory or -d {dir} : specifies current working directory\n";
    OUTPUT << "\t-i {file} : read input and output filename data from file\n";
    OUTPUT << "\t--in_type {pamgen|generated|exodus} : set input type to the argument. Default exodus\n";
    OUTPUT << "\t--out_type {exodus} : set output type to the argument. Default exodus\n";
    OUTPUT << "\t--debug : turn on debugging output\n";
    OUTPUT << "\t--Maximum_Time {time} : maximum time from input mesh to transfer to output mesh\n";
    OUTPUT << "\t--Minimum_Time {time} : minimum time from input mesh to transfer to output mesh\n";
    OUTPUT << "\t--Surface_Split_Scheme {TOPOLOGY|ELEMENT_BLOCK|NO_SPLIT} -- how to split sidesets\n";
    Ioss::NameList db_types;
    Ioss::IOFactory::describe(&db_types);
    OUTPUT << "\nSupports database types:\n\t";
    for (Ioss::NameList::const_iterator IF = db_types.begin(); IF != db_types.end(); ++IF) {
      OUTPUT << *IF << "  ";
    }
    OUTPUT << "\n\n";
  }

  void file_copy(const std::string& inpfile, const std::string& input_type,
		 const std::string& outfile, const std::string& output_type,
		 Globals& globals)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
						    (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == NULL || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    dbi->set_surface_split_type(Ioss::int_to_surface_split(globals.surface_split_type));
    
    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    //========================================================================
    // OUTPUT ...
    //========================================================================
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(output_type, outfile, Ioss::WRITE_RESTART,
						    (MPI_Comm)MPI_COMM_WORLD);
    if (dbo == NULL || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");

    // Set the qa information...
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    if (globals.debug) OUTPUT << "DEFINING MODEL ... \n";
    if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
      OUTPUT << "ERROR: Could not put output region into define model state\n";
      std::exit(EXIT_FAILURE);
    }

    // Get all properties of input database...
    transfer_properties(&region, &output_region);
    transfer_nodeblock(region, output_region, globals.debug);
    transfer_elementblock(region, output_region, globals.debug);
    transfer_nodesets(region, output_region, globals.debug);
    transfer_sidesets(region, output_region, globals.debug);
    transfer_commsets(region, output_region, globals.debug);

    if (globals.debug) OUTPUT << "END STATE_DEFINE_MODEL... " << '\n';
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    if (globals.debug) OUTPUT << "TRANSFERRING MESH FIELD DATA ... " << '\n';
    // Model defined, now fill in the model data...
    output_region.begin_mode(Ioss::STATE_MODEL);

    // Transfer MESH field_data from input to output...
    {
      Ioss::NodeBlock *nb = region.get_node_blocks()[0];
      assert(nb != NULL);

      if (nb != NULL) {
	std::string name = nb->name();
	if (globals.debug) OUTPUT << name << ", ";

	// Find the corresponding output element_block...
	Ioss::NodeBlock *onb = output_region.get_node_block(name);
	assert(onb != NULL);
	transfer_field_data(nb, onb, Ioss::Field::MESH);
	transfer_field_data(nb, onb, Ioss::Field::ATTRIBUTE);
      }
      if (globals.debug) OUTPUT << '\n';
    }

    // Now do the same for element blocks...
    {
      Ioss::ElementBlockContainer ebs = region.get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator I = ebs.begin();

      while (I != ebs.end()) {
	std::string name = (*I)->name();
	if (globals.debug) OUTPUT << name << ", ";

	// Find the corresponding output element_block...
	Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	assert(oeb != NULL);

	if (oeb != NULL) {
	  transfer_field_data(*I, oeb, Ioss::Field::MESH);
	  transfer_field_data(*I, oeb, Ioss::Field::ATTRIBUTE);
	}
	++I;
      }
      if (globals.debug) OUTPUT << '\n';
    }

    // Node Sets
    {
      Ioss::NodeSetContainer nss = region.get_nodesets();
      Ioss::NodeSetContainer::const_iterator I = nss.begin();
      while (I != nss.end()) {
	std::string name     = (*I)->name();
	if (globals.debug) OUTPUT << name << ", ";

	// Find matching output nodeset
	Ioss::NodeSet *ons = output_region.get_nodeset(name);
	if (ons != NULL) {
	  transfer_field_data(*I, ons, Ioss::Field::MESH);
	  transfer_field_data(*I, ons, Ioss::Field::ATTRIBUTE);
	}
	++I;
      }
      if (globals.debug) OUTPUT << '\n';
    }

    // Comm Sets
    {
      Ioss::CommSetContainer css = region.get_commsets();
      Ioss::CommSetContainer::const_iterator I = css.begin();
      while (I != css.end()) {
	std::string name  = (*I)->name();
	if (globals.debug) OUTPUT << name << ", ";

	// Find matching output commset
	Ioss::CommSet *ocs = output_region.get_commset(name);
	if (ocs != NULL) {
	  transfer_field_data(*I, ocs, Ioss::Field::MESH);
	  transfer_field_data(*I, ocs, Ioss::Field::ATTRIBUTE);
	  transfer_field_data(*I, ocs, Ioss::Field::COMMUNICATION);
	}
	++I;
      }
      if (globals.debug) OUTPUT << '\n';
    }

    // Side Sets
    {
      Ioss::SideSetContainer fss = region.get_sidesets();
      Ioss::SideSetContainer::const_iterator I = fss.begin();
      while (I != fss.end()) {
	std::string name     = (*I)->name();
	if (globals.debug) OUTPUT << name << ", ";

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
	    if (globals.debug) OUTPUT << fbname << ", ";
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
      if (globals.debug) OUTPUT << '\n';
    }
    if (globals.debug) OUTPUT << "END STATE_MODEL... " << '\n';
    output_region.end_mode(Ioss::STATE_MODEL);

    if (globals.debug) OUTPUT << "DEFINING TRANSIENT FIELDS ... " << '\n';
    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      if (!globals.debug) {
	OUTPUT << "\n Number of time steps on database     =" << std::setw(9)
				<< region.get_property("state_count").get_int() << "\n\n";
      }

      output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
      {
	// For each 'TRANSIENT' field in the node blocks and element
	// blocks, transfer to the output node and element blocks.

	// Region
	{
	  transfer_fields(&region, &output_region, Ioss::Field::TRANSIENT);
	}

	// Node Blocks...
	{
	  Ioss::NodeBlockContainer ebs = region.get_node_blocks();
	  Ioss::NodeBlockContainer::const_iterator i = ebs.begin();
	  while (i != ebs.end()) {
	    std::string name      = (*i)->name();
	    if (globals.debug) OUTPUT << name << ", ";

	    // Find the corresponding output node_block...
	    Ioss::NodeBlock *oeb = output_region.get_node_block(name);
	    if (oeb != NULL) {
	      transfer_fields(*i, oeb, Ioss::Field::TRANSIENT);
	      if (globals.do_transform_fields)
		transform_fields(*i, oeb, Ioss::Field::TRANSIENT);
	    }
	    ++i;
	  }
	  if (globals.debug) OUTPUT << '\n';
	}

	// Element Blocks...
	{
	  Ioss::ElementBlockContainer ebs = region.get_element_blocks();
	  Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
	  while (i != ebs.end()) {
	    std::string name      = (*i)->name();
	    if (globals.debug) OUTPUT << name << ", ";

	    // Find the corresponding output element_block...
	    Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	    if (oeb != NULL) {
	      transfer_fields(*i, oeb, Ioss::Field::TRANSIENT);
	    }
	    ++i;
	  }
	  if (globals.debug) OUTPUT << '\n';
	}

	// Node Sets
	{
	  Ioss::NodeSetContainer nss = region.get_nodesets();
	  Ioss::NodeSetContainer::const_iterator i = nss.begin();
	  while (i != nss.end()) {
	    std::string name     = (*i)->name();
	    if (globals.debug) OUTPUT << name << ", ";

	    // Find matching output nodeset
	    Ioss::NodeSet *ons = output_region.get_nodeset(name);
	    if (ons != NULL) {
	      transfer_fields(*i, ons, Ioss::Field::TRANSIENT);
	    }
	    ++i;
	  }
	  if (globals.debug) OUTPUT << '\n';
	}

	// Side Sets
	{
	  Ioss::SideSetContainer fss = region.get_sidesets();
	  Ioss::SideSetContainer::const_iterator I = fss.begin();
	  while (I != fss.end()) {
	    std::string name     = (*I)->name();
	    if (globals.debug) OUTPUT << name << ", ";

	    // Find matching output sideset
	    Ioss::SideSet *ofs = output_region.get_sideset(name);

	    if (ofs != NULL) {
	      transfer_fields(*I, ofs, Ioss::Field::TRANSIENT);

	      Ioss::SideBlockContainer fbs = (*I)->get_side_blocks();
	      Ioss::SideBlockContainer::const_iterator J = fbs.begin();
	      while (J != fbs.end()) {

		// Find matching output sideblock
		std::string fbname = (*J)->name();
		if (globals.debug) OUTPUT << fbname << ", ";
		Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

		if (ofb != NULL) {
		  transfer_fields(*J, ofb, Ioss::Field::TRANSIENT);
		}
		++J;
	      }
	    }
	    ++I;
	  }
	  if (globals.debug) OUTPUT << '\n';
	}
      }
      if (globals.debug) OUTPUT << "END STATE_DEFINE_TRANSIENT... " << '\n';
      output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }

    if (globals.debug) OUTPUT << "TRANSFERRING TRANSIENT FIELDS ... " << '\n';
    output_region.begin_mode(Ioss::STATE_TRANSIENT);
    // Get the timesteps from the input database.  Step through them
    // and transfer fields to output database...

    int step_count = region.get_property("state_count").get_int();

    Ioss::NodeBlockContainer nbs = region.get_node_blocks();
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();

    for (int istep = 1; istep <= step_count; istep++) {
      double time = region.get_state_time(istep);
      if (time < globals.minimum_time)
	continue;
      if (globals.maximum_time != 0.0 && time > globals.maximum_time)
	break;
      
      int ostep = output_region.add_state(time);
      show_step(istep, time);

      output_region.begin_state(ostep);
      region.begin_state(istep);
      {
	transfer_field_data(&region, &output_region, Ioss::Field::TRANSIENT);
	{
	  // Output node block TRANSIENT fields...
	  Ioss::NodeBlockContainer::const_iterator I = nbs.begin();

	  while (I != nbs.end()) {
	    Ioss::NodeBlock *nb = *I;
	    ++I;
	    std::string name       = nb->name();

	    // Find the corresponding output node_block...
	    Ioss::NodeBlock *onb = output_region.get_node_block(name);
	    if (onb != NULL) {
	      transfer_field_data(nb, onb, Ioss::Field::TRANSIENT);
	      if (globals.do_transform_fields)
		transform_field_data(nb, onb, Ioss::Field::TRANSIENT);
	    }
	  }
	}

	{
	  // Output element block TRANSIENT fields...
	  Ioss::ElementBlockContainer::const_iterator I = ebs.begin();

	  while (I != ebs.end()) {
	    Ioss::ElementBlock *eb = *I;
	    ++I;
	    std::string name       = eb->name();

	    // Find the corresponding output element_block...
	    Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	    if (oeb != NULL) {
	      transfer_field_data(eb, oeb, Ioss::Field::TRANSIENT);
	    }
	  }
	}

	// Node Sets
	{
	  Ioss::NodeSetContainer nss = region.get_nodesets();
	  Ioss::NodeSetContainer::const_iterator I = nss.begin();
	  while (I != nss.end()) {
	    std::string name     = (*I)->name();
	    if (globals.debug) OUTPUT << name << ", ";

	    // Find matching output nodeset
	    Ioss::NodeSet *ons = output_region.get_nodeset(name);
	    if (ons != NULL) {
	      transfer_field_data(*I, ons, Ioss::Field::TRANSIENT);
	    }
	    ++I;
	  }
	  if (globals.debug) OUTPUT << '\n';
	}

	// Side Sets
	{
	  Ioss::SideSetContainer fss = region.get_sidesets();
	  Ioss::SideSetContainer::const_iterator I = fss.begin();
	  while (I != fss.end()) {
	    std::string name     = (*I)->name();
	    if (globals.debug) OUTPUT << name << ", ";
	  
	    // Find matching output sideset
	    Ioss::SideSet *ofs = output_region.get_sideset(name);
	  
	    if (ofs != NULL) {
	      transfer_field_data(*I, ofs, Ioss::Field::TRANSIENT);
	    
	      Ioss::SideBlockContainer fbs = (*I)->get_side_blocks();
	      Ioss::SideBlockContainer::const_iterator J = fbs.begin();
	      while (J != fbs.end()) {
	      
		// Find matching output sideblock
		std::string fbname = (*J)->name();
		if (globals.debug) OUTPUT << fbname << ", ";
		Ioss::SideBlock *ofb = ofs->get_side_block(fbname);
	      
		if (ofb != NULL) {
		  transfer_field_data(*J, ofb, Ioss::Field::TRANSIENT, "", false);
		}
		++J;
	      }
	    }
	    ++I;
	  }
	}
      }
      region.end_state(istep);
      output_region.end_state(ostep);
    }
    if (globals.debug) OUTPUT << "END STATE_TRANSIENT... " << '\n';
    output_region.end_mode(Ioss::STATE_TRANSIENT);
  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeBlockContainer    nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i = nbs.begin();
    int id = 1;
    while (i != nbs.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
      int    num_nodes = (*i)->get_property("entity_count").get_int();
      int    degree    = (*i)->get_property("component_degree").get_int();
      if (!debug) {
	OUTPUT << " Number of coordinates per node       =" << std::setw(9) << degree << "\n";
	OUTPUT << " Number of nodes                      =" << std::setw(9) << num_nodes << "\n";
      }

      Ioss::NodeBlock *nb = new Ioss::NodeBlock(output_region.get_database(), name, num_nodes, degree);
      output_region.add(nb);

      transfer_properties(*i, nb);
      transfer_fields(*i, nb, Ioss::Field::MESH);
      transfer_fields(*i, nb, Ioss::Field::ATTRIBUTE);
      ++i;
      ++id;
    }
    if (debug) OUTPUT << '\n';
  }

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    int total_elements = 0;
    while (i != ebs.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
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
      OUTPUT << " Number of elements                   =" << std::setw(9) << total_elements << "\n";
      OUTPUT << " Number of element blocks             =" << std::setw(9) << ebs.size() << "\n\n";
    } else {
      OUTPUT << '\n';
    }
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::SideSetContainer      fss = region.get_sidesets();
    Ioss::SideSetContainer::const_iterator i = fss.begin();
    int total_sides = 0;
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
      OUTPUT << " Number of element side sets          =" << std::setw(9) << fss.size() << "\n";
      OUTPUT << "     Number of element sides          =" << std::setw(9) << total_sides << "\n";
    } else {
      OUTPUT << '\n';
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeSetContainer      nss = region.get_nodesets();
    Ioss::NodeSetContainer::const_iterator i = nss.begin();
    int total_nodes = 0;
    while (i != nss.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
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
      OUTPUT << " Number of nodal point sets           =" << std::setw(9) << nss.size() << "\n";
      OUTPUT << "     Length of node list              =" << std::setw(9) << total_nodes << "\n";
    } else {
      OUTPUT << '\n';
    }
  }

  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::CommSetContainer      css = region.get_commsets();
    Ioss::CommSetContainer::const_iterator i = css.begin();
    while (i != css.end()) {
      std::string name      = (*i)->name();
      if (debug) OUTPUT << name << ", ";
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
      if (field_name != "ids" && !oge->field_exists(field_name) &&
	  (prefix.length() == 0 || std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
	// If the field does not already exist, add it to the output node block
	Ioss::Field field = ige->get_field(field_name);
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

      int isize = ige->get_field(field_name).get_size();
      int osize = oge->get_field(out_field_name).get_size();
      assert (isize == osize);

      data.resize(isize);
      ige->get_field_data(field_name,     &data[0], isize);
      oge->put_field_data(out_field_name, &data[0], osize);
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
}
