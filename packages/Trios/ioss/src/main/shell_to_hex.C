/*--------------------------------------------------------------------*/
/*    Copyright 2000, 2008, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_CodeTypes.h>

#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <cstring>

#ifndef IOSS_STANDALONE
#include <Slib_Env.h>
#endif

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#ifndef NO_XDMF_SUPPORT
#include <xdmf/Ioxf_Initializer.h>
#endif

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    bool debug;
    bool do_normals;
    bool reverse_normals;
    double thickness;
  };

  void show_usage(const std::string &prog);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_properties(Ioss::GroupingEntity *ige,
			   Ioss::GroupingEntity *oge);

  void output_normals(Ioss::Region &region, Ioss::Region &output_region,
		      bool reverse_normals, double thickness);
  void calculate_normals(std::vector<double> &node_normal,
			 int num_elem, int num_node_per_elem,
			 std::vector<double> &coord,
			 std::vector<int>    &conn,
			 bool reverse_normals);

  void file_copy(const std::string& inpfile, const std::string& input_type,
		 const std::string& outfile, const std::string& output_type,
		 Globals& globals);
}
// ========================================================================

#ifdef IOSS_STANDALONE
// Global variables, yes, global...
std::string codename;
std::string version = "$Revision$";
#else
namespace {
  std::string codename;
  std::string version = "$Revision$";
}
#endif

#ifndef IOSS_STANDALONE
namespace {

void bootstrap()
{
  // Add my command line options to the option descriptions.
  boost::program_options::options_description desc("Use case options");
  desc.add_options()
    ("input-deck,i",      boost::program_options::value<std::string>(), "Analysis input file")
    ("restart-time,r",    boost::program_options::value<std::string>(), "Restart time")
    ("parser-database,p", boost::program_options::value<std::string>(), "Parser database")
    ("Calculate_Normals", "Calculate surface normal fields")
    ("Reverse_Normals",   "Reverse Normals")
    ("thickness",         boost::program_options::value<std::string>(), "Thickness of generated hex elements");
  
  stk::get_options_description().add(desc);
}

stk::Bootstrap x(&bootstrap);

} // namespace <unnamed>
#endif

int main(int argc, char *argv[])
{
  std::string in_type = "exodusII";
  std::string out_type = "exodusII";
  std::string ss_type = "exodusII";

  Globals globals;

  globals.do_normals = false;
  globals.reverse_normals = false;
  globals.thickness = 1.0;

  codename = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos)
    codename = codename.substr(ind+1, codename.size());

  // Check the program name to see if of the form 'exosaf' or 'safexo'
  // and if it is, set the in_type and out_type accordingly...
  if (std::strncmp(codename.c_str(), "shell_to_hex", 12) == 0) {
    codename = "shell_to_hex";
    in_type  = "exodusII";
    out_type = "exodusII";
    globals.do_normals = true;
  }
  else {
    codename = "SIERRA";
  }

#ifndef IOSS_STANDALONE
  sierra::Env::Startup startup__(&argc, &argv, codename.c_str(), __DATE__ " " __TIME__); //, opts);
#endif

  Ioss::Init::Initializer io;
#ifndef NO_XDMF_SUPPORT
  Ioxf::Initializer ioxf;
#endif

#ifdef IOSS_STANDALONE
  globals.debug = false;
#else
  globals.debug = true;
#endif

  // Skip past any options...
  int i=1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("--reverse", argv[i]) == 0) {
      globals.reverse_normals = true;
      i++;
    }
    else if (std::strcmp("--Reverse_Normals", argv[i]) == 0) {
      globals.reverse_normals = true;
      i++;
    }
    else if (std::strcmp("--Calculate_Normals", argv[i]) == 0) {
      globals.do_normals = true;
      i++;
    }
    else if (std::strcmp("--thickness", argv[i]) == 0) {
      i++;
      globals.thickness = std::strtod(argv[i++], NULL);
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

  if (argc - i == 2) {
    in_file   = Ioss::Utils::local_filename(argv[i++], in_type);
    out_file  = Ioss::Utils::local_filename(argv[i++], out_type);
  }
  else if (argc - i == 1) {
    std::string input_file = Ioss::Utils::local_filename(argv[i++], "text");

    std::ifstream input(input_file.c_str());
    if (!input) {
      std::cerr  << "Error opening file '" << input_file << "'.\n";
      show_usage(codename);
      return (EXIT_FAILURE);
    }
    // First line should be input_file_name input_file_type
    // Last line should be output_file_name output_file_type
    std::string tmp;
    input >> tmp >> in_type;
    in_file = Ioss::Utils::local_filename(tmp, in_type);

    input >> tmp >> out_type;
    out_file = Ioss::Utils::local_filename(tmp, out_type);
  }
#ifndef IOSS_STANDALONE
  else if (argc <= i) {
    const std::string &calc_normals_option    = sierra::Env::get_param("Calculate_Normals");
    const std::string &reverse_normals_option = sierra::Env::get_param("Reverse_Normals");
    const std::string &thickness = sierra::Env::get_param("thickness");
    globals.thickness = std::strtod(thickness.c_str(), NULL);

    if (!calc_normals_option.empty())
      globals.do_normals = true;
    if (!reverse_normals_option.empty())
      globals.reverse_normals = true;

    std::string input_file(sierra::Env::get_param("input-deck"));
    if (input_file == "") {
      std::cerr  << "Error: No input file specified\n";
      show_usage(codename);
      return (EXIT_FAILURE);
    }
    input_file = Ioss::Utils::local_filename(input_file, "text");
    std::ifstream input(input_file.c_str());
    if (!input) {
      std::cerr  << "Error opening file '" << input_file  << "'.\n";
      show_usage(codename);
      return (EXIT_FAILURE);
    }
    // First line should be input_file_name input_file_type
    // Second line should be output_file_name output_file_type

    std::string tmp;
    input >> tmp >> in_type;
    in_file = Ioss::Utils::local_filename(tmp, in_type);

    input >> tmp >> out_type;
    out_file = Ioss::Utils::local_filename(tmp, out_type);
  }
#else
  else {
    show_usage(codename);
    return(EXIT_FAILURE);
  }
#endif

  std::cerr  << "Input:    '" << in_file  << "', Type: " << in_type  << '\n';
  std::cerr  << "Output:   '" << out_file << "', Type: " << out_type << '\n';
  if (globals.reverse_normals)
    std::cerr << "Reversing Normals\n";
  std::cerr << "Thickness:  " << globals.thickness << "\n";
  std::cerr  << '\n';

  if (globals.do_normals)
    globals.debug = false;

  file_copy(in_file, in_type, out_file, out_type, globals);

  std::cerr  << "\n" << codename << " execution successful.\n";
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog)
  {
    std::cerr  << "USAGE: " << prog << " in_file out_file\n";

    std::cerr  << "...or: " << prog << " command_file\n";
    std::cerr  << "       version: " << version << "\n";
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
    if (dbi == NULL || !dbi->ok()) {
      std::cerr  << "ERROR: Could not open database '" << inpfile
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
      std::cerr  << "ERROR: Could not create output database '" << outfile
		 << "' of type '" << output_type << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");

    // Set the qa information...
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    if (globals.debug) std::cerr  << "DEFINING MODEL ... \n";
    if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
      std::cerr  << "ERROR: Could not put output region into define model state\n";
      std::exit(EXIT_FAILURE);
    }

    // Get all properties of input database...
    transfer_properties(&region, &output_region);
    transfer_nodeblock(region, output_region, globals.debug);
    transfer_elementblock(region, output_region, globals.debug);

    if (globals.debug) std::cerr  << "END STATE_DEFINE_MODEL... " << '\n';
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    if (globals.debug) std::cerr  << "TRANSFERRING MESH FIELD DATA ... " << '\n';
    // Model defined, now fill in the model data...
    output_region.begin_mode(Ioss::STATE_MODEL);
    output_normals(region, output_region, globals.reverse_normals, globals.thickness);
    output_region.end_mode(Ioss::STATE_MODEL);

  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::NodeBlockContainer    nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i = nbs.begin();
    int id = 1;
    while (i != nbs.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr  << name << ", ";
      int    num_nodes = (*i)->get_property("entity_count").get_int();
      int    degree    = (*i)->get_property("component_degree").get_int();
      if (!debug) {
	std::cerr  << " Number of coordinates per node       =" << std::setw(9) << degree << "\n";
	std::cerr  << " Number of nodes                      =" << std::setw(9) << num_nodes << "\n";
      }

      Ioss::NodeBlock *nb = new Ioss::NodeBlock(output_region.get_database(), name, 2*num_nodes, degree);
      output_region.add(nb);
      ++i;
      ++id;
    }
    if (debug) std::cerr  << '\n';
  }

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    int total_elements = 0;
    while (i != ebs.end()) {
      std::string name      = (*i)->name();
      if (debug) std::cerr  << name << ", ";
      int    num_elem  = (*i)->get_property("entity_count").get_int();
      int    num_attrib= 0;
      total_elements += num_elem;

      std::string type;
      int num_node_per_elem = (*i)->topology()->number_nodes();
      if (num_node_per_elem == 4)
	type = "hex";
      else
	type = "wedge";
      
      Ioss::ElementBlock *eb = new Ioss::ElementBlock(output_region.get_database(), name, type,
						      num_elem, num_attrib);
      output_region.add(eb);
      ++i;
    }
    if (!debug) {
      std::cerr  << " Number of elements                   =" << std::setw(9) << total_elements << "\n";
      std::cerr  << " Number of element blocks             =" << std::setw(9) << ebs.size() << "\n\n";
    } else {
      std::cerr  << '\n';
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
      if (!oge->property_exists(*I))
	oge->property_add(ige->get_property(*I));
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

  void output_normals(Ioss::Region &region, Ioss::Region &output_region, bool reverse_normals,
		      double thickness)
  {
    Ioss::NodeBlock *nb  = (*region.get_node_blocks().begin());
    Ioss::NodeBlock *nbo = (*output_region.get_node_blocks().begin());

    // Get the nodal coordinates...
    int num_nodes  = nb->get_property("entity_count").get_int();

    {
      std::vector<int>    ids(2*num_nodes);
      for (int i=0; i < num_nodes; i++) {
	ids[i] = i+1;
	ids[num_nodes+i] = num_nodes+i+1;
      }
      nbo->put_field_data("ids", ids);
    }
    
    std::vector<double> coord(3*num_nodes);
    nb->get_field_data("mesh_model_coordinates", coord);

    // Also get an array for the average nodal normal vector...
    std::vector<double> node_normal(3*num_nodes);
    std::fill(node_normal.begin(), node_normal.end(), 0.0);

    // Iterate over the element blocks and calculate node normals
    std::vector<int>    conn;
    std::vector<int>    output_conn;
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer out_ebs = output_region.get_element_blocks();

    Ioss::ElementBlockContainer::const_iterator ib  = ebs.begin();
    Ioss::ElementBlockContainer::const_iterator out_ib  = out_ebs.begin();
    while (ib != ebs.end() && out_ib != out_ebs.end()) {
      Ioss::ElementBlock *eb = *ib; ++ib;
      Ioss::ElementBlock *out_eb = *out_ib; ++out_ib;
      std::string name = (*eb).name();

      int num_elem  = eb->get_property("entity_count").get_int();
      int num_node_per_elem = eb->topology()->number_nodes();

      // Get the connectivity array...
      conn.resize(num_elem * num_node_per_elem);
      output_conn.resize(2*conn.size());

      std::vector<int> eids(num_elem);
      eb->get_field_data("ids", eids);
      out_eb->put_field_data("ids", eids);

      eb->get_field_data("connectivity", conn);

      // Connectivity is in global id space; change to local...
      for (int i=0; i < num_elem * num_node_per_elem; i++) {
	int local = region.node_global_to_local(conn[i]);
	conn[i] = local-1;
      }
      
      calculate_normals(node_normal, num_elem, num_node_per_elem, coord, conn,
			reverse_normals);


      assert(num_node_per_elem == 3 || num_node_per_elem == 4);
      if (reverse_normals) {
	for (int i=0; i < num_elem; i++) {
	  if (num_node_per_elem == 4) {
	    output_conn[8*i+0] = 1+conn[4*i+0] + num_nodes;
	    output_conn[8*i+1] = 1+conn[4*i+1] + num_nodes;
	    output_conn[8*i+2] = 1+conn[4*i+2] + num_nodes;
	    output_conn[8*i+3] = 1+conn[4*i+3] + num_nodes;

	    output_conn[8*i+4] = 1+conn[4*i+0];
	    output_conn[8*i+5] = 1+conn[4*i+1];
	    output_conn[8*i+6] = 1+conn[4*i+2];
	    output_conn[8*i+7] = 1+conn[4*i+3];
	  } else {
	    output_conn[6*i+0] = 1+conn[3*i+0] + num_nodes;
	    output_conn[6*i+1] = 1+conn[3*i+1] + num_nodes;
	    output_conn[6*i+2] = 1+conn[3*i+2] + num_nodes;

	    output_conn[6*i+4] = 1+conn[3*i+0];
	    output_conn[6*i+5] = 1+conn[3*i+1];
	    output_conn[6*i+6] = 1+conn[3*i+2];
	  }
	}
      } else {
	for (int i=0; i < num_elem; i++) {
	  if (num_node_per_elem == 4) {
	    output_conn[8*i+0] = 1+conn[4*i+0];
	    output_conn[8*i+1] = 1+conn[4*i+1];
	    output_conn[8*i+2] = 1+conn[4*i+2];
	    output_conn[8*i+3] = 1+conn[4*i+3];

	    output_conn[8*i+4] = 1+conn[4*i+0] + num_nodes;
	    output_conn[8*i+5] = 1+conn[4*i+1] + num_nodes;
	    output_conn[8*i+6] = 1+conn[4*i+2] + num_nodes;
	    output_conn[8*i+7] = 1+conn[4*i+3] + num_nodes;
	  } else {
	    output_conn[6*i+0] = 1+conn[3*i+0];
	    output_conn[6*i+1] = 1+conn[3*i+1];
	    output_conn[6*i+2] = 1+conn[3*i+2];

	    output_conn[6*i+4] = 1+conn[3*i+0] + num_nodes;
	    output_conn[6*i+5] = 1+conn[3*i+1] + num_nodes;
	    output_conn[6*i+6] = 1+conn[3*i+2] + num_nodes;
	  }
	}
      }
      out_eb->put_field_data("connectivity", output_conn);
    }

    int nsize = node_normal.size();
    for (int i=0; i < nsize; i+=3) {
      Vector3 a(node_normal[i+0], node_normal[i+1], node_normal[i+2]);
      a.normalize();

      node_normal[i+0] = a.x;
      node_normal[i+1] = a.y;
      node_normal[i+2] = a.z;
    }

    // The output created nodes of the new hexes are simply the input nodes
    // translated along the nodal vector the specfied distance.  The node id is
    // simply the input node id + number_of_input_nodes.
    std::vector<double> output_coord(2*num_nodes*3);
    for (int i=0; i < num_nodes; i++) {
      output_coord[3*i+0] = coord[3*i+0];
      output_coord[3*i+1] = coord[3*i+1];
      output_coord[3*i+2] = coord[3*i+2];

      output_coord[3*(num_nodes+i)+0] = coord[3*i+0] + thickness * node_normal[3*i+0];
      output_coord[3*(num_nodes+i)+1] = coord[3*i+1] + thickness * node_normal[3*i+1];
      output_coord[3*(num_nodes+i)+2] = coord[3*i+2] + thickness * node_normal[3*i+2];
    }
    nbo->put_field_data("mesh_model_coordinates", output_coord);
  }

  void calculate_normals(std::vector<double> &node_normal,
			 int num_elem, int num_node_per_elem,
			 std::vector<double> &coord,
			 std::vector<int>    &conn,
			 bool reverse_normals)
  {
    // Iterate the connectivity array and calculate normals...
    // The elements should all be shells with the outward normal
    // defined by the node ordering 1.2.3.

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

	  int node = conn[ioff+i];
	  node_normal[node*3+0] += a.x;
	  node_normal[node*3+1] += a.y;
	  node_normal[node*3+2] += a.z;
	}
      }
    }
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
