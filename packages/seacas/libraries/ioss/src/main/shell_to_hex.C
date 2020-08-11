// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_CodeTypes.h>

#include <Ionit_Initializer.h>
#include <Ioss_Hex8.h>
#include <Ioss_Wedge6.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_State.h"
#include "Ioss_Utils.h"
#include "vector3d.h"

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    bool        debug{};
    bool        do_normals{};
    bool        reverse_normals{};
    double      thickness{};
    std::string working_directory{};
  };

  void show_usage(const std::string &prog);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge);

  void output_normals(Ioss::Region &region, Ioss::Region &output_region, bool reverse_normals,
                      double thickness);
  void calculate_normals(std::vector<double> &node_normal, int num_elem, int num_node_per_elem,
                         std::vector<double> &coord, std::vector<int> &conn, bool reverse_normals);

  void file_copy(const std::string &inpfile, const std::string &input_type,
                 const std::string &outfile, const std::string &output_type, Globals &globals);
} // namespace
// ========================================================================

namespace {
  std::string codename;
  std::string version = "0.9";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  std::string in_type  = "exodusII";
  std::string out_type = "exodusII";

  Globals globals;

  globals.do_normals      = false;
  globals.reverse_normals = false;
  globals.thickness       = 1.0;

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  // Check the program name to see if of the form 'exosaf' or 'safexo'
  // and if it is, set the in_type and out_type accordingly...
  if (Ioss::Utils::str_equal(codename, "shell_to_hex")) {
    codename           = "shell_to_hex";
    in_type            = "exodusII";
    out_type           = "exodusII";
    globals.do_normals = true;
  }
  else {
    codename = "SIERRA";
  }

  Ioss::Init::Initializer io;

  globals.debug = false;

  // Skip past any options...
  int i = 1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("-directory", argv[i]) == 0 || std::strcmp("--directory", argv[i]) == 0 ||
        std::strcmp("-d", argv[i]) == 0) {
      i++;
      globals.working_directory = argv[i++];
    }
    else if (std::strcmp("--reverse", argv[i]) == 0) {
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
      globals.thickness = std::strtod(argv[i++], nullptr);
    }

    // Found an option.  See if it has an argument...
    else if (i + 1 < argc && argv[i + 1][0] == '-') {
      // No argument, another option
      i++;
    }
    else {
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
    in_file  = Ioss::Utils::local_filename(argv[i++], in_type, globals.working_directory);
    out_file = Ioss::Utils::local_filename(argv[i++], out_type, globals.working_directory);
  }
  else if (argc - i == 1) {
    std::string input_file =
        Ioss::Utils::local_filename(argv[i++], "text", globals.working_directory);

    std::ifstream input(input_file);
    if (!input) {
      std::cerr << "Error opening file '" << input_file << "'.\n";
      show_usage(codename);
      return (EXIT_FAILURE);
    }
    // First line should be input_file_name input_file_type
    // Last line should be output_file_name output_file_type
    std::string tmp;
    input >> tmp >> in_type;
    in_file = Ioss::Utils::local_filename(tmp, in_type, globals.working_directory);

    input >> tmp >> out_type;
    out_file = Ioss::Utils::local_filename(tmp, out_type, globals.working_directory);
  }
  else {
    show_usage(codename);
    return (EXIT_FAILURE);
  }

  std::cerr << "Input:    '" << in_file << "', Type: " << in_type << '\n';
  std::cerr << "Output:   '" << out_file << "', Type: " << out_type << '\n';
  if (globals.reverse_normals) {
    std::cerr << "Reversing Normals\n";
  }
  std::cerr << "Thickness:  " << globals.thickness << "\n";
  std::cerr << '\n';

  if (globals.do_normals) {
    globals.debug = false;
  }

  file_copy(in_file, in_type, out_file, out_type, globals);

  std::cerr << "\n" << codename << " execution successful.\n";
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog)
  {
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

  void file_copy(const std::string &inpfile, const std::string &input_type,
                 const std::string &outfile, const std::string &output_type, Globals &globals)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi =
        Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART, (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == nullptr || !dbi->ok()) {
      std::cerr << "ERROR: Could not open database '" << inpfile << "' of type '" << input_type
                << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    //========================================================================
    // OUTPUT ...
    //========================================================================
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(output_type, outfile, Ioss::WRITE_RESTART,
                                                    (MPI_Comm)MPI_COMM_WORLD);
    if (dbo == nullptr || !dbo->ok()) {
      std::cerr << "ERROR: Could not create output database '" << outfile << "' of type '"
                << output_type << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");

    // Set the qa information...
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    if (globals.debug) {
      std::cerr << "DEFINING MODEL ... \n";
    }
    if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
      std::cerr << "ERROR: Could not put output region into define model state\n";
      std::exit(EXIT_FAILURE);
    }

    // Get all properties of input database...
    transfer_properties(&region, &output_region);
    transfer_nodeblock(region, output_region, globals.debug);
    transfer_elementblock(region, output_region, globals.debug);

    if (globals.debug) {
      std::cerr << "END STATE_DEFINE_MODEL... " << '\n';
    }
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    if (globals.debug) {
      std::cerr << "TRANSFERRING MESH FIELD DATA ... " << '\n';
      // Model defined, now fill in the model data...
    }
    output_region.begin_mode(Ioss::STATE_MODEL);
    output_normals(region, output_region, globals.reverse_normals, globals.thickness);
    output_region.end_mode(Ioss::STATE_MODEL);
  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const Ioss::NodeBlockContainer &         nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i   = nbs.begin();
    while (i != nbs.end()) {
      const std::string &name = (*i)->name();
      if (debug) {
        std::cerr << name << ", ";
      }
      int num_nodes = (*i)->entity_count();
      int degree    = (*i)->get_property("component_degree").get_int();
      if (!debug) {
        std::cerr << " Number of coordinates per node       =" << std::setw(9) << degree << "\n";
        std::cerr << " Number of nodes                      =" << std::setw(9) << num_nodes << "\n";
      }

      auto nb = new Ioss::NodeBlock(output_region.get_database(), name, 2 * num_nodes, degree);
      output_region.add(nb);
      ++i;
    }
    if (debug) {
      std::cerr << '\n';
    }
  }

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const Ioss::ElementBlockContainer &         ebs            = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i              = ebs.begin();
    int                                         total_elements = 0;
    while (i != ebs.end()) {
      const std::string &name = (*i)->name();
      if (debug) {
        std::cerr << name << ", ";
      }
      int num_elem = (*i)->entity_count();
      total_elements += num_elem;

      std::string type;
      int         num_node_per_elem = (*i)->topology()->number_nodes();
      if (num_node_per_elem == 4) {
        type = Ioss::Hex8::name;
      }
      else {
        type = Ioss::Wedge6::name;
      }
      auto eb = new Ioss::ElementBlock(output_region.get_database(), name, type, num_elem);
      output_region.add(eb);
      ++i;
    }
    if (!debug) {
      std::cerr << " Number of elements                   =" << std::setw(9) << total_elements
                << "\n";
      std::cerr << " Number of element blocks             =" << std::setw(9) << ebs.size()
                << "\n\n";
    }
    else {
      std::cerr << '\n';
    }
  }

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge)
  {
    Ioss::NameList names;
    ige->property_describe(&names);

    // Iterate through properties and transfer to output database...
    Ioss::NameList::const_iterator I;
    for (I = names.begin(); I != names.end(); ++I) {
      if (!oge->property_exists(*I)) {
        oge->property_add(ige->get_property(*I));
      }
    }
  }

  void output_normals(Ioss::Region &region, Ioss::Region &output_region, bool reverse_normals,
                      double thickness)
  {
    Ioss::NodeBlock *nb  = (*region.get_node_blocks().begin());
    Ioss::NodeBlock *nbo = (*output_region.get_node_blocks().begin());

    // Get the nodal coordinates...
    int num_nodes = nb->entity_count();

    {
      std::vector<int> ids(2 * num_nodes);
      for (int i = 0; i < num_nodes; i++) {
        ids[i]             = i + 1;
        ids[num_nodes + i] = num_nodes + i + 1;
      }
      nbo->put_field_data("ids", ids);
    }

    std::vector<double> coord(3 * num_nodes);
    nb->get_field_data("mesh_model_coordinates", coord);

    // Also get an array for the average nodal normal vector...
    std::vector<double> node_normal(3 * num_nodes);
    std::fill(node_normal.begin(), node_normal.end(), 0.0);

    // Iterate over the element blocks and calculate node normals
    std::vector<int>                   conn;
    std::vector<int>                   output_conn;
    const Ioss::ElementBlockContainer &ebs     = region.get_element_blocks();
    const Ioss::ElementBlockContainer &out_ebs = output_region.get_element_blocks();

    Ioss::ElementBlockContainer::const_iterator ib     = ebs.begin();
    Ioss::ElementBlockContainer::const_iterator out_ib = out_ebs.begin();
    while (ib != ebs.end() && out_ib != out_ebs.end()) {
      Ioss::ElementBlock *eb = *ib;
      ++ib;
      Ioss::ElementBlock *out_eb = *out_ib;
      ++out_ib;

      int num_elem          = eb->entity_count();
      int num_node_per_elem = eb->topology()->number_nodes();

      // Get the connectivity array...
      conn.resize(num_elem * num_node_per_elem);
      output_conn.resize(2 * conn.size());

      std::vector<int> eids(num_elem);
      eb->get_field_data("ids", eids);
      out_eb->put_field_data("ids", eids);

      eb->get_field_data("connectivity", conn);

      // Connectivity is in global id space; change to local...
      for (int i = 0; i < num_elem * num_node_per_elem; i++) {
        int local = region.node_global_to_local(conn[i]);
        conn[i]   = local - 1;
      }

      calculate_normals(node_normal, num_elem, num_node_per_elem, coord, conn, reverse_normals);

      assert(num_node_per_elem == 3 || num_node_per_elem == 4);
      if (reverse_normals) {
        for (int i = 0; i < num_elem; i++) {
          if (num_node_per_elem == 4) {
            output_conn[8 * i + 0] = 1 + conn[4 * i + 0] + num_nodes;
            output_conn[8 * i + 1] = 1 + conn[4 * i + 1] + num_nodes;
            output_conn[8 * i + 2] = 1 + conn[4 * i + 2] + num_nodes;
            output_conn[8 * i + 3] = 1 + conn[4 * i + 3] + num_nodes;

            output_conn[8 * i + 4] = 1 + conn[4 * i + 0];
            output_conn[8 * i + 5] = 1 + conn[4 * i + 1];
            output_conn[8 * i + 6] = 1 + conn[4 * i + 2];
            output_conn[8 * i + 7] = 1 + conn[4 * i + 3];
          }
          else {
            output_conn[6 * i + 0] = 1 + conn[3 * i + 0] + num_nodes;
            output_conn[6 * i + 1] = 1 + conn[3 * i + 1] + num_nodes;
            output_conn[6 * i + 2] = 1 + conn[3 * i + 2] + num_nodes;

            output_conn[6 * i + 4] = 1 + conn[3 * i + 0];
            output_conn[6 * i + 5] = 1 + conn[3 * i + 1];
            output_conn[6 * i + 6] = 1 + conn[3 * i + 2];
          }
        }
      }
      else {
        for (int i = 0; i < num_elem; i++) {
          if (num_node_per_elem == 4) {
            output_conn[8 * i + 0] = 1 + conn[4 * i + 0];
            output_conn[8 * i + 1] = 1 + conn[4 * i + 1];
            output_conn[8 * i + 2] = 1 + conn[4 * i + 2];
            output_conn[8 * i + 3] = 1 + conn[4 * i + 3];

            output_conn[8 * i + 4] = 1 + conn[4 * i + 0] + num_nodes;
            output_conn[8 * i + 5] = 1 + conn[4 * i + 1] + num_nodes;
            output_conn[8 * i + 6] = 1 + conn[4 * i + 2] + num_nodes;
            output_conn[8 * i + 7] = 1 + conn[4 * i + 3] + num_nodes;
          }
          else {
            output_conn[6 * i + 0] = 1 + conn[3 * i + 0];
            output_conn[6 * i + 1] = 1 + conn[3 * i + 1];
            output_conn[6 * i + 2] = 1 + conn[3 * i + 2];

            output_conn[6 * i + 4] = 1 + conn[3 * i + 0] + num_nodes;
            output_conn[6 * i + 5] = 1 + conn[3 * i + 1] + num_nodes;
            output_conn[6 * i + 6] = 1 + conn[3 * i + 2] + num_nodes;
          }
        }
      }
      out_eb->put_field_data("connectivity", output_conn);
    }

    int nsize = node_normal.size();
    for (int i = 0; i < nsize; i += 3) {
      vector3d a(node_normal[i + 0], node_normal[i + 1], node_normal[i + 2]);
      a.normalize();

      node_normal[i + 0] = a.x;
      node_normal[i + 1] = a.y;
      node_normal[i + 2] = a.z;
    }

    // The output created nodes of the new hexes are simply the input nodes
    // translated along the nodal vector the specified distance.  The node id is
    // simply the input node id + number_of_input_nodes.
    std::vector<double> output_coord(2 * num_nodes * 3);
    for (int i = 0; i < num_nodes; i++) {
      output_coord[3 * i + 0] = coord[3 * i + 0];
      output_coord[3 * i + 1] = coord[3 * i + 1];
      output_coord[3 * i + 2] = coord[3 * i + 2];

      output_coord[3 * (num_nodes + i) + 0] = coord[3 * i + 0] + thickness * node_normal[3 * i + 0];
      output_coord[3 * (num_nodes + i) + 1] = coord[3 * i + 1] + thickness * node_normal[3 * i + 1];
      output_coord[3 * (num_nodes + i) + 2] = coord[3 * i + 2] + thickness * node_normal[3 * i + 2];
    }
    nbo->put_field_data("mesh_model_coordinates", output_coord);
  }

  void calculate_normals(std::vector<double> &node_normal, int num_elem, int num_node_per_elem,
                         std::vector<double> &coord, std::vector<int> &conn, bool reverse_normals)
  {
    // Iterate the connectivity array and calculate normals...
    // The elements should all be shells with the outward normal
    // defined by the node ordering 1.2.3.

    for (int iel = 0; iel < num_elem; iel++) {
      int ioff = iel * num_node_per_elem;

      // Triangular faces...
      if (num_node_per_elem == 3) {
        vector3d local[3];
        for (int i = 0; i < 3; i++) {
          int node = conn[ioff + i];
          local[i].set(coord[node * 3 + 0], coord[node * 3 + 1], coord[node * 3 + 2]);
        }

        vector3d plnorm = vector3d::plane_normal(local[0], local[1], local[2]);
        plnorm.normalize();
        if (reverse_normals) {
          plnorm.reverse();
        }

        for (int i = 0; i < 3; i++) {
          int node = conn[ioff + i];
          node_normal[node * 3 + 0] += plnorm.x;
          node_normal[node * 3 + 1] += plnorm.y;
          node_normal[node * 3 + 2] += plnorm.z;
        }
      }
      else {

        // Quadrilateral faces...
        assert(num_node_per_elem == 4);
        vector3d local[4];
        for (int i = 0; i < 4; i++) {
          int node = conn[ioff + i];
          local[i].set(&coord[node * 3]);
        }

        for (int i = 0; i < 4; i++) {
          // at node 0 -- vector from 3-0 X 0-1
          // at node 1 -- vector from 0-1 X 1-2
          // at node 2 -- vector from 1-2 X 2-3
          // at node 3 -- vector from 2-3 X 3-0
          int nb = (i + 3) % 4;
          int na = (i + 1) % 4;

          vector3d a = vector3d::plane_normal(local[nb], local[i], local[na]);
          a.normalize();
          if (reverse_normals) {
            a.reverse();
          }

          int node = conn[ioff + i];
          node_normal[node * 3 + 0] += a.x;
          node_normal[node * 3 + 1] += a.y;
          node_normal[node * 3 + 2] += a.z;
        }
      }
    }
  }
} // namespace
