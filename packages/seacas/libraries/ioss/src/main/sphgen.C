// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_Hex8.h>
#include <Ioss_Utils.h>
#include <Ioss_Wedge6.h>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_Field.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_State.h"

#define DO_OUTPUT std::cerr

// ========================================================================

namespace {
  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    std::string working_directory{};
    double      scale_factor{};
  };

  void hex_volume(Ioss::ElementBlock *block, const std::vector<double> &coordinates,
                  double *centroids, std::vector<double> &volume, std::vector<double> &radius,
                  double scale_factor);

  void show_usage(const std::string &prog);

  void create_sph(const std::string &inpfile, const std::string &input_type,
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

  std::string in_type  = "exodus";
  std::string out_type = "exodus";

  Globals globals;

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  Ioss::Init::Initializer io;

  std::string input_file;
  globals.scale_factor = 8.0;

  // Skip past any options...
  int i = 1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("-directory", argv[i]) == 0 || std::strcmp("--directory", argv[i]) == 0 ||
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
    else if (std::strcmp("-i", argv[i]) == 0 || std::strcmp("--input", argv[i]) == 0) {
      i++;
      input_file = argv[i++];
    }
    else if (std::strcmp("--scale_factor", argv[i]) == 0) {
      i++;
      globals.scale_factor = std::strtod(argv[i++], nullptr);
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

  // If two or more arguments are specified at end of command line, they are
  // the input and output files to be converted.
  // The file types are assumed to be as 'hardwired' above...

  std::string cwd = globals.working_directory;

  if (input_file.empty()) {
    if (argc - i == 2) {
      in_file  = Ioss::Utils::local_filename(argv[i++], in_type, cwd);
      out_file = Ioss::Utils::local_filename(argv[i++], out_type, cwd);
    }
  }

  if (in_file.empty() || out_file.empty()) {
    show_usage(codename);
    return (EXIT_FAILURE);
  }

  DO_OUTPUT << "Input:    '" << in_file << "', Type: " << in_type << '\n';
  DO_OUTPUT << "Output:   '" << out_file << "', Type: " << out_type << '\n';
  DO_OUTPUT << '\n';

  create_sph(in_file, in_type, out_file, out_type, globals);

  DO_OUTPUT << "\n" << codename << " execution successful.\n";
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog)
  {
    DO_OUTPUT
        << "\nUSAGE: " << prog << " in_file out_file\n"
        << "...or: " << prog << " command_file\n"
        << "       version: " << version << "\n\n"
        << "\tConverts all HEX element blocks to SPHERE element blocks\n"
        << "\tand creates a nodeset for each element block containing the node at the center of "
           "the sphere.\n"
        << "\tignores all other element block types and deletes all existing nodesets.\n"
        << "Options:\n"
        << "\t--directory or -d {dir}  : specifies current working directory\n"
        << "\t--input     or -i {file} : read input and output filename data from file\n"
        << "\t--in_type {pamgen|generated|exodus} : set input type to the argument. Default "
           "exodus\n"
        << "\t--out_type {exodus} : set output type to the argument. Default exodus\n"
        << "\t--scale_factor : radius = (volume/scale_factor)^1/3; default = 8.0\n"
        << "\t\t which gives the half-side-length of a cube with that volume\n"
        << "\t\t use 4*pi/3 = 4.18879 for the radius of a sphere with that volume.\n";

    Ioss::NameList db_types;
    Ioss::IOFactory::describe(&db_types);
    DO_OUTPUT << "\nSupports database types:\n\t";
    for (Ioss::NameList::const_iterator IF = db_types.begin(); IF != db_types.end(); ++IF) {
      DO_OUTPUT << *IF << "  ";
    }
    DO_OUTPUT << "\n\n";
  }

  void create_sph(const std::string &inpfile, const std::string &input_type,
                  const std::string &outfile, const std::string &output_type, Globals &globals)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi =
        Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART, (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    if (!region.get_nodesets().empty()) {
      DO_OUTPUT << " *** All nodesets will be replaced by element block nodesets\n";
    }

    //========================================================================
    // OUTPUT ...
    //========================================================================
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(output_type, outfile, Ioss::WRITE_RESTART,
                                                    (MPI_Comm)MPI_COMM_WORLD);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");

    // Set the qa information...
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    // Get all elementblocks.
    // ... For all elementblocks of type 'HEX8'
    //    ... Calculate the element volume and the radius of a sphere
    //        which would have that same volume.
    //    ... Output a new element block of type 'SPHERE' with
    //        attributes volume and radius
    //    ... Output a nodeset with same id containing all nodes
    //        that are in the element block.
    // ...

    // Create the output mesh...
    output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    // Iterate through all element blocks...
    // The first time, just count the number of hex elements since
    // that will be the number of output element and nodes in the
    // sphere mesh.
    size_t                                      sph_node_count = 0;
    const Ioss::ElementBlockContainer &         ebs            = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator I              = ebs.begin();
    while (I != ebs.end()) {
      Ioss::ElementBlock *eb = *I;
      ++I;
      std::string type = eb->get_property("topology_type").get_string();
      if (type == Ioss::Hex8::name) {
        sph_node_count += eb->entity_count();

        // Add the element block...
        int                 num_elem = eb->entity_count();
        std::string         name     = eb->name();
        Ioss::ElementBlock *ebn =
            new Ioss::ElementBlock(output_region.get_database(), name, "sphere", num_elem);
        int id = eb->get_property("id").get_int();
        ebn->property_add(Ioss::Property("id", id));
        ebn->field_add(Ioss::Field("radius", Ioss::Field::REAL, "scalar", Ioss::Field::ATTRIBUTE,
                                   num_elem, 1));
        ebn->field_add(Ioss::Field("volume", Ioss::Field::REAL, "scalar", Ioss::Field::ATTRIBUTE,
                                   num_elem, 2));
        output_region.add(ebn);

        // Add a nodeset for each sph element block...
        name += "_nodes";
        auto ns = new Ioss::NodeSet(output_region.get_database(), name, num_elem);
        ns->property_add(Ioss::Property("id", id));
        output_region.add(ns);
      }
    }

    // Define a node block...
    std::string block_name        = "nodeblock_1";
    int         spatial_dimension = region.get_property("spatial_dimension").get_int();
    auto block = new Ioss::NodeBlock(output_region.get_database(), block_name, sph_node_count,
                                     spatial_dimension);
    block->property_add(Ioss::Property("id", 1));
    output_region.add(block);

    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    output_region.begin_mode(Ioss::STATE_MODEL);

    Ioss::NodeBlock *   nb = region.get_node_blocks()[0];
    std::vector<double> coordinates;
    nb->get_field_data("mesh_model_coordinates", coordinates);

    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    std::vector<int> node_ids(sph_node_count);
    for (size_t i = 0; i < sph_node_count; i++) {
      node_ids[i] = i + 1;
    }
    onb->put_field_data("ids", node_ids);

    std::vector<double> centroids(spatial_dimension * sph_node_count);

    I             = ebs.begin();
    size_t offset = 0;
    while (I != ebs.end()) {
      if ((*I)->get_property("topology_type").get_string() == Ioss::Hex8::name) {

        std::vector<double> volume;
        std::vector<double> radius;

        hex_volume(*I, coordinates, &centroids[3 * offset], volume, radius, globals.scale_factor);

        // Find corresponding output element block...
        Ioss::ElementBlock *output_eb = output_region.get_element_block((*I)->name());
        if (output_eb == nullptr) {
          std::cerr << "ERROR: Could not put find element block " << (*I)->name() << "\n";
          std::exit(EXIT_FAILURE);
        }

        Ioss::NodeSet *output_ns = output_region.get_nodeset((*I)->name() + "_nodes");
        if (output_ns == nullptr) {
          std::cerr << "ERROR: Could not put find node set " << (*I)->name() + "_nodes"
                    << "\n";
          std::exit(EXIT_FAILURE);
        }

        // Connectivity for the sphere element blocks is just their local element location
        size_t           num_elem = output_eb->entity_count();
        std::vector<int> connectivity(num_elem);
        for (size_t i = 0; i < num_elem; i++) {
          connectivity[i] = offset + 1 + i;
        }
        output_eb->put_field_data("ids", connectivity);
        output_eb->put_field_data("connectivity_raw", connectivity);
        output_eb->put_field_data("volume", volume);
        output_eb->put_field_data("radius", radius);

        output_ns->put_field_data("ids", connectivity);

        offset += output_eb->entity_count();
      }
      ++I;
    }

    onb->put_field_data("mesh_model_coordinates", centroids);

    output_region.end_mode(Ioss::STATE_MODEL);
  }

  void comp_grad12x(double *const grad_ptr, const double *const x, const double *const y,
                    const double *const z)
  {
    double R42 = (z[3] - z[1]);
    double R52 = (z[4] - z[1]);
    double R54 = (z[4] - z[3]);

    double R63 = (z[5] - z[2]);
    double R83 = (z[7] - z[2]);
    double R86 = (z[7] - z[5]);

    double R31 = (z[2] - z[0]);
    double R61 = (z[5] - z[0]);
    double R74 = (z[6] - z[3]);

    double R72 = (z[6] - z[1]);
    double R75 = (z[6] - z[4]);
    double R81 = (z[7] - z[0]);

    double t1 = (R63 + R54);
    double t2 = (R61 + R74);
    double t3 = (R72 + R81);

    double t4 = (R86 + R42);
    double t5 = (R83 + R52);
    double t6 = (R75 + R31);

    grad_ptr[0] =
        (y[1] * t1) - (y[2] * R42) - (y[3] * t5) + (y[4] * t4) + (y[5] * R52) - (y[7] * R54);
    grad_ptr[1] =
        (y[2] * t2) + (y[3] * R31) - (y[0] * t1) - (y[5] * t6) + (y[6] * R63) - (y[4] * R61);
    grad_ptr[2] =
        (y[3] * t3) + (y[0] * R42) - (y[1] * t2) - (y[6] * t4) + (y[7] * R74) - (y[5] * R72);
    grad_ptr[3] =
        (y[0] * t5) - (y[1] * R31) - (y[2] * t3) + (y[7] * t6) + (y[4] * R81) - (y[6] * R83);
    grad_ptr[4] =
        (y[5] * t3) + (y[6] * R86) - (y[7] * t2) - (y[0] * t4) - (y[3] * R81) + (y[1] * R61);
    grad_ptr[5] =
        (y[6] * t5) - (y[4] * t3) - (y[7] * R75) + (y[1] * t6) - (y[0] * R52) + (y[2] * R72);
    grad_ptr[6] =
        (y[7] * t1) - (y[5] * t5) - (y[4] * R86) + (y[2] * t4) - (y[1] * R63) + (y[3] * R83);
    grad_ptr[7] =
        (y[4] * t2) - (y[6] * t1) + (y[5] * R75) - (y[3] * t6) - (y[2] * R74) + (y[0] * R54);

    R42 = (x[3] - x[1]);
    R52 = (x[4] - x[1]);
    R54 = (x[4] - x[3]);

    R63 = (x[5] - x[2]);
    R83 = (x[7] - x[2]);
    R86 = (x[7] - x[5]);

    R31 = (x[2] - x[0]);
    R61 = (x[5] - x[0]);
    R74 = (x[6] - x[3]);

    R72 = (x[6] - x[1]);
    R75 = (x[6] - x[4]);
    R81 = (x[7] - x[0]);

    t1 = (R63 + R54);
    t2 = (R61 + R74);
    t3 = (R72 + R81);

    t4 = (R86 + R42);
    t5 = (R83 + R52);
    t6 = (R75 + R31);

    grad_ptr[8] =
        (z[1] * t1) - (z[2] * R42) - (z[3] * t5) + (z[4] * t4) + (z[5] * R52) - (z[7] * R54);
    grad_ptr[9] =
        (z[2] * t2) + (z[3] * R31) - (z[0] * t1) - (z[5] * t6) + (z[6] * R63) - (z[4] * R61);
    grad_ptr[10] =
        (z[3] * t3) + (z[0] * R42) - (z[1] * t2) - (z[6] * t4) + (z[7] * R74) - (z[5] * R72);
    grad_ptr[11] =
        (z[0] * t5) - (z[1] * R31) - (z[2] * t3) + (z[7] * t6) + (z[4] * R81) - (z[6] * R83);
    grad_ptr[12] =
        (z[5] * t3) + (z[6] * R86) - (z[7] * t2) - (z[0] * t4) - (z[3] * R81) + (z[1] * R61);
    grad_ptr[13] =
        (z[6] * t5) - (z[4] * t3) - (z[7] * R75) + (z[1] * t6) - (z[0] * R52) + (z[2] * R72);
    grad_ptr[14] =
        (z[7] * t1) - (z[5] * t5) - (z[4] * R86) + (z[2] * t4) - (z[1] * R63) + (z[3] * R83);
    grad_ptr[15] =
        (z[4] * t2) - (z[6] * t1) + (z[5] * R75) - (z[3] * t6) - (z[2] * R74) + (z[0] * R54);

    R42 = (y[3] - y[1]);
    R52 = (y[4] - y[1]);
    R54 = (y[4] - y[3]);

    R63 = (y[5] - y[2]);
    R83 = (y[7] - y[2]);
    R86 = (y[7] - y[5]);

    R31 = (y[2] - y[0]);
    R61 = (y[5] - y[0]);
    R74 = (y[6] - y[3]);

    R72 = (y[6] - y[1]);
    R75 = (y[6] - y[4]);
    R81 = (y[7] - y[0]);

    t1 = (R63 + R54);
    t2 = (R61 + R74);
    t3 = (R72 + R81);

    t4 = (R86 + R42);
    t5 = (R83 + R52);
    t6 = (R75 + R31);

    grad_ptr[16] =
        (x[1] * t1) - (x[2] * R42) - (x[3] * t5) + (x[4] * t4) + (x[5] * R52) - (x[7] * R54);
    grad_ptr[17] =
        (x[2] * t2) + (x[3] * R31) - (x[0] * t1) - (x[5] * t6) + (x[6] * R63) - (x[4] * R61);
    grad_ptr[18] =
        (x[3] * t3) + (x[0] * R42) - (x[1] * t2) - (x[6] * t4) + (x[7] * R74) - (x[5] * R72);
    grad_ptr[19] =
        (x[0] * t5) - (x[1] * R31) - (x[2] * t3) + (x[7] * t6) + (x[4] * R81) - (x[6] * R83);
    grad_ptr[20] =
        (x[5] * t3) + (x[6] * R86) - (x[7] * t2) - (x[0] * t4) - (x[3] * R81) + (x[1] * R61);
    grad_ptr[21] =
        (x[6] * t5) - (x[4] * t3) - (x[7] * R75) + (x[1] * t6) - (x[0] * R52) + (x[2] * R72);
    grad_ptr[22] =
        (x[7] * t1) - (x[5] * t5) - (x[4] * R86) + (x[2] * t4) - (x[1] * R63) + (x[3] * R83);
    grad_ptr[23] =
        (x[4] * t2) - (x[6] * t1) + (x[5] * R75) - (x[3] * t6) - (x[2] * R74) + (x[0] * R54);
  }

  inline double dot8(const double *const x1, const double *const x2)
  {
    double d1 = x1[0] * x2[0] + x1[1] * x2[1];
    double d2 = x1[2] * x2[2] + x1[3] * x2[3];
    double d4 = x1[4] * x2[4] + x1[5] * x2[5];
    double d5 = x1[6] * x2[6] + x1[7] * x2[7];
    return d1 + d2 + d4 + d5;
  }

  void hex_volume(Ioss::ElementBlock *block, const std::vector<double> &coordinates,
                  double *centroids, std::vector<double> &volume, std::vector<double> &radius,
                  double scale_factor)
  {
    const double     one12th = 1.0 / 12.0;
    size_t           nelem   = block->entity_count();
    std::vector<int> connectivity;
    block->get_field_data("connectivity_raw", connectivity);

    volume.resize(nelem);
    radius.resize(nelem);

    double gradop12x[24];
    double x[8];
    double y[8];
    double z[8];

    for (size_t ielem = 0; ielem < nelem; ++ielem) {
      for (size_t j = 0; j < 8; j++) {
        size_t node = connectivity[ielem * 8 + j] - 1;
        x[j]        = coordinates[node * 3 + 0];
        y[j]        = coordinates[node * 3 + 1];
        z[j]        = coordinates[node * 3 + 2];
      }

      centroids[ielem * 3 + 0] = .125 * (x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7]);
      centroids[ielem * 3 + 1] = .125 * (y[0] + y[1] + y[2] + y[3] + y[4] + y[5] + y[6] + y[7]);
      centroids[ielem * 3 + 2] = .125 * (z[0] + z[1] + z[2] + z[3] + z[4] + z[5] + z[6] + z[7]);

      comp_grad12x(&gradop12x[0], x, y, z);
      const double volume12x = dot8(x, &gradop12x[0]);
      volume[ielem]          = volume12x * one12th;
      radius[ielem]          = std::pow(volume[ielem] / scale_factor, (1. / 3.));
    }
  }
} // namespace
