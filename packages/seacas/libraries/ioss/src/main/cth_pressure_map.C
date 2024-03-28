// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ionit_Initializer.h"
#include "Ioss_Utils.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_EntityType.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_State.h"
#include "Ioss_VariableType.h"
#include "vector3d.h"

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    enum FinalPressureType { FINAL, INITIAL, ZERO, OFFSET };
    bool              debug{};
    bool              do_normals{};
    bool              reverse_normals{};
    bool              add_sset{};
    bool              convert_gage{};
    FinalPressureType final_pressure{FINAL};
    double            delta_time{};
    double            maximum_time{};
    double            minimum_time{};
    double            offset_time{};
    double            offset_pressure{};
    std::string       working_directory{};
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
  void add_sideset_mesh_fields(Ioss::Region &ss_region, Ioss::Region &region,
                               Ioss::Region &output_region, Globals &globals);
  void add_sideset_transient_fields(Ioss::Region &ss_region, Ioss::Region &region,
                                    Ioss::Region &output_region, Globals &globals);

  void transfer_sideset_field_data(Ioss::Region &ss_region, Ioss::Region &region,
                                   Ioss::Region &output_region, Globals &globals);

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix = "");

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const std::string &prefix = "",
                           bool transfer_connectivity = true);

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge);

  void output_normals(Ioss::Region &region, Ioss::Region &output_region, bool reverse_normals);
  void calculate_normals(std::vector<double> &node_normal, std::vector<double> &face_normal,
                         int num_elem, int num_node_per_elem, std::vector<double> &coord,
                         std::vector<int> &conn, bool reverse_normals);
  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &field_name);

  void file_copy(const std::string &inpfile, const std::string &input_type,
                 const std::string &outfile, const std::string &output_type,
                 const std::string &ss_file, Globals &globals);
} // namespace
// ========================================================================

namespace {
  std::string codename;
  std::string version = "$Revision$";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  std::string in_type  = "exodusII";
  std::string out_type = "exodusII";
  std::string ss_type  = "exodusII";

  Globals globals;

  globals.do_normals      = false;
  globals.reverse_normals = false;
  globals.add_sset        = true; // This is the CTH Pressure map executable
  globals.convert_gage    = false;
  globals.final_pressure  = Globals::FINAL;
  globals.delta_time      = 0.0;
  globals.maximum_time    = 0.0;
  globals.minimum_time    = 0.0;
  globals.offset_time     = 0.0;
  globals.offset_pressure = 0.0;

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  Ioss::Init::Initializer io;

  globals.debug = false;

  std::string input_file;
  // Skip past any options...
  int i = 1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("-directory", argv[i]) == 0 || std::strcmp("--directory", argv[i]) == 0 ||
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
      globals.add_sset   = false;
      i++;
    }
    else if (std::strcmp("--Convert_Gage", argv[i]) == 0) {
      globals.convert_gage = true;
      i++;
    }
    else if (std::strcmp("--Final_Pressure", argv[i]) == 0) {
      i++;
      if (std::strcmp("ZERO", argv[i]) == 0) {
        globals.final_pressure = Globals::ZERO;
      }
      else if (std::strcmp("INITIAL", argv[i]) == 0) {
        globals.final_pressure = Globals::INITIAL;
      }
      else if (std::strcmp("FINAL", argv[i]) == 0) {
        globals.final_pressure = Globals::FINAL;
      }
      else if (std::strcmp("OFFSET", argv[i]) == 0) {
        globals.final_pressure = Globals::OFFSET;
      }
      i++;
    }
    else if (std::strcmp("--Final_Time_Delta", argv[i]) == 0) {
      i++;
      globals.delta_time = std::strtod(argv[i++], nullptr);
    }
    else if (std::strcmp("--Maximum_Time", argv[i]) == 0) {
      i++;
      globals.maximum_time = std::strtod(argv[i++], nullptr);
    }
    else if (std::strcmp("--Minimum_Time", argv[i]) == 0) {
      i++;
      globals.minimum_time = std::strtod(argv[i++], nullptr);
    }
    else if (std::strcmp("--Offset_Time", argv[i]) == 0) {
      i++;
      globals.offset_time = std::strtod(argv[i++], nullptr);
    }
    else if (std::strcmp("--Offset_Pressure", argv[i]) == 0) {
      i++;
      globals.offset_pressure = std::strtod(argv[i++], nullptr);
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
      in_file   = Ioss::Utils::local_filename(argv[i++], in_type, cwd);
      sset_file = Ioss::Utils::local_filename(argv[i++], ss_type, cwd);
      out_file  = Ioss::Utils::local_filename(argv[i++], out_type, cwd);
    }
    else if (argc - i == 2) {
      in_file  = Ioss::Utils::local_filename(argv[i++], in_type, cwd);
      out_file = Ioss::Utils::local_filename(argv[i++], out_type, cwd);
    }
    else if (argc - i == 1) {
      input_file = Ioss::Utils::local_filename(argv[i++], "text", cwd);
    }
  }

  if (!input_file.empty()) {
    std::ifstream input(input_file);
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
    return (EXIT_FAILURE);
  }

  std::cerr << "Input:    '" << in_file << "', Type: " << in_type << '\n';
  if (globals.add_sset) {
    std::cerr << "Pressure: '" << sset_file << "', Type: " << in_type << '\n';
  }
  std::cerr << "Output:   '" << out_file << "', Type: " << out_type << '\n';
  std::cerr << '\n';

  if (globals.add_sset || globals.do_normals) {
    globals.debug = false;
  }

  if (!(globals.add_sset || globals.do_normals)) {
    std::cerr << "\n" << codename << "ERROR: Either add sideset or do normals must be selected.\n";
    return EXIT_FAILURE;
  }

  file_copy(in_file, in_type, out_file, out_type, sset_file, globals);

  std::cerr << "\n" << codename << " execution successful.\n";
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog, bool add_sset)
  {
    if (add_sset) {
      std::cerr << "\nUSAGE: " << prog << " in_file sset_file out_file\n";
    }
    else {
      std::cerr << "\nUSAGE: " << prog << " in_file out_file\n";
    }
    std::cerr << "...or: " << prog << " command_file\n";
    std::cerr << "       version: " << version << "\n";
    Ioss::NameList db_types = Ioss::IOFactory::describe();
    std::cerr << "\nSupports database types:\n\t";
    for (const auto &db : db_types) {
      std::cerr << db << "  ";
    }
    std::cerr << "\n\n";
  }

  void file_copy(const std::string &inpfile, const std::string &input_type,
                 const std::string &outfile, const std::string &output_type,
                 const std::string &ss_file, Globals &globals)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
                                                    Ioss::ParallelUtils::comm_world());
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
                                                    Ioss::ParallelUtils::comm_world());
    if (dbo == nullptr || !dbo->ok()) {
      std::cerr << "ERROR: Could not create output database '" << outfile << "' of type '"
                << output_type << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");

    //========================================================================
    // Optional sideset file...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbs       = nullptr;
    Ioss::Region     *ss_region = nullptr;
    if (globals.add_sset) {
      // Open the file containing the data which will be the new sideset...
      dbs = Ioss::IOFactory::create(input_type, ss_file, Ioss::READ_RESTART,
                                    Ioss::ParallelUtils::comm_world());
      if (dbs == nullptr || !dbs->ok()) {
        std::cerr << "ERROR: Could not open database '" << ss_file << "' of type '" << input_type
                  << "'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      ss_region = new Ioss::Region(dbs, "sideset_region");
    }

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
    transfer_nodesets(region, output_region, globals.debug);
    if (!globals.do_normals) {
      transfer_sidesets(region, output_region, globals.debug);
    }
    transfer_commsets(region, output_region, globals.debug);

    // If adding the sideset record from an external file, do it here...
    if (globals.add_sset) {
      add_sideset(*ss_region, region, output_region, globals);
    }

    if (globals.debug) {
      std::cerr << "END STATE_DEFINE_MODEL... " << '\n';
    }
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    if (globals.debug) {
      std::cerr << "TRANSFERRING MESH FIELD DATA ... " << '\n';
      // Model defined, now fill in the model data...
    }
    output_region.begin_mode(Ioss::STATE_MODEL);

    // Transfer MESH field_data from input to output...
    {
      Ioss::NodeBlock *nb = region.get_node_blocks()[0];
      assert(nb != nullptr);

      if (nb != nullptr) {
        const std::string &name = nb->name();
        if (globals.debug) {
          std::cerr << name << ", ";

          // Find the corresponding output element_block...
        }
        Ioss::NodeBlock *onb = output_region.get_node_block(name);
        assert(onb != nullptr);
        transfer_field_data(nb, onb, Ioss::Field::MESH);
        transfer_field_data(nb, onb, Ioss::Field::ATTRIBUTE);
      }
      if (globals.debug) {
        std::cerr << '\n';
      }
    }

    // Now do the same for element blocks...
    {
      const auto &ebs = region.get_element_blocks();
      auto        I   = ebs.cbegin();

      while (I != ebs.end()) {
        const std::string &name = (*I)->name();
        if (globals.debug) {
          std::cerr << name << ", ";

          // Find the corresponding output element_block...
        }
        Ioss::ElementBlock *oeb = output_region.get_element_block(name);
        assert(oeb != nullptr);

        if (oeb != nullptr) {
          transfer_field_data(*I, oeb, Ioss::Field::MESH);
          transfer_field_data(*I, oeb, Ioss::Field::ATTRIBUTE);
        }
        ++I;
      }
      if (globals.debug) {
        std::cerr << '\n';
      }
    }

    // Node Sets
    {
      const auto &nss = region.get_nodesets();
      auto        I   = nss.cbegin();
      while (I != nss.end()) {
        const std::string &name = (*I)->name();
        if (globals.debug) {
          std::cerr << name << ", ";

          // Find matching output nodeset
        }
        Ioss::NodeSet *ons = output_region.get_nodeset(name);
        if (ons != nullptr) {
          transfer_field_data(*I, ons, Ioss::Field::MESH);
          transfer_field_data(*I, ons, Ioss::Field::ATTRIBUTE);
        }
        ++I;
      }
      if (globals.debug) {
        std::cerr << '\n';
      }
    }

    // Comm Sets
    {
      const auto &css = region.get_commsets();
      auto        I   = css.cbegin();
      while (I != css.end()) {
        const std::string &name = (*I)->name();
        if (globals.debug) {
          std::cerr << name << ", ";

          // Find matching output commset
        }
        Ioss::CommSet *ocs = output_region.get_commset(name);
        if (ocs != nullptr) {
          transfer_field_data(*I, ocs, Ioss::Field::MESH);
          transfer_field_data(*I, ocs, Ioss::Field::ATTRIBUTE);
          transfer_field_data(*I, ocs, Ioss::Field::COMMUNICATION);
        }
        ++I;
      }
      if (globals.debug) {
        std::cerr << '\n';
      }
    }

    // Side Sets
    if (!globals.do_normals) {
      const auto &fss = region.get_sidesets();
      auto        I   = fss.cbegin();
      while (I != fss.end()) {
        const std::string &name = (*I)->name();
        if (globals.debug) {
          std::cerr << name << ", ";

          // Find matching output sideset
        }
        Ioss::SideSet *ofs = output_region.get_sideset(name);

        if (ofs != nullptr) {
          transfer_field_data(*I, ofs, Ioss::Field::MESH);
          transfer_field_data(*I, ofs, Ioss::Field::ATTRIBUTE);

          const auto &fbs = (*I)->get_side_blocks();
          auto        J   = fbs.cbegin();
          while (J != fbs.end()) {

            // Find matching output sideblock
            const std::string &fbname = (*J)->name();
            if (globals.debug) {
              std::cerr << fbname << ", ";
            }
            Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

            if (ofb != nullptr) {
              transfer_field_data(*J, ofb, Ioss::Field::MESH, "", false);
              transfer_field_data(*J, ofb, Ioss::Field::ATTRIBUTE, "", false);
            }
            ++J;
          }
        }
        ++I;
      }
      if (globals.debug) {
        std::cerr << '\n';
      }
    }

    if (globals.add_sset) {
      add_sideset_mesh_fields(*ss_region, region, output_region, globals);
    }

    if (globals.debug) {
      std::cerr << "END STATE_MODEL... " << '\n';
    }
    output_region.end_mode(Ioss::STATE_MODEL);

    if (globals.do_normals) {
      output_normals(region, output_region, globals.reverse_normals);
    }
    else if (globals.add_sset) {
      if (ss_region->property_exists("state_count") &&
          ss_region->get_property("state_count").get_int() > 0) {
        output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        add_sideset_transient_fields(*ss_region, region, output_region, globals);
        output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
      }
    }
    else {
      std::cerr << "Internal Error\n";
      std::abort();
    }

    if (globals.debug) {
      std::cerr << "TRANSFERRING TRANSIENT FIELDS ... " << '\n';
    }
    output_region.begin_mode(Ioss::STATE_TRANSIENT);
    // Get the timesteps from the input database.  Step through them
    // and transfer fields to output database...

    if (globals.do_normals) {
      // Do nothing, normals were already output above...
    }
    else if (globals.add_sset) {
      transfer_sideset_field_data(*ss_region, region, output_region, globals);
    }
    else {
      std::cerr << "Internal Error\n";
      std::abort();
    }
    if (globals.debug) {
      std::cerr << "END STATE_TRANSIENT... " << '\n';
    }
    output_region.end_mode(Ioss::STATE_TRANSIENT);
    delete ss_region;
  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nbs = region.get_node_blocks();
    auto        i   = nbs.cbegin();
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

      auto nb = new Ioss::NodeBlock(output_region.get_database(), name, num_nodes, degree);
      output_region.add(nb);

      transfer_properties(*i, nb);
      transfer_fields(*i, nb, Ioss::Field::MESH);
      transfer_fields(*i, nb, Ioss::Field::ATTRIBUTE);
      ++i;
    }
    if (debug) {
      std::cerr << '\n';
    }
  }

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &ebs            = region.get_element_blocks();
    auto        i              = ebs.cbegin();
    int         total_elements = 0;
    while (i != ebs.end()) {
      const std::string &name = (*i)->name();
      if (debug) {
        std::cerr << name << ", ";
      }
      std::string type     = (*i)->topology()->name();
      int         num_elem = (*i)->entity_count();
      total_elements += num_elem;

      auto eb = new Ioss::ElementBlock(output_region.get_database(), name, type, num_elem);
      output_region.add(eb);
      transfer_properties(*i, eb);
      transfer_fields(*i, eb, Ioss::Field::MESH);
      transfer_fields(*i, eb, Ioss::Field::ATTRIBUTE);
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

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &fss         = region.get_sidesets();
    auto        i           = fss.cbegin();
    int         total_sides = 0;
    while (i != fss.end()) {
      const std::string &name = (*i)->name();
      if (debug) {
        std::cerr << name << ", ";
      }
      auto surf = new Ioss::SideSet(output_region.get_database(), name);

      const auto &fbs = (*i)->get_side_blocks();
      auto        j   = fbs.cbegin();
      while (j != fbs.end()) {
        const std::string &fbname = (*j)->name();
        if (debug) {
          std::cerr << fbname << ", ";
        }
        std::string fbtype   = (*j)->topology()->name();
        std::string partype  = (*j)->parent_element_topology()->name();
        int         num_side = (*j)->entity_count();
        total_sides += num_side;

        auto block =
            new Ioss::SideBlock(output_region.get_database(), fbname, fbtype, partype, num_side);
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
    }
    else {
      std::cerr << '\n';
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nss         = region.get_nodesets();
    auto        i           = nss.cbegin();
    int         total_nodes = 0;
    while (i != nss.end()) {
      const std::string &name = (*i)->name();
      if (debug) {
        std::cerr << name << ", ";
      }
      int count = (*i)->entity_count();
      total_nodes += count;
      auto ns = new Ioss::NodeSet(output_region.get_database(), name, count);
      output_region.add(ns);
      transfer_properties(*i, ns);
      transfer_fields(*i, ns, Ioss::Field::MESH);
      transfer_fields(*i, ns, Ioss::Field::ATTRIBUTE);
      ++i;
    }
    if (!debug) {
      std::cerr << " Number of nodal point sets           =" << std::setw(9) << nss.size() << "\n";
      std::cerr << "     Length of node list              =" << std::setw(9) << total_nodes << "\n";
    }
    else {
      std::cerr << '\n';
    }
  }

  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &css = region.get_commsets();
    auto        i   = css.cbegin();
    while (i != css.end()) {
      const std::string &name = (*i)->name();
      if (debug) {
        std::cerr << name << ", ";
      }
      std::string type  = (*i)->get_property("entity_type").get_string();
      int         count = (*i)->entity_count();
      auto        cs    = new Ioss::CommSet(output_region.get_database(), name, type, count);
      output_region.add(cs);
      transfer_properties(*i, cs);
      transfer_fields(*i, cs, Ioss::Field::MESH);
      transfer_fields(*i, cs, Ioss::Field::ATTRIBUTE);
      transfer_fields(*i, cs, Ioss::Field::COMMUNICATION);
      ++i;
    }
    if (debug) {
      std::cerr << '\n';
    }
  }

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix)
  {
    // Check for transient fields...
    Ioss::NameList fields = ige->field_describe(role);

    // Iterate through results fields and transfer to output
    // database...  If a prefix is specified, only transfer fields
    // whose names begin with the prefix
    for (const auto &field_name : fields) {
      if (field_name != "ids" && !oge->field_exists(field_name) &&
          Ioss::Utils::substr_equal(prefix, field_name)) {
        // If the field does not already exist, add it to the output node block
        Ioss::Field field = ige->get_field(field_name);
        oge->field_add(std::move(field));
      }
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const std::string &prefix,
                           bool transfer_connectivity)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields = ige->field_describe(role);

    // Complication here is that if the 'role' is 'Ioss::Field::MESH',
    // then the 'ids' field must be transferred first...
    if (role == Ioss::Field::MESH) {
      for (const auto &field_name : state_fields) {
        assert(oge->field_exists(field_name));
        if (field_name == "ids") {
          transfer_field_data_internal(ige, oge, field_name);
          break;
        }
      }
    }

    for (const auto &field_name : state_fields) {
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (!transfer_connectivity && field_name == "connectivity") {
        continue;
      }

      if (field_name != "ids" && Ioss::Utils::substr_equal(prefix, field_name)) {
        assert(oge->field_exists(field_name));
        transfer_field_data_internal(ige, oge, field_name);
      }
    }
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &field_name)
  {

    size_t isize = ige->get_field(field_name).get_size();
    assert(isize == oge->get_field(field_name).get_size());

    data.resize(isize);

    if (field_name == "mesh_model_coordinates_x") {
      return;
    }
    if (field_name == "mesh_model_coordinates_y") {
      return;
    }
    if (field_name == "mesh_model_coordinates_z") {
      return;
    }
    if (field_name == "connectivity_raw") {
      return;
    }
    if (field_name == "element_side_raw") {
      return;
    }
    if (field_name == "ids_raw") {
      return;
    }
    if (field_name == "node_connectivity_status") {
      return;
    }
    if (field_name == "entity_processor_raw") {
      return;
    }
    if (field_name == "owning_processor") {
      return;
    }
    if (ige->type() == Ioss::SIDEBLOCK && field_name == "ids") {
      return;
    }
    ige->get_field_data(field_name, Data(data), isize);
    oge->put_field_data(field_name, Data(data), isize);
  }

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge)
  {
    Ioss::NameList names = ige->property_describe();

    // Iterate through properties and transfer to output database...
    Ioss::NameList::const_iterator I;
    for (I = names.begin(); I != names.end(); ++I) {
      if (!oge->property_exists(*I)) {
        oge->property_add(ige->get_property(*I));
      }
    }
  }

  void add_sideset(Ioss::Region &ss_region, Ioss::Region & /* region */,
                   Ioss::Region &output_region, Globals &globals)
  {
    auto *pressures = new Ioss::SideSet(output_region.get_database(), "cth_pressures");
    // Each element block in the sset file will be a surside in the mesh file...
    {
      const auto &ebs = ss_region.get_element_blocks();
      auto        i   = ebs.begin();
      while (i != ebs.end()) {
        std::string name = (*i)->name();
        name             = "ss" + name;

        if (globals.debug) {
          std::cerr << name << ", ";
        }
        int         num_elem = (*i)->entity_count();
        std::string type     = (*i)->topology()->name();

        // Should be able to get this from the input mesh element blocks...
        std::string partype = "unknown";
        auto fb = new Ioss::SideBlock(output_region.get_database(), name, type, partype, num_elem);
        pressures->add(fb);
        ++i;
      }
      output_region.add(pressures);
    }
  }

  void add_sideset_mesh_fields(Ioss::Region &ss_region, Ioss::Region & /* region */,
                               Ioss::Region &output_region, Globals &globals)
  {
    // Each element block in the sset file will be a surface in the mesh file...
    {
      const auto &ebs = ss_region.get_element_blocks();
      auto        i   = ebs.begin();
      while (i != ebs.end()) {
        std::string name    = (*i)->name();
        name                = "ss" + name;
        Ioss::SideBlock *fb = output_region.get_sideblock(name);
        if (fb == nullptr) {
          std::cerr << "INTERNAL ERROR: Could not find sideblock named '" << name << "'\n";
          std::exit(EXIT_FAILURE);
        }

        if (globals.debug) {
          std::cerr << name << ", ";

          // The "skin" field should exist on each element block.  It
          // contains the map back to the element/local face in the
          // original (region/output_region) mesh.  It is in the exact
          // same format as a sidesets "element_side" field...
          //
          // See if the field exists...
        }
        if ((*i)->field_exists("skin")) {
          Ioss::Field skin = (*i)->get_field("skin");

          int isize = skin.get_size();
          int osize = fb->get_field("element_side").get_size();
          assert(isize == osize);

          data.resize(isize);

          (*i)->get_field_data("skin", Data(data), isize);
          fb->put_field_data("element_side", Data(data), osize);
          ++i;
        }
        else {
          std::cerr << "ERROR: Field 'skin' does not exist on element block '" << name << "'.\n";
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }

  void transfer_sideset_field_data(Ioss::Region &ss_region, Ioss::Region & /* region */,
                                   Ioss::Region &output_region, Globals &globals)
  {
    int step_count = ss_region.get_property("state_count").get_int();
    int istep      = 1; // Need this to exist after end of step_count loop.

    int initial_state = 1;
    // Find first state >= globals.minimum_time;
    for (; initial_state <= step_count; initial_state++) {
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
    std::map<std::string, double *> time_zero_field_data;
    if (globals.convert_gage) {
      ss_region.begin_state(1);

      const auto &ebs = ss_region.get_element_blocks();
      auto        i   = ebs.begin();
      while (i != ebs.end()) {

        // The gage pressure conversion is currently only applied to the field "cth_pressure"
        const std::string &name = (*i)->name();
        if ((*i)->field_exists(cth_pressure)) {
          int   isize                = (*i)->get_field(cth_pressure).get_size();
          void *zdata                = new char[isize];
          time_zero_field_data[name] = reinterpret_cast<double *>(zdata);
          (*i)->get_field_data(cth_pressure, zdata, isize);
        }
        else {
          time_zero_field_data[name] = (double *)nullptr;
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
      const auto &ebs = ss_region.get_element_blocks();
      auto        i   = ebs.begin();
      while (i != ebs.end()) {
        const std::string &eb_name = (*i)->name();
        std::string        name    = "ss" + eb_name;
        Ioss::SideBlock   *fb      = output_region.get_sideblock(name);
        if (fb == nullptr) {
          std::cerr << "INTERNAL ERROR: Could not find sideblock named '" << name << "'\n";
          std::exit(EXIT_FAILURE);
        }

        auto state_fields = (*i)->field_describe(Ioss::Field::TRANSIENT);
        for (const auto &field_name : state_fields) {
          // NOTE: Only dealing with the "cth_" fields here.
          // If there are other fields, we probably have an invalid
          // output database...
          if (Ioss::Utils::substr_equal("cth_", field_name)) {
            int isize = (*i)->get_field(field_name).get_size();
            int count = (*i)->get_field(field_name).raw_count();
            data.resize(isize);
            auto *rdata = reinterpret_cast<double *>(Data(data));
            for (int ii = 0; ii < count; ii++) {
              rdata[ii] = globals.offset_pressure;
            }
            fb->put_field_data(field_name, Data(data), isize);
            break;
          }
        }
        ++i;
      }
      output_region.end_state(ostep);
    }

    for (istep = initial_state; istep <= step_count; istep++) {
      double time = ss_region.get_state_time(istep);
      if (time < globals.minimum_time) {
        continue;
      }

      if (globals.maximum_time > 0.0 && time > globals.maximum_time) {
        break;
      }

      int ostep = output_region.add_state(time - globals.minimum_time + globals.offset_time);
      show_step(istep, time);

      output_region.begin_state(ostep);
      ss_region.begin_state(istep);
      const auto &ebs = ss_region.get_element_blocks();
      auto        i   = ebs.begin();
      while (i != ebs.end()) {
        const std::string &eb_name = (*i)->name();
        std::string        name    = "ss" + eb_name;
        auto              *fb      = output_region.get_sideblock(name);
        if (fb == nullptr) {
          std::cerr << "INTERNAL ERROR: Could not find side block named '" << name << "'\n";
          std::exit(EXIT_FAILURE);
        }

        auto state_fields = (*i)->field_describe(Ioss::Field::TRANSIENT);
        for (const auto &field_name : state_fields) {
          if (globals.convert_gage && field_name == cth_pressure) {
            // Subtract the time zero pressures (stored in
            // time_zero_field_data) from each time step
            double *zdata = time_zero_field_data[eb_name];
            assert(zdata != nullptr);

            int isize = (*i)->get_field(field_name).get_size();
            int count = (*i)->get_field(field_name).raw_count();
            data.resize(isize);
            auto *rdata = reinterpret_cast<double *>(Data(data));

            (*i)->get_field_data(field_name, Data(data), isize);
            for (int ii = 0; ii < count; ii++) {
              rdata[ii] -= zdata[ii];
            }

            if (globals.offset_pressure != 0.0) {
              for (int ii = 0; ii < count; ii++) {
                rdata[ii] += globals.offset_pressure;
              }
            }
            fb->put_field_data(field_name, Data(data), isize);
          }
          else if (globals.offset_pressure != 0.0 && field_name == cth_pressure) {
            int isize = (*i)->get_field(field_name).get_size();
            int count = (*i)->get_field(field_name).raw_count();
            data.resize(isize);
            auto *rdata = reinterpret_cast<double *>(Data(data));

            (*i)->get_field_data(field_name, Data(data), isize);
            for (int ii = 0; ii < count; ii++) {
              rdata[ii] += globals.offset_pressure;
            }
            fb->put_field_data(field_name, Data(data), isize);
          }
          else if (Ioss::Utils::substr_equal("cth_", field_name)) {
            assert(fb->field_exists(field_name));
            transfer_field_data_internal(*i, fb, field_name);
          }
        }
        ++i;
      }
      ss_region.end_state(istep);
      output_region.end_state(ostep);
    }

    step_count = istep - 1;

    // See if special handling of the last step is specified...
    if (globals.final_pressure == Globals::INITIAL || globals.final_pressure == Globals::ZERO ||
        globals.final_pressure == Globals::OFFSET) {
      double time = ss_region.get_state_time(step_count);
      double delta_time;
      if (globals.delta_time > 0.0) {
        delta_time = globals.delta_time;
      }
      else {
        delta_time = time - ss_region.get_state_time(step_count - 1);
      }
      time += delta_time;
      int ostep = output_region.add_state(time - globals.minimum_time + globals.offset_time);

      output_region.begin_state(ostep);
      ss_region.begin_state(initial_state);

      // Repeat the data from the first step at the last step.  This is to bring the model back to
      // an
      // "equilibrium" state in case the CTH analysis was not run out to an equilibrium
      // state. If ZERO was specified, then it simply zeros out the pressure field at the last step.
      const auto &ebs = ss_region.get_element_blocks();
      auto        i   = ebs.begin();
      while (i != ebs.end()) {
        const std::string &eb_name = (*i)->name();
        std::string        name    = "ss" + eb_name;
        Ioss::SideBlock   *fb      = output_region.get_sideblock(name);
        if (fb == nullptr) {
          std::ostringstream msg;
          msg << " INTERNAL_ERROR: Could not find sideblock '" << name << "'\n";
          throw std::runtime_error(msg.str());
        }

        auto state_fields = (*i)->field_describe(Ioss::Field::TRANSIENT);
        for (const auto &field_name : state_fields) {
          if (Ioss::Utils::substr_equal("cth_", field_name)) {
            if (field_name == cth_pressure &&
                (globals.final_pressure == Globals::ZERO ||
                 globals.final_pressure == Globals::OFFSET || globals.convert_gage)) {

              // If convert_gage is true and we are outputting the
              // cth_pressure variable, then for both the INITIAL and
              // ZERO cases, we output a zero-filled field...

              double value = 0.0;
              if (globals.final_pressure == Globals::OFFSET) {
                value = globals.offset_pressure;
              }

              int isize = (*i)->get_field(field_name).get_size();
              int count = (*i)->get_field(field_name).raw_count();
              data.resize(isize);
              auto *rdata = reinterpret_cast<double *>(Data(data));
              for (int ii = 0; ii < count; ii++) {
                rdata[ii] = value;
              }
              fb->put_field_data(field_name, Data(data), isize);
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
      auto i  = time_zero_field_data.begin();
      auto ie = time_zero_field_data.end();
      while (i != ie) {
        delete[] reinterpret_cast<char *>((*i).second);
        ++i;
      }
    }
  }

  void add_sideset_transient_fields(Ioss::Region &ss_region, Ioss::Region & /* region */,
                                    Ioss::Region &output_region, Globals &globals)
  {
    // Each element block in the sset file will be a surface in the mesh file...
    {
      const auto &ebs = ss_region.get_element_blocks();
      for (auto *eb : ebs) {
        std::string name = eb->name();
        name             = "ss" + name;
        auto *fb         = output_region.get_sideblock(name);
        assert(fb != nullptr);

        if (globals.debug) {
          std::cerr << name << ", ";

          // Each element variable in the sset file which begins with
          // "cth_" will be a sideset variable in the outptut file...
        }
        transfer_fields(eb, fb, Ioss::Field::TRANSIENT, "cth_");
      }
    }
  }

  void output_normals(Ioss::Region &region, Ioss::Region &output_region, bool reverse_normals)
  {
    // Define output fields...
    {
      const auto *v3d = Ioss::VariableType::factory("vector_3d");
      output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
      auto       *nb        = (*output_region.get_node_blocks().begin());
      int         num_nodes = nb->entity_count();
      Ioss::Field node_normal("node_normal", Ioss::Field::REAL, v3d, Ioss::Field::TRANSIENT,
                              num_nodes);
      nb->field_add(std::move(node_normal));

      // Iterate over the element blocks and calculate both node normals and face normals...
      const auto &ebs = output_region.get_element_blocks();
      for (auto *eb : ebs) {
        int         num_elem = eb->entity_count();
        Ioss::Field face_normal("face_normal", Ioss::Field::REAL, v3d, Ioss::Field::TRANSIENT,
                                num_elem);
        eb->field_add(std::move(face_normal));
      }
      output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }

    output_region.begin_mode(Ioss::STATE_TRANSIENT);
    int ostep = output_region.add_state(0.0);
    output_region.begin_state(ostep);

    auto *nb  = (*region.get_node_blocks().begin());
    auto *nbo = (*output_region.get_node_blocks().begin());

    // Get the nodal coordinates...
    int num_nodes  = nb->entity_count();
    int coord_size = nb->get_field("mesh_model_coordinates").get_size();

    std::vector<double> coord(3 * num_nodes);
    assert(3 * num_nodes * sizeof(double) == (size_t)coord_size);

    nb->get_field_data("mesh_model_coordinates", Data(coord), coord_size);

    // Also get an array for the average nodal normal vector...
    std::vector<double> node_normal(3 * num_nodes);
    std::fill(node_normal.begin(), node_normal.end(), 0.0);

    // Iterate over the element blocks and calculate both node normals and face normals...
    std::vector<int>    conn;
    std::vector<double> face_normal;
    const auto         &ebs = region.get_element_blocks();
    for (const auto *eb : ebs) {
      const std::string &name = (*eb).name();

      auto *ebo = output_region.get_element_block(name);
      if (ebo == nullptr) {
        std::cerr << "INTERNAL ERROR: Could not find element block named '" << name << "'\n";
        std::exit(EXIT_FAILURE);
      }

      int num_elem          = eb->entity_count();
      int num_node_per_elem = eb->topology()->number_nodes();

      // Get the connectivity array...
      conn.resize(num_elem * num_node_per_elem);
      eb->get_field_data("connectivity", conn);

      // Connectivity is in global id space; change to local...
      for (int i = 0; i < num_elem * num_node_per_elem; i++) {
        int local = region.node_global_to_local(conn[i]);
        conn[i]   = local - 1;
      }
      // Get an array for storing the face normals for this block...
      face_normal.resize(3 * num_elem);
      std::fill(face_normal.begin(), face_normal.end(), 0.0);

      calculate_normals(node_normal, face_normal, num_elem, num_node_per_elem, coord, conn,
                        reverse_normals);
      ebo->put_field_data("face_normal", face_normal);
    }

    int nsize = node_normal.size();
    for (int i = 0; i < nsize; i += 3) {
      vector3d a(node_normal[i + 0], node_normal[i + 1], node_normal[i + 2]);
      a.normalize();

      node_normal[i + 0] = a.x;
      node_normal[i + 1] = a.y;
      node_normal[i + 2] = a.z;
    }
    nbo->put_field_data("node_normal", node_normal);
    output_region.end_state(ostep);
    output_region.end_mode(Ioss::STATE_TRANSIENT);
  }

  void calculate_normals(std::vector<double> &node_normal, std::vector<double> &face_normal,
                         int num_elem, int num_node_per_elem, std::vector<double> &coord,
                         std::vector<int> &conn, bool reverse_normals)
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

        face_normal[iel * 3 + 0] = plnorm.x;
        face_normal[iel * 3 + 1] = plnorm.y;
        face_normal[iel * 3 + 2] = plnorm.z;

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
          local[i].set(coord[node * 3 + 0], coord[node * 3 + 1], coord[node * 3 + 2]);
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

          face_normal[iel * 3 + 0] += a.x;
          face_normal[iel * 3 + 1] += a.y;
          face_normal[iel * 3 + 2] += a.z;

          int node = conn[ioff + i];
          node_normal[node * 3 + 0] += a.x;
          node_normal[node * 3 + 1] += a.y;
          node_normal[node * 3 + 2] += a.z;
        }
      }
      vector3d a(face_normal[iel * 3 + 0], face_normal[iel * 3 + 1], face_normal[iel * 3 + 2]);
      a.normalize();
      face_normal[iel * 3 + 0] = a.x;
      face_normal[iel * 3 + 1] = a.y;
      face_normal[iel * 3 + 2] = a.z;
    }
  }

  void show_step(int istep, double time)
  {
    std::cerr.setf(std::ios::scientific);
    std::cerr.setf(std::ios::showpoint);
    std::cerr << "     Time step " << std::setw(5) << istep << " at time " << std::setprecision(5)
              << time << '\n';
  }
} // namespace
