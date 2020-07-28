/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "EML_CodeTypes.h"
#include "EML_SystemInterface.h"

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_SubSystem.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Utils.h>

#include "add_to_log.h"
#include "fmt/ostream.h"

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

// ========================================================================
namespace {
  bool file_info(const std::string &inpfile, const std::string &input_type,
                 SystemInterface &interFace);

  void output_names(const std::string &type, const Ioss::NameList &fields,
                    Ioss::GroupingEntity *entity)
  {
    fmt::print("\n{} variables on exodus data base:\n", type);
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string               field_name = *IF;
      const Ioss::VariableType *var_type   = entity->get_field(field_name).raw_storage();
      fmt::print("\t{}\t{}\n", field_name, var_type->name());
    }
  }
} // namespace
// ========================================================================

namespace {
  std::string codename;
  std::string version = "1.1";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  time_t      begin_time = time(nullptr);
  std::string in_type    = "exodusII";

  bool ok    = false;
  codename   = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  try {
    SystemInterface::show_version();
    Ioss::Init::Initializer io;

    SystemInterface interFace;
    ok = interFace.parse_options(argc, argv);

    if (ok) {
      std::string in_file     = interFace.input_file();
      std::string output_file = interFace.output_file();

      fmt::print("Input:    '{}', Type: {}\n", in_file, in_type);
      fmt::print("Output:   '{}', Type: matlab script\n\n", output_file);

      ok = file_info(in_file, in_type, interFace);
    }
    std::string success = ok ? "successful" : "unsuccessful";
    fmt::print("\n{} execution {}.\n", codename, success);
  }
  catch (std::exception &e) {
    fmt::print("ERROR: (EXOMATLAB) Standard exception: {}\n", e.what());
  }
  time_t end_time = time(nullptr);
  add_to_log(codename.c_str(), (int)(end_time - begin_time));
#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

namespace {
  bool file_info(const std::string &inpfile, const std::string &input_type,
                 SystemInterface &interFace)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi =
        Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART, (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == nullptr || !dbi->ok(true)) {
      return false;
    }

    dbi->set_field_separator(interFace.field_suffix());
    dbi->set_lower_case_variable_names(false);

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    if (interFace.list_vars()) {
      StringIdVector types_to_list = interFace.vars_to_list();
      for (auto types : types_to_list) {
        std::string type = types.first;

        if (type == "all" || type == "global") {
          Ioss::NameList fields;
          region.field_describe(Ioss::Field::TRANSIENT, &fields);
          output_names("Global", fields, &region);
        }
        if (type == "all" || type == "nodal") {
          Ioss::NameList   fields;
          Ioss::NodeBlock *nb = region.get_node_blocks()[0];
          nb->field_describe(Ioss::Field::TRANSIENT, &fields);
          output_names("Nodal", fields, nb);
        }
      }
      return true;
    }

    Ioss::NameList fields;
    StringIdVector global_vars = interFace.global_var_names();
    if (!global_vars.empty()) {
      if (global_vars[0].first == "all") {
        region.field_describe(Ioss::Field::TRANSIENT, &fields);
      }
      else if (global_vars[0].first == "none") {
        ; // do nothing.  This will be used when nodal, element, ... supported
      }
      else {
        for (auto &global_var : global_vars) {
          std::string field_name = global_var.first;
          if (region.field_exists(field_name)) {
            fields.push_back(field_name);
          }
          else {
            fmt::print("WARNING: Global variable named '{}' does not exist; it will be skipped.\n",
                       field_name);
            ;
          }
        }
      }
    }
    else {
      region.field_describe(Ioss::Field::TRANSIENT, &fields);
    }

    if (fields.empty()) {
      fmt::print("No variables selected; no output will be written\n");
      return false;
    }

    std::ofstream out_stream;
    out_stream.open(interFace.output_file().c_str());

    out_stream.setf(std::ios::scientific);
    out_stream.setf(std::ios::showpoint);

    fmt::print(out_stream, "% number of curves\nnvars = {}\n", fields.size() + 1);

    size_t                         namelen = 4; // length of 'time'
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      if (field_name.length() > namelen) {
        namelen = field_name.length();
      }
    }

    fmt::print(out_stream, "names= [\n'{:<{}}';\n", "TIME", namelen);

    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      fmt::print(out_stream, "'{:<{}}';\n", field_name, namelen);
    }
    fmt::print(out_stream, "];\n");

    // Get number of timesteps...
    int num_steps = 0;
    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      num_steps = region.get_property("state_count").get_int();
    }
    else {
      fmt::print(out_stream, "GENESIS file -- no time steps written\n");
      return false;
    }

    // ========================================================================
    // Calculate min and max times to extract data...
    int st_min = 1;
    int st_max = num_steps;

    double tmin = interFace.minimum_time();
    double tmax = interFace.maximum_time();
    if (tmax == -1.0) {
      tmax = region.get_max_time().second;
    }

    if (tmin > 0.0 || tmax != region.get_max_time().second) {
      st_min = 0;

      // Count number of active steps...
      for (int i = 1; i <= num_steps; i++) {
        double time = region.get_state_time(i);
        if (time < tmin) {
          st_min = i;
        }
        if (time <= tmax) {
          st_max = i;
        }
      }
      st_min++;
      num_steps = st_max - st_min + 1;
    }

    // ========================================================================
    // Output time values...

    fmt::print(out_stream, "TIME = zeros(1, {});\n", num_steps);
    fmt::print(out_stream, "TIME= [\n");
    for (int i = st_min; i <= st_max; i++) {
      double time = region.get_state_time(i);
      fmt::print(out_stream, " {:13.6e}\n", time);
    }
    fmt::print(out_stream, "];\n");

    // ========================================================================
    // Output field values...

    std::vector<double> data(1);
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      int         comp_count = region.get_field(field_name).raw_storage()->component_count();
      fmt::print(out_stream, "{} = zeros({}, {});\n", field_name, comp_count, num_steps);
      fmt::print(out_stream, "{} = [\n", field_name);
      for (int istep = st_min; istep <= st_max; istep++) {
        region.begin_state(istep);
        region.get_field_data(field_name, data);
        for (int i = 0; i < comp_count; i++) {
          fmt::print(out_stream, " {:13.6e}", data[i]);
        }
        fmt::print(out_stream, ";\n");
        region.end_state(istep);
      }
      fmt::print(out_stream, "];\n");
    }
    fmt::print("Wrote data for {} variables at {} timesteps.\n", fields.size(),
               st_max - st_min + 1);
    return true;
  }
} // namespace
