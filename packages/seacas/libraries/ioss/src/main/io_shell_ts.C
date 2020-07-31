// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_FileInfo.h>
#include <Ioss_MeshType.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_ScopeGuard.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_SubSystem.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <stddef.h>
#include <stdlib.h>
#include <string>
#ifndef _MSC_VER
#include <sys/times.h>
#endif
#include <unistd.h>
#include <vector>

#include "shell_interface.h"

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

#ifdef SEACAS_HAVE_KOKKOS
#include <Kokkos_Core.hpp> // for Kokkos::View
#endif

#define DO_OUTPUT                                                                                  \
  if (rank == 0)                                                                                   \
  std::cerr

// ========================================================================

namespace {

  struct my_numpunct : std::numpunct<char>
  {
  protected:
    char        do_thousands_sep() const { return ','; }
    std::string do_grouping() const { return "\3"; }
  };

  struct DataPool
  {
    // Data space shared by most field input/output routines...
    std::vector<char>    data;
    std::vector<int>     data_int;
    std::vector<int64_t> data_int64;
    std::vector<double>  data_double;
    std::vector<Complex> data_complex;
#ifdef SEACAS_HAVE_KOKKOS
    Kokkos::View<char *>    data_view_char;
    Kokkos::View<int *>     data_view_int;
    Kokkos::View<int64_t *> data_view_int64;
    Kokkos::View<double *>  data_view_double;
    // Kokkos::View<Kokkos_Complex *> data_view_complex cannot be a global variable,
    // Since Kokkos::initialize() has not yet been called. Also, a Kokkos:View cannot
    // have type std::complex entities.
    Kokkos::View<char **>    data_view_2D_char;
    Kokkos::View<int **>     data_view_2D_int;
    Kokkos::View<int64_t **> data_view_2D_int64;
    Kokkos::View<double **>  data_view_2D_double;
    // Kokkos::View<Kokkos_Complex **> data_view_2D_complex cannot be a global variable,
    // Since Kokkos::initialize() has not yet been called. Also, a Kokkos:View cannot
    // have type std::complex entities.
    Kokkos::View<char **, Kokkos::LayoutRight, Kokkos::HostSpace> data_view_2D_char_layout_space;
    Kokkos::View<int **, Kokkos::LayoutRight, Kokkos::HostSpace>  data_view_2D_int_layout_space;
    Kokkos::View<int64_t **, Kokkos::LayoutRight, Kokkos::HostSpace>
        data_view_2D_int64_layout_space;
    Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace>
        data_view_2D_double_layout_space;
    // Kokkos::View<Kokkos_Complex **, Kokkos::LayoutRight, Kokkos::HostSpace>
    // data_view_2D_complex_layout_space cannot be a global variable,
    // Since Kokkos::initialize() has not yet been called. Also, a Kokkos:View cannot
    // have type std::complex entities.
#endif
  };

  int  rank      = 0;
  bool mem_stats = false;

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
  void transfer_fields(const std::vector<T *> &entities, Ioss::Region &output_region,
                       Ioss::Field::RoleType role, const IOShell::Interface &interFace);

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix = "");

  template <typename T>
  void transfer_field_data(const std::vector<T *> &entities, Ioss::Region &output_region,
                           Ioss::Field::RoleType role, const IOShell::Interface &interFace);

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const IOShell::Interface &interFace,
                           const std::string &prefix = "");

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge);

  void transfer_qa_info(Ioss::Region &in, Ioss::Region &out);

  void transform_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                        Ioss::Field::RoleType role);

  void transform_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                            Ioss::Field::RoleType role, const IOShell::Interface &interFace);
  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &       field_name,
                                    const IOShell::Interface &interFace);

  void file_copy(IOShell::Interface &interFace);

  template <typename INT>
  void set_owned_node_count(Ioss::Region &region, int my_processor, INT dummy);
} // namespace
// ========================================================================

namespace {
  std::string codename;
  std::string version = "4.7";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  std::cerr.imbue(std::locale(std::locale(), new my_numpunct));

#ifdef SEACAS_HAVE_KOKKOS
  Kokkos::ScopeGuard kokkos(argc, argv);
#endif

  IOShell::Interface interFace;
  bool               success = interFace.parse_options(argc, argv);
  if (!success) {
    exit(EXIT_FAILURE);
  }

  Ioss::SerializeIO::setGroupFactor(interFace.serialize_io_size);
  mem_stats = interFace.memory_statistics;

  Ioss::Init::Initializer io;

  std::string in_file  = interFace.inputFile[0];
  std::string out_file = interFace.outputFile;

  DO_OUTPUT << "Input:    '" << in_file << "', Type: " << interFace.inFiletype << '\n';
  DO_OUTPUT << "Output:   '" << out_file << "', Type: " << interFace.outFiletype << '\n';
  DO_OUTPUT << '\n';

#ifdef SEACAS_HAVE_KOKKOS
  DO_OUTPUT << "Kokkos default execution space configuration:\n";
  Kokkos::DefaultExecutionSpace::print_configuration(std::cerr, false);
  DO_OUTPUT << "\n";
#endif

  file_copy(interFace);

  if (mem_stats) {
    int64_t MiB = 1024 * 1024;
#ifdef SEACAS_HAVE_MPI
    int64_t             min, max, avg;
    Ioss::ParallelUtils parallel(MPI_COMM_WORLD);
    parallel.memory_stats(min, max, avg);
    DO_OUTPUT << "\n\tCurrent Memory: " << min / MiB << "M  " << max / MiB << "M  " << avg / MiB
              << "M\n";

    parallel.hwm_memory_stats(min, max, avg);
    DO_OUTPUT << "\n\tHigh Water Memory: " << min / MiB << "M  " << max / MiB << "M  " << avg / MiB
              << "M\n";
#else
    int64_t mem = Ioss::Utils::get_memory_info();
    int64_t hwm = Ioss::Utils::get_hwm_memory_info();
    DO_OUTPUT << "\n\tCurrent Memory:    " << mem / MiB << "M\n"
              << "\n\tHigh Water Memory: " << hwm / MiB << "M\n";
#endif
  }
  DO_OUTPUT << "\n" << codename << " execution successful.\n";

  return EXIT_SUCCESS;
}

namespace {
  void file_copy(IOShell::Interface &interFace)
  {
    Ioss::PropertyManager properties;
    if (interFace.ints_64_bit) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interFace.ints_32_bit) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 4));
    }

    if (interFace.reals_32_bit) {
      properties.add(Ioss::Property("REAL_SIZE_DB", 4));
    }

    if (interFace.in_memory_read) {
      properties.add(Ioss::Property("MEMORY_READ", 1));
    }

    if (interFace.in_memory_write) {
      properties.add(Ioss::Property("MEMORY_WRITE", 1));
    }

    if (interFace.compression_level > 0 || interFace.shuffle) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", static_cast<int>(interFace.shuffle)));
    }

    if (interFace.compose_output == "default") {
      if (interFace.outFiletype == "cgns") {
        properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
        properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
      }
      else {
        properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
        properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
      }
    }
    else if (interFace.compose_output == "external") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
    }
    else if (interFace.compose_output != "none") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
    }

    if (interFace.netcdf4) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    if (interFace.inputFile.size() > 1) {
      properties.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
    }

    if (interFace.debug) {
      properties.add(Ioss::Property("LOGGING", 1));
    }

    if (interFace.memory_statistics) {
      properties.add(Ioss::Property("ENABLE_TRACING", true));
    }

    if (!interFace.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interFace.decomp_method));
    }
    bool first = true;
    for (const auto &inpfile : interFace.inputFile) {
      Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(
          interFace.inFiletype, inpfile, Ioss::READ_MODEL, (MPI_Comm)MPI_COMM_WORLD, properties);
      if (dbi == nullptr || !dbi->ok(true)) {
        std::exit(EXIT_FAILURE);
      }
      dbi->set_parallel_consistency(false);

      if (mem_stats) {
        dbi->progress("Database Creation");
      }
      if (!interFace.lower_case_variable_names) {
        dbi->set_lower_case_variable_names(false);
      }
      dbi->set_surface_split_type(Ioss::int_to_surface_split(interFace.surface_split_type));
      dbi->set_field_separator(interFace.fieldSuffixSeparator);
      if (interFace.ints_64_bit) {
        dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
      }

      if (!interFace.groupName.empty()) {
        bool success = dbi->open_group(interFace.groupName);
        if (!success) {
          DO_OUTPUT << "ERROR: Unable to open group '" << interFace.groupName << "' in file '"
                    << inpfile << "\n";
          return;
        }
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region region(dbi, "region_1");

      if (region.mesh_type() != Ioss::MeshType::UNSTRUCTURED) {
        DO_OUTPUT << "\nERROR: io_shell does not support '" << region.mesh_type_string()
                  << "' meshes.  Only 'Unstructured' mesh is supported at this time.\n";
        return;
      }

      // Get length of longest name on input file...
      int max_name_length = dbi->maximum_symbol_length();
      if (max_name_length > 0) {
        properties.add(Ioss::Property("MAXIMUM_NAME_LENGTH", max_name_length));
      }

      // Get integer size being used on the input file and propgate
      // to output file...
      int int_byte_size_api = dbi->int_byte_size_api();
      if (!properties.exists("INTEGER_SIZE_API")) {
        if (interFace.ints_32_bit) {
          properties.add(Ioss::Property("INTEGER_SIZE_DB", 4));
        }
        else {
          properties.add(Ioss::Property("INTEGER_SIZE_DB", int_byte_size_api));
        }
        properties.add(Ioss::Property("INTEGER_SIZE_API", int_byte_size_api));
      }
      if (int_byte_size_api == 8) {
        interFace.ints_64_bit = true;
      }
      //========================================================================
      // OUTPUT ...
      //========================================================================
      Ioss::DatabaseIO *dbo =
          Ioss::IOFactory::create(interFace.outFiletype, interFace.outputFile, Ioss::WRITE_RESTART,
                                  (MPI_Comm)MPI_COMM_WORLD, properties);
      if (dbo == nullptr || !dbo->ok(true)) {
        std::exit(EXIT_FAILURE);
      }
      dbo->set_parallel_consistency(false);

      // NOTE: 'output_region' owns 'dbo' pointer at this time
      Ioss::Region output_region(dbo, "region_2");
      // Set the qa information...
      output_region.property_add(Ioss::Property(std::string("code_name"), codename));
      output_region.property_add(Ioss::Property(std::string("code_version"), version));

      if (interFace.inputFile.size() > 1) {
        properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));

        if (!first) {
          // Putting each file into its own output group...
          // The name of the group will be the basename portion of the filename...
          Ioss::FileInfo file(inpfile);
          dbo->create_subgroup(file.tailname());
        }
        else {
          first = false;
        }
      }

      if (interFace.debug) {
        DO_OUTPUT << "DEFINING MODEL ... \n";
      }
      if (mem_stats) {
        dbi->progress("DEFINING MODEL");
      }
      if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
        DO_OUTPUT << "ERROR: Could not put output region into define model state\n";
        std::exit(EXIT_FAILURE);
      }

      // Get all properties of input database...
      transfer_properties(&region, &output_region);
      transfer_qa_info(region, output_region);

      transfer_nodeblock(region, output_region, interFace.debug);

#ifdef SEACAS_HAVE_MPI
      // This also assumes that the node order and count is the same for input
      // and output regions... (This is checked during nodeset output)
      if (output_region.get_database()->needs_shared_node_information()) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (interFace.ints_64_bit)
          set_owned_node_count(region, rank, (int64_t)0);
        else
          set_owned_node_count(region, rank, (int)0);
      }
#endif

      transfer_edgeblocks(region, output_region, interFace.debug);
      transfer_faceblocks(region, output_region, interFace.debug);
      transfer_elementblocks(region, output_region, interFace.debug);

      transfer_nodesets(region, output_region, interFace.debug);
      transfer_edgesets(region, output_region, interFace.debug);
      transfer_facesets(region, output_region, interFace.debug);
      transfer_elemsets(region, output_region, interFace.debug);

      transfer_sidesets(region, output_region, interFace.debug);
      transfer_commsets(region, output_region, interFace.debug);

      transfer_coordinate_frames(region, output_region, interFace.debug);

      if (interFace.debug) {
        DO_OUTPUT << "END STATE_DEFINE_MODEL... " << '\n';
      }
      if (mem_stats) {
        dbi->progress("END STATE_DEFINE_MODEL");
      }

      output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

      if (interFace.debug) {
        DO_OUTPUT << "TRANSFERRING MESH FIELD DATA ... " << '\n';
      }
      if (mem_stats) {
        dbi->progress("TRANSFERRING MESH FIELD DATA ... ");
      }

      // Model defined, now fill in the model data...
      output_region.begin_mode(Ioss::STATE_MODEL);

      // Transfer MESH field_data from input to output...
      transfer_field_data(region.get_node_blocks(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_node_blocks(), output_region, Ioss::Field::ATTRIBUTE,
                          interFace);

      transfer_field_data(region.get_edge_blocks(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_edge_blocks(), output_region, Ioss::Field::ATTRIBUTE,
                          interFace);

      transfer_field_data(region.get_face_blocks(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_face_blocks(), output_region, Ioss::Field::ATTRIBUTE,
                          interFace);

      transfer_field_data(region.get_element_blocks(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_element_blocks(), output_region, Ioss::Field::ATTRIBUTE,
                          interFace);

      transfer_field_data(region.get_nodesets(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_nodesets(), output_region, Ioss::Field::ATTRIBUTE, interFace);

      transfer_field_data(region.get_edgesets(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_edgesets(), output_region, Ioss::Field::ATTRIBUTE, interFace);

      transfer_field_data(region.get_facesets(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_facesets(), output_region, Ioss::Field::ATTRIBUTE, interFace);

      transfer_field_data(region.get_elementsets(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_elementsets(), output_region, Ioss::Field::ATTRIBUTE,
                          interFace);

      transfer_field_data(region.get_commsets(), output_region, Ioss::Field::MESH, interFace);
      transfer_field_data(region.get_commsets(), output_region, Ioss::Field::ATTRIBUTE, interFace);
      transfer_field_data(region.get_commsets(), output_region, Ioss::Field::COMMUNICATION,
                          interFace);

      // Side Sets
      {
        const auto &fss = region.get_sidesets();
        for (const auto &ifs : fss) {
          const std::string &name = ifs->name();
          if (interFace.debug) {
            DO_OUTPUT << name << ", ";
          }
          // Find matching output sideset
          Ioss::SideSet *ofs = output_region.get_sideset(name);

          if (ofs != nullptr) {
            transfer_field_data(ifs, ofs, Ioss::Field::MESH, interFace);
            transfer_field_data(ifs, ofs, Ioss::Field::ATTRIBUTE, interFace);

            const auto &fbs = ifs->get_side_blocks();
            for (const auto &ifb : fbs) {

              // Find matching output sideblock
              const std::string &fbname = ifb->name();
              if (interFace.debug) {
                DO_OUTPUT << fbname << ", ";
              }
              Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

              if (ofb != nullptr) {
                transfer_field_data(ifb, ofb, Ioss::Field::MESH, interFace);
                transfer_field_data(ifb, ofb, Ioss::Field::ATTRIBUTE, interFace);
              }
            }
          }
        }
        if (interFace.debug) {
          DO_OUTPUT << '\n';
        }
      }
      if (interFace.debug) {
        DO_OUTPUT << "END STATE_MODEL... " << '\n';
      }
      if (mem_stats) {
        dbi->progress("END STATE_MODEL... ");
      }
      output_region.end_mode(Ioss::STATE_MODEL);

      if (interFace.delete_timesteps) {
        if (mem_stats) {
          dbi->progress("Prior to Memory Released... ");
          dbi->release_memory();
          dbo->release_memory();
          dbi->progress("Memory Released... ");
        }
        return;
      }

      if (interFace.debug) {
        DO_OUTPUT << "DEFINING TRANSIENT FIELDS ... " << '\n';
      }
      if (mem_stats) {
        dbi->progress("DEFINING TRANSIENT FIELDS ... ");
      }

      if (region.property_exists("state_count") &&
          region.get_property("state_count").get_int() > 0) {
        if (!interFace.debug) {
          DO_OUTPUT << "\n Number of time steps on database     =" << std::setw(12)
                    << region.get_property("state_count").get_int() << "\n\n";
        }

        output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

        // For each 'TRANSIENT' field in the node blocks and element
        // blocks, transfer to the output node and element blocks.
        transfer_fields(&region, &output_region, Ioss::Field::TRANSIENT);

        transfer_fields(region.get_node_blocks(), output_region, Ioss::Field::TRANSIENT, interFace);
        transfer_fields(region.get_edge_blocks(), output_region, Ioss::Field::TRANSIENT, interFace);
        transfer_fields(region.get_face_blocks(), output_region, Ioss::Field::TRANSIENT, interFace);
        transfer_fields(region.get_element_blocks(), output_region, Ioss::Field::TRANSIENT,
                        interFace);

        transfer_fields(region.get_nodesets(), output_region, Ioss::Field::TRANSIENT, interFace);
        transfer_fields(region.get_edgesets(), output_region, Ioss::Field::TRANSIENT, interFace);
        transfer_fields(region.get_facesets(), output_region, Ioss::Field::TRANSIENT, interFace);
        transfer_fields(region.get_elementsets(), output_region, Ioss::Field::TRANSIENT, interFace);

        // Side Sets
        {
          const auto &fss = region.get_sidesets();
          for (const auto &ifs : fss) {
            const std::string &name = ifs->name();
            if (interFace.debug) {
              DO_OUTPUT << name << ", ";
            }

            // Find matching output sideset
            Ioss::SideSet *ofs = output_region.get_sideset(name);
            if (ofs != nullptr) {
              transfer_fields(ifs, ofs, Ioss::Field::TRANSIENT);

              const auto &fbs = ifs->get_side_blocks();
              for (const auto &ifb : fbs) {

                // Find matching output sideblock
                const std::string &fbname = ifb->name();
                if (interFace.debug) {
                  DO_OUTPUT << fbname << ", ";
                }

                Ioss::SideBlock *ofb = ofs->get_side_block(fbname);
                if (ofb != nullptr) {
                  transfer_fields(ifb, ofb, Ioss::Field::TRANSIENT);
                }
              }
            }
          }
          if (interFace.debug) {
            DO_OUTPUT << '\n';
          }
        }
        if (interFace.debug) {
          DO_OUTPUT << "END STATE_DEFINE_TRANSIENT... " << '\n';
        }
        if (mem_stats) {
          dbi->progress("END STATE_DEFINE_TRANSIENT... ");
        }
        output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
      }

      if (interFace.debug) {
        DO_OUTPUT << "TRANSFERRING TRANSIENT FIELDS ... " << '\n';
      }
      if (mem_stats) {
        dbi->progress("TRANSFERRING TRANSIENT FIELDS... ");
      }

      output_region.begin_mode(Ioss::STATE_TRANSIENT);
      // Get the timesteps from the input database.  Step through them
      // and transfer fields to output database...

      int step_count = region.get_property("state_count").get_int();

      for (int istep = 1; istep <= step_count; istep++) {
        double time = region.get_state_time(istep);
        if (time < interFace.minimum_time) {
          continue;
        }
        if (interFace.maximum_time != 0.0 && time > interFace.maximum_time) {
          break;
        }

        int ostep = output_region.add_state(time);
        show_step(istep, time);

        output_region.begin_state(ostep);
        region.begin_state(istep);

        transfer_field_data(&region, &output_region, Ioss::Field::TRANSIENT, interFace);

        transfer_field_data(region.get_node_blocks(), output_region, Ioss::Field::TRANSIENT,
                            interFace);
        transfer_field_data(region.get_edge_blocks(), output_region, Ioss::Field::TRANSIENT,
                            interFace);
        transfer_field_data(region.get_face_blocks(), output_region, Ioss::Field::TRANSIENT,
                            interFace);
        transfer_field_data(region.get_element_blocks(), output_region, Ioss::Field::TRANSIENT,
                            interFace);

        transfer_field_data(region.get_nodesets(), output_region, Ioss::Field::TRANSIENT,
                            interFace);
        transfer_field_data(region.get_edgesets(), output_region, Ioss::Field::TRANSIENT,
                            interFace);
        transfer_field_data(region.get_facesets(), output_region, Ioss::Field::TRANSIENT,
                            interFace);
        transfer_field_data(region.get_elementsets(), output_region, Ioss::Field::TRANSIENT,
                            interFace);

        // Side Sets
        {
          const auto &fss = region.get_sidesets();
          for (const auto &ifs : fss) {
            const std::string &name = ifs->name();
            if (interFace.debug) {
              DO_OUTPUT << name << ", ";
            }

            // Find matching output sideset
            Ioss::SideSet *ofs = output_region.get_sideset(name);
            if (ofs != nullptr) {
              transfer_field_data(ifs, ofs, Ioss::Field::TRANSIENT, interFace);

              const auto &fbs = ifs->get_side_blocks();
              for (const auto &ifb : fbs) {

                // Find matching output sideblock
                const std::string &fbname = ifb->name();
                if (interFace.debug) {
                  DO_OUTPUT << fbname << ", ";
                }

                Ioss::SideBlock *ofb = ofs->get_side_block(fbname);
                if (ofb != nullptr) {
                  transfer_field_data(ifb, ofb, Ioss::Field::TRANSIENT, interFace);
                }
              }
            }
          }
        }
        region.end_state(istep);
        output_region.end_state(ostep);
      }
      if (interFace.debug) {
        DO_OUTPUT << "END STATE_TRANSIENT... " << '\n';
      }
      if (mem_stats) {
        dbi->progress("END STATE_TRANSIENT ... ");
      }
      output_region.end_mode(Ioss::STATE_TRANSIENT);

      if (mem_stats) {
        dbi->progress("Prior to Memory Released... ");
        dbi->release_memory();
        dbo->release_memory();
        dbi->progress("Memory Released... ");
      }
    }
  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nbs = region.get_node_blocks();
    size_t      id  = 1;
    for (const auto &inb : nbs) {
      const std::string &name = inb->name();
      if (debug) {
        DO_OUTPUT << name << ", ";
      }
      size_t num_nodes = inb->entity_count();
      size_t degree    = inb->get_property("component_degree").get_int();
      if (!debug) {
        DO_OUTPUT << " Number of coordinates per node       =" << std::setw(12) << degree << "\n";
        DO_OUTPUT << " Number of nodes                      =" << std::setw(12) << num_nodes
                  << "\n";
      }

      auto nb = new Ioss::NodeBlock(output_region.get_database(), name, num_nodes, degree);
      output_region.add(nb);

      transfer_properties(inb, nb);

      if (output_region.get_database()->needs_shared_node_information()) {
        // If the "owning_processor" field exists on the input
        // nodeblock, transfer it and the "ids" field to the output
        // nodeblock at this time since it is used to determine
        // per-processor sizes of nodeblocks and nodesets.
        if (inb->field_exists("owning_processor")) {
          size_t            isize = inb->get_field("ids").get_size();
          std::vector<char> data(isize);
          inb->get_field_data("ids", data.data(), isize);
          nb->put_field_data("ids", data.data(), isize);
          isize = inb->get_field("owning_processor").get_size();
          data.resize(isize);
          inb->get_field_data("owning_processor", data.data(), isize);
          nb->put_field_data("owning_processor", data.data(), isize);
        }
      }

      transfer_fields(inb, nb, Ioss::Field::MESH);
      transfer_fields(inb, nb, Ioss::Field::ATTRIBUTE);
      ++id;
    }
    if (debug) {
      DO_OUTPUT << '\n';
    }
  }

#if 0
  template <typename T>
  void transfer_fields(const std::vector<T *> &entities, Ioss::Region &output_region,
                       Ioss::Field::RoleType role, const IOShell::Interface &interFace)
  {
    for (const auto &entity : entities) {
      const std::string &name = entity->name();
      if (interFace.debug) {
        DO_OUTPUT << name << ", ";
      }

      // Find the corresponding output entity...
      Ioss::GroupingEntity *oeb = output_region.get_entity(name, entity->type());
      if (oeb != nullptr) {
        transfer_fields(entity, oeb, role);
        if (interFace.do_transform_fields) {
          transform_fields(entity, oeb, role);
        }
      }
    }
    if (interFace.debug) {
      DO_OUTPUT << '\n';
    }
  }
#endif

  struct param
  {
    Ioss::GroupingEntity *    entity;
    Ioss::Region *            output_region;
    Ioss::Field::RoleType     role;
    const IOShell::Interface *interFace;
  };

  void *transfer_fields_ts(void *varg)
  {
    param *arg           = static_cast<param *>(varg);
    auto   entity        = arg->entity;
    auto   output_region = arg->output_region;
    auto   interFace     = arg->interFace;
    auto   role          = arg->role;

    const std::string &name = entity->name();
    if (interFace->debug) {
      DO_OUTPUT << name << ", ";
    }

    // Find the corresponding output entity...
    Ioss::GroupingEntity *oeb = output_region->get_entity(name, entity->type());
    if (oeb != nullptr) {
      transfer_fields(entity, oeb, role);
      if (interFace->do_transform_fields) {
        transform_fields(entity, oeb, role);
      }
    }
    if (interFace->debug) {
      DO_OUTPUT << '\n';
    }
    return varg;
  }

  template <typename T>
  void transfer_fields(const std::vector<T *> &entities, Ioss::Region &output_region,
                       Ioss::Field::RoleType role, const IOShell::Interface &interFace)
  {
    std::vector<pthread_t> threads(entities.size());
    std::vector<param>     params(entities.size());

    int t = 0;
    for (const auto &entity : entities) {
      params[t].entity        = entity;
      params[t].output_region = &output_region;
      params[t].interFace     = &interFace;
      params[t].role          = role;
      pthread_create(&threads[t], nullptr, transfer_fields_ts, (void *)(params.data() + t));
      t++;
    }

    for (t = 0; t < (int)threads.size(); t++) {
      pthread_join(threads[t], nullptr);
    }
  }

  void *transfer_field_data_ts(void *varg)
  {
    param *arg           = static_cast<param *>(varg);
    auto   entity        = arg->entity;
    auto   output_region = arg->output_region;
    auto   interFace     = arg->interFace;
    auto   role          = arg->role;

    const std::string &name = entity->name();

    // Find the corresponding output block...
    Ioss::GroupingEntity *output = output_region->get_entity(name, entity->type());
    if (output != nullptr) {
      transfer_field_data(entity, output, role, *interFace);
      if (interFace->do_transform_fields) {
        transform_field_data(entity, output, role, *interFace);
      }
    }
    return arg;
  }

  template <typename T>
  void transfer_field_data(const std::vector<T *> &entities, Ioss::Region &output_region,
                           Ioss::Field::RoleType role, const IOShell::Interface &interFace)
  {
    std::vector<pthread_t> threads(entities.size());
    std::vector<param>     params(entities.size());

    int t = 0;
    for (const auto &entity : entities) {
      params[t].entity        = entity;
      params[t].output_region = &output_region;
      params[t].interFace     = &interFace;
      params[t].role          = role;
      pthread_create(&threads[t], nullptr, transfer_field_data_ts, (void *)(params.data() + t));
      t++;
    }

    for (t = 0; t < (int)threads.size(); t++) {
      pthread_join(threads[t], nullptr);
    }
  }

  template <typename T>
  void transfer_blocks(const std::vector<T *> &blocks, Ioss::Region &output_region, bool debug)
  {
    if (!blocks.empty()) {
      size_t total_entities = 0;
      for (const auto &iblock : blocks) {
        const std::string &name = iblock->name();
        if (debug) {
          DO_OUTPUT << name << ", ";
        }
        std::string type  = iblock->get_property("topology_type").get_string();
        size_t      count = iblock->entity_count();
        total_entities += count;

        auto block = new T(output_region.get_database(), name, type, count);
        output_region.add(block);
        transfer_properties(iblock, block);
        transfer_fields(iblock, block, Ioss::Field::MESH);
        transfer_fields(iblock, block, Ioss::Field::ATTRIBUTE);
      }
      if (!debug) {
        DO_OUTPUT << " Number of " << std::setw(14) << (*blocks.begin())->type_string()
                  << "s            =" << std::setw(12) << blocks.size() << "\t"
                  << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      }
      else {
        DO_OUTPUT << '\n';
      }
    }
  }

  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &ebs = region.get_element_blocks();
    transfer_blocks(ebs, output_region, debug);
  }

  void transfer_edgeblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &ebs = region.get_edge_blocks();
    transfer_blocks(ebs, output_region, debug);
  }

  void transfer_faceblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &ebs = region.get_face_blocks();
    transfer_blocks(ebs, output_region, debug);
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &fss         = region.get_sidesets();
    size_t      total_sides = 0;
    for (const auto &ss : fss) {
      const std::string &name = ss->name();
      if (debug) {
        DO_OUTPUT << name << ", ";
      }

      auto        surf = new Ioss::SideSet(output_region.get_database(), name);
      const auto &fbs  = ss->get_side_blocks();
      for (const auto &fb : fbs) {
        const std::string &fbname = fb->name();
        if (debug) {
          DO_OUTPUT << fbname << ", ";
        }
        std::string fbtype   = fb->get_property("topology_type").get_string();
        std::string partype  = fb->get_property("parent_topology_type").get_string();
        size_t      num_side = fb->entity_count();
        total_sides += num_side;

        auto block =
            new Ioss::SideBlock(output_region.get_database(), fbname, fbtype, partype, num_side);
        surf->add(block);
        transfer_properties(fb, block);
        transfer_fields(fb, block, Ioss::Field::MESH);
        transfer_fields(fb, block, Ioss::Field::ATTRIBUTE);
      }
      transfer_properties(ss, surf);
      transfer_fields(ss, surf, Ioss::Field::MESH);
      transfer_fields(ss, surf, Ioss::Field::ATTRIBUTE);
      output_region.add(surf);
    }
    if (!debug) {
      DO_OUTPUT << " Number of        SideSets            =" << std::setw(12) << fss.size() << "\t"
                << "Number of element sides =" << std::setw(12) << total_sides << "\n";
    }
    else {
      DO_OUTPUT << '\n';
    }
  }

  template <typename T>
  void transfer_sets(const std::vector<T *> &sets, Ioss::Region &output_region, bool debug)
  {
    if (!sets.empty()) {
      size_t total_entities = 0;
      for (const auto &set : sets) {
        const std::string &name = set->name();
        if (debug) {
          DO_OUTPUT << name << ", ";
        }
        size_t count = set->entity_count();
        total_entities += count;
        auto o_set = new T(output_region.get_database(), name, count);
        output_region.add(o_set);
        transfer_properties(set, o_set);
        transfer_fields(set, o_set, Ioss::Field::MESH);
        transfer_fields(set, o_set, Ioss::Field::ATTRIBUTE);
      }

      if (!debug) {
        DO_OUTPUT << " Number of " << std::setw(14) << (*sets.begin())->type_string()
                  << "s            =" << std::setw(12) << sets.size() << "\t"
                  << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      }
      else {
        DO_OUTPUT << '\n';
      }
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nss = region.get_nodesets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_edgesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nss = region.get_edgesets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_facesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nss = region.get_facesets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_elemsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &nss = region.get_elementsets();
    transfer_sets(nss, output_region, debug);
  }

  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const auto &css = region.get_commsets();
    for (const auto &ics : css) {
      const std::string &name = ics->name();
      if (debug) {
        DO_OUTPUT << name << ", ";
      }
      std::string type  = ics->get_property("entity_type").get_string();
      size_t      count = ics->entity_count();
      auto        cs    = new Ioss::CommSet(output_region.get_database(), name, type, count);
      output_region.add(cs);
      transfer_properties(ics, cs);
      transfer_fields(ics, cs, Ioss::Field::MESH);
      transfer_fields(ics, cs, Ioss::Field::ATTRIBUTE);
      transfer_fields(ics, cs, Ioss::Field::COMMUNICATION);
    }
    if (debug) {
      DO_OUTPUT << '\n';
    }
  }

  void transfer_coordinate_frames(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const Ioss::CoordinateFrameContainer &cf = region.get_coordinate_frames();
    for (const auto &frame : cf) {
      output_region.add(frame);
    }
    if (debug) {
      DO_OUTPUT << '\n';
    }
  }

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output
    // database...  If a prefix is specified, only transfer fields
    // whose names begin with the prefix
    for (const auto &field_name : fields) {
      Ioss::Field field = ige->get_field(field_name);
      if (field_name != "ids" && !oge->field_exists(field_name) &&
          Ioss::Utils::substr_equal(prefix, field_name)) {
        // If the field does not already exist, add it to the output node block
        oge->field_add(field);
      }
    }
  }

  void transform_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                        Ioss::Field::RoleType role)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output database...
    for (const auto &field_name : fields) {
      std::string out_field_name = field_name + "_mag";
      if (!oge->field_exists(out_field_name)) {
        // If the field does not already exist, add it to the output node block
        Ioss::Field field = ige->get_field(field_name);
        Ioss::Field tr_field(out_field_name, field.get_type(), field.raw_storage(),
                             field.get_role(), field.raw_count());

        Ioss::Transform *transform = Iotr::Factory::create("vector magnitude");
        assert(transform != nullptr);
        tr_field.add_transform(transform);

        Ioss::Transform *max_transform = Iotr::Factory::create("absolute_maximum");
        assert(max_transform != nullptr);
        tr_field.add_transform(max_transform);

        oge->field_add(tr_field);
      }
    }
  }

  void transform_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                            Ioss::Field::RoleType role, const IOShell::Interface &interFace)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    ige->field_describe(role, &state_fields);
    // Iterate through mesh description fields and transfer to
    // output database...
    for (const auto &field_name : state_fields) {
      std::string out_field_name = field_name + "_mag";

      assert(oge->field_exists(out_field_name));

      int basic_type = ige->get_field(field_name).get_type();

      size_t isize = ige->get_field(field_name).get_size();
      size_t osize = oge->get_field(out_field_name).get_size();
      assert(isize == osize);

      DataPool pool;
      pool.data.resize(isize);
      switch (interFace.data_storage_type) {
      case 1: ige->get_field_data(field_name, pool.data.data(), isize); break;
      case 2:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          ige->get_field_data(field_name, pool.data.data(), isize);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          ige->get_field_data(field_name, pool.data_int);
        }
        else if (basic_type == Ioss::Field::INT64) {
          ige->get_field_data(field_name, pool.data_int64);
        }
        else if (basic_type == Ioss::Field::REAL) {
          ige->get_field_data(field_name, pool.data_double);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          ige->get_field_data(field_name, pool.data_complex);
        }
        else {
        }
        break;
#ifdef SEACAS_HAVE_KOKKOS
      case 3:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          ige->get_field_data<char>(field_name, pool.data_view_char);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          ige->get_field_data<int>(field_name, pool.data_view_int);
        }
        else if (basic_type == Ioss::Field::INT64) {
          ige->get_field_data<int64_t>(field_name, pool.data_view_int64);
        }
        else if (basic_type == Ioss::Field::REAL) {
          ige->get_field_data<double>(field_name, pool.data_view_double);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          // Since data_view_complex cannot be a global variable.
          ige->get_field_data(field_name, pool.data.data(), isize);
        }
        else {
        }
        break;
      case 4:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          ige->get_field_data<char>(field_name, pool.data_view_2D_char);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          ige->get_field_data<int>(field_name, pool.data_view_2D_int);
        }
        else if (basic_type == Ioss::Field::INT64) {
          ige->get_field_data<int64_t>(field_name, pool.data_view_2D_int64);
        }
        else if (basic_type == Ioss::Field::REAL) {
          ige->get_field_data<double>(field_name, pool.data_view_2D_double);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          // Since data_view_complex cannot be a global variable.
          ige->get_field_data(field_name, pool.data.data(), isize);
        }
        else {
        }
        break;
      case 5:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          ige->get_field_data<char, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_char_layout_space);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          ige->get_field_data<int, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_int_layout_space);
        }
        else if (basic_type == Ioss::Field::INT64) {
          ige->get_field_data<int64_t, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_int64_layout_space);
        }
        else if (basic_type == Ioss::Field::REAL) {
          ige->get_field_data<double, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_double_layout_space);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          // Since data_view_complex cannot be a global variable.
          ige->get_field_data(field_name, pool.data.data(), isize);
        }
        else {
        }
        break;
#endif
      default:
        if (field_name == "mesh_model_coordinates") {
          std::cerr << "data_storage option not recognized.";
        }
        return;
      }

      switch (interFace.data_storage_type) {
      case 1: oge->put_field_data(out_field_name, pool.data.data(), osize); break;
      case 2:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          oge->put_field_data(field_name, pool.data.data(), osize);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          oge->put_field_data(field_name, pool.data_int);
        }
        else if (basic_type == Ioss::Field::INT64) {
          oge->put_field_data(field_name, pool.data_int64);
        }
        else if (basic_type == Ioss::Field::REAL) {
          oge->put_field_data(field_name, pool.data_double);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          oge->put_field_data(field_name, pool.data_complex);
        }
        else {
        }
        break;
#ifdef SEACAS_HAVE_KOKKOS
      case 3:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          oge->put_field_data<char>(field_name, pool.data_view_char);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          oge->put_field_data<int>(field_name, pool.data_view_int);
        }
        else if (basic_type == Ioss::Field::INT64) {
          oge->put_field_data<int64_t>(field_name, pool.data_view_int64);
        }
        else if (basic_type == Ioss::Field::REAL) {
          oge->put_field_data<double>(field_name, pool.data_view_double);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          // Since data_view_complex cannot be a global variable.
          oge->put_field_data(out_field_name, pool.data.data(), osize);
        }
        else {
        }
        break;
      case 4:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          oge->put_field_data<char>(field_name, pool.data_view_2D_char);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          oge->put_field_data<int>(field_name, pool.data_view_2D_int);
        }
        else if (basic_type == Ioss::Field::INT64) {
          oge->put_field_data<int64_t>(field_name, pool.data_view_2D_int64);
        }
        else if (basic_type == Ioss::Field::REAL) {
          oge->put_field_data<double>(field_name, pool.data_view_2D_double);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          // Since data_view_complex cannot be a global variable.
          oge->put_field_data(out_field_name, pool.data.data(), osize);
        }
        else {
        }
        break;
      case 5:
        if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
          oge->put_field_data<char, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_char_layout_space);
        }
        else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
          oge->put_field_data<int, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_int_layout_space);
        }
        else if (basic_type == Ioss::Field::INT64) {
          oge->put_field_data<int64_t, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_int64_layout_space);
        }
        else if (basic_type == Ioss::Field::REAL) {
          oge->put_field_data<double, Kokkos::LayoutRight, Kokkos::HostSpace>(
              field_name, pool.data_view_2D_double_layout_space);
        }
        else if (basic_type == Ioss::Field::COMPLEX) {
          // Since data_view_complex cannot be a global variable.
          oge->put_field_data(out_field_name, pool.data.data(), osize);
        }
        else {
        }
        break;
#endif
      default:
        if (field_name == "mesh_model_coordinates") {
          std::cerr << "data_storage option not recognized.";
        }
        return;
      }
    }
    return;
  }

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const IOShell::Interface &interFace,
                           const std::string &prefix)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    ige->field_describe(role, &state_fields);

    // Complication here is that if the 'role' is 'Ioss::Field::MESH',
    // then the 'ids' field must be transferred first...
    if (role == Ioss::Field::MESH) {
      for (const auto &field_name : state_fields) {
        assert(oge->field_exists(field_name));
        if (field_name == "ids") {
          transfer_field_data_internal(ige, oge, field_name, interFace);
          break;
        }
      }
    }

    for (const auto &field_name : state_fields) {
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (field_name == "connectivity" && ige->type() != Ioss::ELEMENTBLOCK) {
        continue;
      }

      if (field_name != "ids" && Ioss::Utils::substr_equal(prefix, field_name)) {
        assert(oge->field_exists(field_name));
        transfer_field_data_internal(ige, oge, field_name, interFace);
      }
    }
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &       field_name,
                                    const IOShell::Interface &interFace)
  {

    size_t isize = ige->get_field(field_name).get_size();
    if (isize != oge->get_field(field_name).get_size()) {
      assert(isize == oge->get_field(field_name).get_size());
    }
    int basic_type = ige->get_field(field_name).get_type();

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
    if (field_name == "implicit_ids") {
      return;
    }
    if (field_name == "node_connectivity_status") {
      return;
    }
    if (field_name == "owning_processor") {
      return;
    }
    if (field_name == "entity_processor_raw") {
      return;
    }
    if (field_name == "ids" && ige->type() == Ioss::SIDEBLOCK) {
      return;
    }

    DataPool pool;
    pool.data.resize(isize);
    switch (interFace.data_storage_type) {
    case 1: ige->get_field_data(field_name, pool.data.data(), isize); break;
    case 2:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        ige->get_field_data(field_name, pool.data_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data(field_name, pool.data_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data(field_name, pool.data_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        ige->get_field_data(field_name, pool.data_complex);
      }
      else {
      }
      break;
#ifdef SEACAS_HAVE_KOKKOS
    case 3:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data<char>(field_name, pool.data_view_char);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        ige->get_field_data<int>(field_name, pool.data_view_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data<int64_t>(field_name, pool.data_view_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data<double>(field_name, pool.data_view_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 4:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data<char>(field_name, pool.data_view_2D_char);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        ige->get_field_data<int>(field_name, pool.data_view_2D_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data<int64_t>(field_name, pool.data_view_2D_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data<double>(field_name, pool.data_view_2D_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 5:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data<char, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_char_layout_space);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        ige->get_field_data<int, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int_layout_space);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data<int64_t, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int64_layout_space);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data<double, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_double_layout_space);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
#endif
    default:
      if (field_name == "mesh_model_coordinates") {
        std::cerr << "data_storage option not recognized.";
      }
      return;
    }

    switch (interFace.data_storage_type) {
    case 1: oge->put_field_data(field_name, pool.data.data(), isize); break;
    case 2:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        oge->put_field_data(field_name, pool.data_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data(field_name, pool.data_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data(field_name, pool.data_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        oge->put_field_data(field_name, pool.data_complex);
      }
      else {
      }
      break;
#ifdef SEACAS_HAVE_KOKKOS
    case 3:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data<char>(field_name, pool.data_view_char);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        oge->put_field_data<int>(field_name, pool.data_view_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data<int64_t>(field_name, pool.data_view_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data<double>(field_name, pool.data_view_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 4:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data<char>(field_name, pool.data_view_2D_char);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        oge->put_field_data<int>(field_name, pool.data_view_2D_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data<int64_t>(field_name, pool.data_view_2D_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data<double>(field_name, pool.data_view_2D_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 5:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data<char, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_char_layout_space);
      }
      else if ((basic_type == Ioss::Field::INTEGER) || (basic_type == Ioss::Field::INT32)) {
        oge->put_field_data<int, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int_layout_space);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data<int64_t, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int64_layout_space);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data<double, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_double_layout_space);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
#endif
    default: return;
    }
    return;
  }

  void transfer_qa_info(Ioss::Region &in, Ioss::Region &out)
  {
    out.add_information_records(in.get_information_records());

    const std::vector<std::string> &qa = in.get_qa_records();
    for (size_t i = 0; i < qa.size(); i += 4) {
      out.add_qa_record(qa[i + 0], qa[i + 1], qa[i + 2], qa[i + 3]);
    }
  }

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge)
  {
    Ioss::NameList properties;
    ige->property_describe(&properties);

    // Iterate through properties and transfer to output database...
    for (const auto &property : properties) {
      if (!oge->property_exists(property)) {
        oge->property_add(ige->get_property(property));
      }
    }
  }

  void show_step(int istep, double time)
  {
    DO_OUTPUT.setf(std::ios::scientific);
    DO_OUTPUT.setf(std::ios::showpoint);
    DO_OUTPUT << "     Time step " << std::setw(5) << istep << " at time " << std::setprecision(5)
              << time << '\n';
  }

  template <typename INT>
  void set_owned_node_count(Ioss::Region &region, int my_processor, INT /*dummy*/)
  {
    Ioss::NodeBlock *nb = region.get_node_block("nodeblock_1");
    if (nb->field_exists("owning_processor")) {
      std::vector<int> my_data;
      nb->get_field_data("owning_processor", my_data);

      INT owned = std::count(my_data.begin(), my_data.end(), my_processor);
      nb->property_add(Ioss::Property("locally_owned_count", owned));

      // Set locally_owned_count property on all nodesets...
      const Ioss::NodeSetContainer &nss = region.get_nodesets();
      for (auto ns : nss) {

        std::vector<INT> ids;
        ns->get_field_data("ids_raw", ids);
        owned = 0;
        for (size_t n = 0; n < ids.size(); n++) {
          INT id = ids[n];
          if (my_data[id - 1] == my_processor) {
            owned++;
          }
        }
        ns->property_add(Ioss::Property("locally_owned_count", owned));
      }
    }
  }
} // namespace
