// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ionit_Initializer.h"
#include "Ioss_CodeTypes.h"
#include "Ioss_DataPool.h"
#include "Ioss_FileInfo.h"
#include "Ioss_MemoryUtils.h"
#include "Ioss_MeshType.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_SerializeIO.h"
#include "Ioss_SubSystem.h"
#include "Ioss_SurfaceSplit.h"
#include "Ioss_Transform.h"
#include "Ioss_Utils.h"
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "shell_interface.h"

#ifdef SEACAS_HAVE_KOKKOS
#include <Kokkos_Core.hpp> // for Kokkos::View
#endif

#define DO_OUTPUT                                                                                  \
  if (rank == 0)                                                                                   \
  std::cerr

// ========================================================================

namespace {
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
                                    const std::string        &field_name,
                                    const IOShell::Interface &interFace);

  void file_copy(IOShell::Interface &interFace, int rank);

  Ioss::PropertyManager set_properties(IOShell::Interface &interFace);

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
  ON_BLOCK_EXIT(MPI_Finalize);
#endif
  Ioss::ParallelUtils pu{};
  rank         = pu.parallel_rank();
  int num_proc = pu.parallel_size();

#ifdef SEACAS_HAVE_KOKKOS
  Kokkos::ScopeGuard kokkos(argc, argv);
#endif

  IOShell::Interface interFace(version);
  bool               success = interFace.parse_options(argc, argv, rank);
  if (!success) {
    exit(EXIT_FAILURE);
  }

  Ioss::SerializeIO::setGroupFactor(interFace.serialize_io_size);
  mem_stats = interFace.memory_statistics;

  Ioss::Init::Initializer io;

  std::string in_file  = interFace.inputFile[0];
  std::string out_file = interFace.outputFile;

  if (rank == 0 && !interFace.quiet) {
    fmt::print(stderr,
               "Input:    '{}', Type: {}\n"
               "Output:   '{}', Type: {}\n\n",
               in_file, interFace.inFiletype, out_file, interFace.outFiletype);
  }

#ifdef SEACAS_HAVE_KOKKOS
  if (rank == 0)
    fmt::print(stderr, "Kokkos default execution space configuration:\n");
  Kokkos::DefaultExecutionSpace().print_configuration(std::cerr, false);
  if (rank == 0)
    fmt::print(stderr, "\n");
#endif

  double begin = Ioss::Utils::timer();
  file_copy(interFace, rank);
  double end = Ioss::Utils::timer();

  if (rank == 0 && !interFace.quiet) {
    if (num_proc > 1) {
      fmt::print(stderr, "\n\n\tTotal Execution Time = {:.5} seconds on {} processors.\n",
                 end - begin, num_proc);
    }
    else {
      fmt::print(stderr, "\n\n\tTotal Execution Time = {:.5} seconds.\n", end - begin);
    }
  }

  if (mem_stats) {
    int64_t MiB = 1024 * 1024;
#ifdef SEACAS_HAVE_MPI
    int64_t min, max, avg;
    int64_t hwmin, hwmax, hwavg;
    pu.memory_stats(min, max, avg);
    pu.hwm_memory_stats(hwmin, hwmax, hwavg);
    if (rank == 0) {
      fmt::print(stderr, "\n\tCurrent Memory: {}M  {}M  {}M\n", fmt::group_digits(min / MiB),
                 fmt::group_digits(max / MiB), fmt::group_digits(avg / MiB));
      fmt::print(stderr, "\tHigh Water Memory: {}M  {}M  {}M\n", fmt::group_digits(hwmin / MiB),
                 fmt::group_digits(hwmax / MiB), fmt::group_digits(hwavg / MiB));
    }
#else
    int64_t mem = Ioss::MemoryUtils::get_memory_info();
    int64_t hwm = Ioss::MemoryUtils::get_hwm_memory_info();
    if (rank == 0) {
      fmt::print(stderr,
                 "\n\tCurrent Memory:    {}M\n"
                 "\tHigh Water Memory: {}M\n",
                 fmt::group_digits(mem / MiB), fmt::group_digits(hwm / MiB));
    }
#endif
  }
  if (rank == 0) {
    fmt::print(stderr, "\n{} execution successful.\n", codename);
  }

  return EXIT_SUCCESS;
}

namespace {
  void file_copy(IOShell::Interface &interFace, int rank)
  {
    Ioss::PropertyManager properties = set_properties(interFace);

    bool first = true;
    for (const auto &inpfile : interFace.inputFile) {
      Ioss::DatabaseIO *dbi =
          Ioss::IOFactory::create(interFace.inFiletype, inpfile, Ioss::READ_MODEL,
                                  Ioss::ParallelUtils::comm_world(), properties);
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
          if (rank == 0) {
            fmt::print(stderr, "ERROR: Unable to open group '{}' in file '{}'\n",
                       interFace.groupName, inpfile);
          }
          return;
        }
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region region(dbi, "region_1");

      if (region.mesh_type() != Ioss::MeshType::UNSTRUCTURED) {
        if (rank == 0) {
          fmt::print(stderr,
                     "\nERROR: io_shell does not support '{}' meshes. Only 'Unstructured' mesh is "
                     "supported at this time.\n",
                     region.mesh_type_string());
        }
        return;
      }

      // Get length of longest name on input file...
      int max_name_length = dbi->maximum_symbol_length();
      if (max_name_length > 0) {
        properties.add(Ioss::Property("MAXIMUM_NAME_LENGTH", max_name_length));
      }

      // Get integer size being used on the input file and propagate
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
                                  Ioss::ParallelUtils::comm_world(), properties);
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

      // This also assumes that the node order and count is the same for input
      // and output regions... (This is checked during nodeset output)
      if (output_region.get_database()->needs_shared_node_information()) {
        if (interFace.ints_64_bit)
          set_owned_node_count(region, rank, (int64_t)0);
        else
          set_owned_node_count(region, rank, (int)0);
      }

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

      int step_count = region.get_optional_property("state_count", 0);
      if (step_count > 0) {
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

      auto *nb = new Ioss::NodeBlock(*inb);
      output_region.add(nb);

      if (output_region.get_database()->needs_shared_node_information()) {
        // If the "owning_processor" field exists on the input
        // nodeblock, transfer it and the "ids" field to the output
        // nodeblock at this time since it is used to determine
        // per-processor sizes of nodeblocks and nodesets.
        if (inb->field_exists("owning_processor")) {
          size_t            isize = inb->get_field("ids").get_size();
          std::vector<char> data(isize);
          inb->get_field_data("ids", Data(data), isize);
          nb->put_field_data("ids", Data(data), isize);
          isize = inb->get_field("owning_processor").get_size();
          data.resize(isize);
          inb->get_field_data("owning_processor", Data(data), isize);
          nb->put_field_data("owning_processor", Data(data), isize);
        }
      }
    }
    if (debug) {
      DO_OUTPUT << '\n';
    }
  }

  struct param
  {
    Ioss::GroupingEntity     *entity;
    Ioss::Region             *output_region;
    Ioss::Field::RoleType     role;
    const IOShell::Interface *interFace;
  };

  void *transfer_fields_ts(void *varg)
  {
    auto *arg           = static_cast<param *>(varg);
    auto  entity        = arg->entity;
    auto  output_region = arg->output_region;
    auto  interFace     = arg->interFace;
    auto  role          = arg->role;

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
      pthread_create(&threads[t], nullptr, transfer_fields_ts, (void *)(Data(params) + t));
      t++;
    }

    for (t = 0; t < (int)threads.size(); t++) {
      pthread_join(threads[t], nullptr);
    }
  }

  void *transfer_field_data_ts(void *varg)
  {
    auto *arg           = static_cast<param *>(varg);
    auto  entity        = arg->entity;
    auto  output_region = arg->output_region;
    auto  interFace     = arg->interFace;
    auto  role          = arg->role;

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
      pthread_create(&threads[t], nullptr, transfer_field_data_ts, (void *)(Data(params) + t));
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
        size_t count = iblock->entity_count();
        total_entities += count;

        auto *block = new T(*iblock);
        output_region.add(block);
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
    const auto &fss = region.get_sidesets();
    for (const auto &ss : fss) {
      const std::string &name = ss->name();
      if (debug) {
        DO_OUTPUT << name << ", ";
      }

      auto *surf = new Ioss::SideSet(*ss);
      output_region.add(surf);

      // Fix up the optional 'owner_block' in copied SideBlocks...
      const auto &fbs = ss->get_side_blocks();
      for (const auto &ifb : fbs) {
        if (ifb->parent_block() != nullptr) {
          const auto &fb_name = ifb->parent_block()->name();
          auto       *parent  = dynamic_cast<Ioss::EntityBlock *>(
              output_region.get_entity(fb_name, Ioss::ELEMENTBLOCK));
          if (parent == nullptr) {
            parent = dynamic_cast<Ioss::EntityBlock *>(
                output_region.get_entity(fb_name, Ioss::STRUCTUREDBLOCK));
          }

          auto *ofb = surf->get_side_block(ifb->name());
          ofb->set_parent_block(parent);
        }
      }
    }
    if (!debug) {
      DO_OUTPUT << " Number of        SideSets            =" << std::setw(12) << fss.size() << "\n";
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
        auto *o_set = new T(*set);
        output_region.add(o_set);
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
      if (debug) {
        const std::string &name = ics->name();
        DO_OUTPUT << name << ", ";
      }
      auto *cs = new Ioss::CommSet(*ics);
      output_region.add(cs);
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
    Ioss::NameList fields = ige->field_describe(role);

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
    Ioss::NameList fields = ige->field_describe(role);

    // Iterate through results fields and transfer to output database...
    for (const auto &field_name : fields) {
      std::string out_field_name = field_name + "_mag";
      if (!oge->field_exists(out_field_name)) {
        // If the field does not already exist, add it to the output node block
        Ioss::Field field = ige->get_field(field_name);
        Ioss::Field tr_field(out_field_name, field.get_type(), field.raw_storage(),
                             field.get_role(), field.raw_count());

        auto *transform = Ioss::Transform::create("vector magnitude");
        assert(transform != nullptr);
        tr_field.add_transform(transform);

        auto *max_transform = Ioss::Transform::create("absolute_maximum");
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
    Ioss::NameList state_fields = ige->field_describe(role);
    // Iterate through mesh description fields and transfer to
    // output database...
    for (const auto &field_name : state_fields) {
      std::string out_field_name = field_name + "_mag";

      assert(oge->field_exists(out_field_name));

      int basic_type = ige->get_field(field_name).get_type();

      size_t isize = ige->get_field(field_name).get_size();
      size_t osize = oge->get_field(out_field_name).get_size();
      assert(isize == osize);

      Ioss::DataPool pool;
      pool.data.resize(isize);
      switch (interFace.data_storage_type) {
      case 1: ige->get_field_data(field_name, Data(pool.data), isize); break;
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
          ige->get_field_data(field_name, Data(pool.data), isize);
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
          ige->get_field_data(field_name, Data(pool.data), isize);
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
          ige->get_field_data(field_name, Data(pool.data), isize);
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
      case 1: oge->put_field_data(out_field_name, Data(pool.data), osize); break;
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
          oge->put_field_data(out_field_name, Data(pool.data), osize);
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
          oge->put_field_data(out_field_name, Data(pool.data), osize);
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
          oge->put_field_data(out_field_name, Data(pool.data), osize);
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
    Ioss::NameList state_fields = ige->field_describe(role);

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
                                    const std::string        &field_name,
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

    Ioss::DataPool pool;
    pool.data.resize(isize);
    switch (interFace.data_storage_type) {
    case 1: ige->get_field_data(field_name, Data(pool.data), isize); break;
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
        ige->get_field_data(field_name, Data(pool.data), isize);
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
        ige->get_field_data(field_name, Data(pool.data), isize);
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
        ige->get_field_data(field_name, Data(pool.data), isize);
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
    case 1: oge->put_field_data(field_name, Data(pool.data), isize); break;
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
        oge->put_field_data(field_name, Data(pool.data), isize);
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
        oge->put_field_data(field_name, Data(pool.data), isize);
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
        oge->put_field_data(field_name, Data(pool.data), isize);
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

    const Ioss::NameList &qa = in.get_qa_records();
    for (size_t i = 0; i < qa.size(); i += 4) {
      out.add_qa_record(qa[i + 0], qa[i + 1], qa[i + 2], qa[i + 3]);
    }
  }

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge)
  {
    Ioss::NameList properties = ige->property_describe();

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

  Ioss::PropertyManager set_properties(IOShell::Interface &interFace)
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

    if (interFace.compression_level > 0 || interFace.shuffle || interFace.szip) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", static_cast<int>(interFace.shuffle)));

      if (interFace.szip) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "szip"));
      }
      else if (interFace.zlib) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "zlib"));
      }
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

    if (interFace.file_per_state) {
      properties.add(Ioss::Property("FILE_PER_STATE", "YES"));
    }

    if (interFace.netcdf4) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    if (interFace.netcdf5) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf5"));
    }

    if (interFace.inputFile.size() > 1) {
      properties.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
    }

    if (interFace.debug) {
      properties.add(Ioss::Property("LOGGING", 1));
    }

    if (interFace.memory_statistics) {
      properties.add(Ioss::Property("ENABLE_TRACING", 1));
    }

    if (!interFace.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interFace.decomp_method));
    }

    if (interFace.retain_empty_blocks) {
      properties.add(Ioss::Property("RETAIN_EMPTY_BLOCKS", "YES"));
    }
    return properties;
  }
} // namespace
