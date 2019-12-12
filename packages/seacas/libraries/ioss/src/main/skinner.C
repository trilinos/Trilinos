// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <functional>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "Ionit_Initializer.h"
#include "Ioss_CodeTypes.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_FaceGenerator.h"
#include "Ioss_FileInfo.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_Utils.h"

#include <cassert>

#include "skinner_interface.h"

// ========================================================================

namespace {
  template <typename INT> void skinner(Skinner::Interface &interface, INT /*dummy*/);
  std::string                  codename;
  std::string                  version = "0.99";
} // namespace

int main(int argc, char *argv[])
{
  int my_rank = 0;
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  codename = Ioss::FileInfo(argv[0]).basename();

  Skinner::Interface interface;
  bool               success = interface.parse_options(argc, argv);
  if (!success) {
    return EXIT_FAILURE;
  }

  Ioss::Init::Initializer io;

  if (my_rank == 0) {
    fmt::print("\nInput:    '{}', Type: {}\n", interface.input_filename(), interface.input_type());
    if (!interface.no_output()) {
      fmt::print("Output:   '{}', Type: {}\n", interface.output_filename(),
                 interface.output_type());
    }
  }

  double begin = Ioss::Utils::timer();
  try {
    if (interface.ints_64_bit()) {
      skinner(interface, static_cast<int64_t>(0));
    }
    else {
      skinner(interface, 0);
    }
  }
  catch (std::exception &e) {
    fmt::print(stderr, "\n{}\n\nskinner terminated due to exception\n", e.what());
    exit(EXIT_FAILURE);
  }
#ifdef SEACAS_HAVE_MPI
  Ioss::ParallelUtils parallel(MPI_COMM_WORLD);
  parallel.barrier();
#endif
  double end = Ioss::Utils::timer();

  if (my_rank == 0) {
    fmt::print("\n\tTotal Execution time = {} seconds\n", end - begin);
    fmt::print("\n{} execution successful.\n\n", codename);
  }
  return EXIT_SUCCESS;
}

namespace {
  Ioss::PropertyManager set_properties(Skinner::Interface &interface);

  template <typename INT> void skinner(Skinner::Interface &interface, INT /*dummy*/)
  {
    std::string inpfile    = interface.input_filename();
    std::string input_type = interface.input_type();

    Ioss::PropertyManager properties = set_properties(interface);

    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
                                                    (MPI_Comm)MPI_COMM_WORLD, properties);
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    if (interface.ints_64_bit()) {
      dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");
    int my_rank = region.get_database()->util().parallel_rank();

    if (my_rank == 0) {
      region.output_summary(std::cerr, false);
    }

    // Generate the faces...
    Ioss::FaceGenerator face_generator(region);
    face_generator.generate_faces((INT)0, true);

    if (interface.no_output()) {
      return;
    }

    // Get vector of all boundary faces which will be output as the skin...
    std::map<std::string, std::vector<Ioss::Face>> boundary_faces;
    const Ioss::ElementBlockContainer &            ebs = region.get_element_blocks();
    for (auto eb : ebs) {
      const std::string &name     = eb->name();
      auto &             boundary = boundary_faces[name];
      auto &             faces    = face_generator.faces(name);
      for (auto &face : faces) {
        if (face.elementCount_ == 1) {
          boundary.push_back(face);
        }
      }
    }

    // Iterate the boundary faces and determine which nodes are referenced...
    size_t           my_node_count = region.get_property("node_count").get_int();
    std::vector<int> ref_nodes(my_node_count);
    for (const auto &boundaries : boundary_faces) {
      for (const auto &face : boundaries.second) {
        for (auto &gnode : face.connectivity_) {
          if (gnode > 0) {
            auto node       = region.get_database()->node_global_to_local(gnode, true) - 1;
            ref_nodes[node] = 1;
          }
        }
      }
    }

    size_t ref_count = std::accumulate(ref_nodes.begin(), ref_nodes.end(), 0);

    // Map ids from total mesh down to skin mesh...
    // Also get coordinates...
    std::vector<double> coord_in;
    Ioss::NodeBlock *   nb = region.get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", coord_in);

    std::vector<INT> ids;
    nb->get_field_data("ids", ids);


    std::vector<int> owner;
    nb->get_field_data("owning_processor", owner);

    std::vector<INT>    ref_ids(ref_count);
    std::vector<double> coord_out(3 * ref_count);
    std::vector<int>    owner_out(ref_count);

    size_t j = 0;
    for (size_t i = 0; i < ref_nodes.size(); i++) {
      if (ref_nodes[i] == 1) {
        coord_out[3 * j + 0] = coord_in[3 * i + 0];
        coord_out[3 * j + 1] = coord_in[3 * i + 1];
        coord_out[3 * j + 2] = coord_in[3 * i + 2];
        owner_out[j]         = owner[i];
        ref_ids[j++]         = ids[i];
      }
    }
    Ioss::Utils::clear(coord_in);
    Ioss::Utils::clear(ids);
    Ioss::Utils::clear(owner);

    std::string       file = interface.output_filename();
    std::string       type = interface.output_type();
    Ioss::DatabaseIO *dbo  = Ioss::IOFactory::create(type, file, Ioss::WRITE_RESTART,
                                                    (MPI_Comm)MPI_COMM_WORLD, properties);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "skin");
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    Ioss::NodeBlock *nbo =
        new Ioss::NodeBlock(output_region.get_database(), "nodeblock_1", ref_count, 3);

    // Count number of nodes owned by this processor (owner_out[i] == myProcessor);
    size_t owned = std::count_if(owner_out.begin(), owner_out.end(),
                                 [my_rank](int i) { return i == my_rank; });
    nbo->property_add(Ioss::Property("locally_owned_count", (INT)owned));

    output_region.add(nbo);

    // Output element blocks: will have a "skin" block for each input element block.
    //
    // NOTE: currently only handle a single skin face topology per
    // input element block. Eventually may need to handle the case
    // where there are two skin face topologies per input element
    // block (wedge -> tri and quad).

    // Count faces per element block and create output element block...
    for (auto eb : ebs) {
      const std::string &name      = eb->name();
      auto &             boundary  = boundary_faces[name];
      auto               face_topo = eb->topology()->face_type(0);
      std::string        topo      = "shell";
      if (face_topo == nullptr) {
	std::ostringstream errmsg;
	fmt::print(errmsg, "ERROR: Block '{}' with topology '{}' does not have"
		   " a unique face topology.\nThis is not supported at this time.\n",
		   name, eb->topology()->name());
	IOSS_ERROR(errmsg);
      }
      if (face_topo->name() == "tri3") {
        topo = "trishell";
      }
      auto block =
          new Ioss::ElementBlock(output_region.get_database(), name, topo, boundary.size());
      output_region.add(block);
    }
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    output_region.begin_mode(Ioss::STATE_MODEL);
    nbo->put_field_data("ids", ref_ids);
    nbo->put_field_data("owning_processor", owner_out);
    nbo->put_field_data("mesh_model_coordinates", coord_out);
    Ioss::Utils::clear(coord_out);
    Ioss::Utils::clear(owner_out);
    Ioss::Utils::clear(ref_ids);

    bool use_face_hash_ids = interface.useFaceHashIds_;
    INT  fid               = 0;
    for (auto eb : ebs) {
      const std::string &name       = eb->name();
      auto &             boundary   = boundary_faces[name];
      auto *             block      = output_region.get_element_block(name);
      size_t             node_count = block->topology()->number_corner_nodes();

      std::vector<INT> conn;
      std::vector<INT> elids;
      conn.reserve(node_count * boundary.size());
      elids.reserve(boundary.size());

      for (auto &face : boundary) {
        if (use_face_hash_ids) {
          fid = face.hashId_;
          if (fid < 0) { // Due to (unsigned)size_t -> (signed)INT conversion
            fid = -fid;
          }
        }
        else {
          fid = face.element[0];
        }

        for (size_t i = 0; i < node_count; i++) {
          conn.push_back(face.connectivity_[i]);
        }
        elids.push_back(fid);
      }
      block->put_field_data("ids", elids);
      block->put_field_data("connectivity", conn);
    }
    output_region.end_mode(Ioss::STATE_MODEL);
    if (my_rank == 0) {
      output_region.output_summary(std::cerr, false);
    }
  }

  Ioss::PropertyManager set_properties(Skinner::Interface &interface)
  {
    Ioss::PropertyManager properties;
    if (interface.ints_64_bit()) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interface.debug) {
      properties.add(Ioss::Property("LOGGING", 1));
    }

    if (!interface.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interface.decomp_method));
    }

    if (interface.compression_level > 0 || interface.shuffle) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interface.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", interface.shuffle));
    }

    if (interface.compose_output == "default") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
    }
    else if (interface.compose_output == "external") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
    }
    else if (interface.compose_output != "none") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
    }

    if (interface.netcdf4_) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    return properties;
  }
} // namespace
