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
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Ionit_Initializer.h"
#include "Ioss_CodeTypes.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_FaceGenerator.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_Utils.h"

#include <cassert>

#include "skinner_interface.h"

#define OUTPUT std::cout

// ========================================================================

namespace {
  template <typename INT> void skinner(Skinner::Interface &interface, INT /*dummy*/);
  std::string codename;
  std::string version = "0.6";
} // namespace

int main(int argc, char *argv[])
{
  int my_rank = 0;
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  Skinner::Interface interface;
  interface.parse_options(argc, argv);

  std::string in_type = "exodusII";

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  Ioss::Init::Initializer io;

  if (my_rank == 0) {
    OUTPUT << "Input:    '" << interface.input_filename() << "', Type: " << interface.input_type()
           << '\n';
    if (!interface.no_output()) {
      OUTPUT << "Output:   '" << interface.output_filename()
             << "', Type: " << interface.output_type() << '\n';
    }
  }

  if (interface.ints_64_bit()) {
    skinner(interface, static_cast<int64_t>(0));
  }
  else {
    skinner(interface, 0);
  }

  if (my_rank == 0) {
    OUTPUT << "\n" << codename << " execution successful.\n";
  }
#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
  template <typename INT> void skinner(Skinner::Interface &interface, INT /*dummy*/)
  {
    std::string inpfile    = interface.input_filename();
    std::string input_type = interface.input_type();

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

    Ioss::FaceGenerator face_generator(region);
#ifdef SEACAS_HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    auto start = std::chrono::steady_clock::now();
    face_generator.generate_faces((INT)0);
#ifdef SEACAS_HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    auto duration = std::chrono::steady_clock::now() - start;

    Ioss::FaceUnorderedSet &faces = face_generator.faces();

    // Faces have been generated at this point.
    // Categorize (boundary/interior)
    size_t interior  = 0;
    size_t boundary  = 0;
    size_t error     = 0;
    size_t pboundary = 0;

    for (auto &face : faces) {
      if (face.elementCount_ == 2) {
        interior++;
        if (face.sharedWithProc_ != -1) {
          pboundary++;
        }
      }
      else if (face.elementCount_ == 1) {
        boundary++;
      }
      else {
        error++;
      }
    }

#ifdef SEACAS_HAVE_MPI
    Ioss::Int64Vector counts(3), global(3);
    counts[0] = interior;
    counts[1] = boundary;
    counts[2] = pboundary;
    region.get_database()->util().global_count(counts, global);
    interior  = global[0];
    boundary  = global[1];
    pboundary = global[2];
#endif

    size_t my_rank = region.get_database()->parallel_rank();
    if (my_rank == 0) {
      OUTPUT << "Face count = " << interior + boundary - pboundary / 2
             << "\tInterior = " << interior - pboundary / 2 << "\tBoundary = " << boundary
             << "\tShared   = " << pboundary << "\tError = " << error << "\n"
             << "Total Time = " << std::chrono::duration<double, std::milli>(duration).count()
             << " ms\t"
             << (interior + boundary - pboundary / 2) /
                    std::chrono::duration<double>(duration).count()
             << " faces/second\n\n";

      OUTPUT << "Hash Statistics: Bucket Count = " << faces.bucket_count()
             << "\tLoad Factor = " << faces.load_factor() << "\n";
      size_t numel = region.get_property("element_count").get_int();
      OUTPUT << "Faces/Element ratio = " << static_cast<double>(faces.size()) / numel << "\n";
    }

    if (interface.no_output()) {
      return;
    }

    // Get vector of all boundary faces which will be output as the skin...
    std::vector<Ioss::Face> boundary_faces;
    boundary_faces.reserve(boundary);
    for (auto &face : faces) {
      if (face.elementCount_ == 1) {
        boundary_faces.push_back(face);
      }
    }

    // Iterate the boundary faces and determine which nodes are referenced...
    size_t           node_count = region.get_property("node_count").get_int();
    std::vector<int> ref_nodes(node_count);
    size_t           quad = 0;
    size_t           tri  = 0;
    for (const auto &face : boundary_faces) {
      size_t face_node_count = 0;
      for (auto &gnode : face.connectivity_) {
        if (gnode > 0) {
          face_node_count++;
          auto node       = region.get_database()->node_global_to_local(gnode, true) - 1;
          ref_nodes[node] = 1;
        }
      }
      if (face_node_count == 3) {
        tri++;
      }
      else if (face_node_count == 4) {
        quad++;
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

    std::vector<INT>    ref_ids(ref_count);
    std::vector<double> coord_out(3 * ref_count);

    size_t j = 0;
    for (size_t i = 0; i < ref_nodes.size(); i++) {
      if (ref_nodes[i] == 1) {
        coord_out[3 * j + 0] = coord_in[3 * i + 0];
        coord_out[3 * j + 1] = coord_in[3 * i + 1];
        coord_out[3 * j + 2] = coord_in[3 * i + 2];
        ref_ids[j++]         = ids[i];
      }
    }

    // Create output file...
    if (interface.compression_level > 0 || interface.shuffle) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interface.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", interface.shuffle));
    }

    if (interface.compose_output != "none") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
      if (interface.compose_output != "default") {
        properties.add(Ioss::Property("PARALLEL_IO_MODE", interface.compose_output));
      }
    }

    if (interface.netcdf4) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    std::string       file = interface.output_filename();
    std::string       type = interface.output_type();
    Ioss::DatabaseIO *dbo  = Ioss::IOFactory::create(type, file, Ioss::WRITE_RESTART,
                                                    (MPI_Comm)MPI_COMM_WORLD, properties);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "skin");
    output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    Ioss::NodeBlock *nbo =
        new Ioss::NodeBlock(output_region.get_database(), "nodeblock_1", ref_count, 3);
    output_region.add(nbo);

    Ioss::ElementBlock *quadblock = nullptr;
    Ioss::ElementBlock *triblock  = nullptr;
    if (quad > 0) {
      quadblock = new Ioss::ElementBlock(output_region.get_database(), "quad", "shell", quad);
      output_region.add(quadblock);
    }

    if (tri > 0) {
      triblock = new Ioss::ElementBlock(output_region.get_database(), "tri", "trishell", tri);
      output_region.add(triblock);
    }

    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    output_region.begin_mode(Ioss::STATE_MODEL);
    nbo->put_field_data("ids", ref_ids);
    nbo->put_field_data("mesh_model_coordinates", coord_out);

    std::vector<INT> quad_conn;
    std::vector<INT> quad_ids;
    std::vector<INT> tri_conn;
    std::vector<INT> tri_ids;
    quad_conn.reserve(4 * quad);
    quad_ids.reserve(quad);
    tri_conn.reserve(3 * tri);
    tri_ids.reserve(tri);

    bool use_face_ids = !interface.ignoreFaceIds_;
    INT  fid          = 1;
    for (auto &face : boundary_faces) {
      if (use_face_ids) {
        fid = face.id_;
        if (fid < 0) {
          fid = -fid;
        }
      }
      else {
        fid++;
      }

      if (face.connectivity_[3] != 0) {
        for (int i = 0; i < 4; i++) {
          quad_conn.push_back(face.connectivity_[i]);
        }
        quad_ids.push_back(fid);
      }
      else {
        for (int i = 0; i < 3; i++) {
          tri_conn.push_back(face.connectivity_[i]);
        }
        tri_ids.push_back(fid);
      }
    }

    if (quad > 0) {
      quadblock->put_field_data("ids", quad_ids);
      quadblock->put_field_data("connectivity", quad_conn);
    }
    if (tri > 0) {
      triblock->put_field_data("ids", tri_ids);
      triblock->put_field_data("connectivity", tri_conn);
    }

    output_region.end_mode(Ioss::STATE_MODEL);
  }
} // namespace
