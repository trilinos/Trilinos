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

#include <Ioss_CodeTypes.h>
#include <Ioss_CommSet.h>
#include <Ioss_DBUsage.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_FaceGenerator.h>
#include <Ioss_IOFactory.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Property.h>
#include <Ioss_Region.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <random>
#include <utility>

// Options for generating hash function key...
#define USE_MURMUR
//#define USE_RANDOM

namespace {
  template <typename T> void generate_index(std::vector<T> &index)
  {
    T sum = 0;
    for (auto &idx : index) {
      T cnt = idx;
      idx   = sum;
      sum += cnt;
    }
  }

#if defined(USE_MURMUR)
  uint64_t MurmurHash64A(const void *key, int len, uint64_t seed);
#endif

  size_t id_hash(size_t id)
  {
#if defined(USE_RANDOM)
    std::mt19937_64 rng;
    rng.seed(id);
    return rng();
#elif defined(USE_MURMUR)
    return MurmurHash64A(&id, sizeof(size_t), 24713);
#else
    return id;
#endif
  }

  void create_face(Ioss::FaceUnorderedSet &faces, size_t id, std::array<size_t, 4> &conn,
                   size_t element)
  {
    Ioss::Face face(id, conn);
    auto       face_iter = faces.insert(face);

    (*(face_iter.first)).add_element(element);
  }

  template <typename INT>
  void resolve_parallel_faces(Ioss::Region &region, Ioss::FaceUnorderedSet &faces,
                              const std::vector<size_t> &hash_ids, INT /*dummy*/)
  {
#ifdef SEACAS_HAVE_MPI
    size_t proc_count = region.get_database()->util().parallel_size();

    if (proc_count > 1) {
      // If parallel, resolve faces on processor boundaries.
      // For each boundary face, need to check whether all of the nodes on
      // the face are shared with the same processor.  If so, then that face
      // is *possibly* shared with that processor.
      //
      // With the current continuum element only restriction, then a face
      // can only be shared with one other processor...

      // get nodal communication data CommSet...
      Ioss::CommSet *css = region.get_commset("commset_node");

      std::vector<std::pair<INT, INT>> proc_entity;
      {
        // entity_processor consists of node,proc, node,proc, entries.
        std::vector<INT> entity_processor;
        css->get_field_data("entity_processor_raw", entity_processor);

        proc_entity.reserve(entity_processor.size() / 2);
        for (size_t i = 0; i < entity_processor.size(); i += 2) {
          // Converts from 1-based to 0-based local nodes.
          proc_entity.push_back(std::make_pair(entity_processor[i + 1], entity_processor[i] - 1));
        }
      }

      // 'id_span' gives index into proc_entity for all nodes.
      // 'id_span[local_node_id] .. id_span[local_node_id+1]' gives
      // the location in 'proc_entity' of the sharing information
      // for node 'local_node_id'
      std::vector<size_t> id_span(hash_ids.size() + 1);
      for (const auto &pe : proc_entity) {
        INT node = pe.second;
        assert(node >= 0 && node < (INT)id_span.size() - 1);
        id_span[node]++;
      }
      generate_index(id_span);

      // Each boundary face ...
      // .. See if all of its nodes are shared with same processor.
      //  .. Iterate face nodes
      //  .. Determine shared proc.
      //  .. (for now, use a map of <proc,count>
      //   .. if potentially shared with 'proc', then count == num_nodes_face

      std::vector<INT> potential_count(proc_count);
      for (auto &face : faces) {
        if (face.elementCount_ == 1) {
          // On 'boundary' -- try to determine whether on processor or exterior
          // boundary
          std::map<int, int> shared_nodes;
          int                face_node_count = 0;
          for (auto &gnode : face.connectivity_) {
            if (gnode > 0) {
              auto node = region.get_database()->node_global_to_local(gnode, true) - 1;
              face_node_count++;
              size_t begin = id_span[node];
              size_t end   = id_span[node + 1];
              for (size_t j = begin; j < end; j++) {
                assert(proc_entity[j].second == node);
                int proc = proc_entity[j].first;
                shared_nodes[proc]++;
              }
            }
          }
          for (auto &node : shared_nodes) {
            if (node.second == face_node_count) {
              potential_count[node.first]++;
            }
          }
        }
      }

      std::vector<INT> potential_offset(potential_count.begin(), potential_count.end());
      generate_index(potential_offset);

      size_t potential = potential_offset[proc_count - 1] + potential_count[proc_count - 1];
      std::vector<int64_t> potential_faces(6 * potential);

      for (auto &face : faces) {
        if (face.elementCount_ == 1) {
          // On 'boundary' -- try to determine whether on processor or exterior
          // boundary
          std::map<int, int> shared_nodes;
          int                face_node_count = 0;
          for (auto &gnode : face.connectivity_) {
            if (gnode > 0) {
              auto node = region.get_database()->node_global_to_local(gnode, true) - 1;
              face_node_count++;
              size_t begin = id_span[node];
              size_t end   = id_span[node + 1];
              for (size_t j = begin; j < end; j++) {
                assert(proc_entity[j].second == node);
                int proc = proc_entity[j].first;
                shared_nodes[proc]++;
              }
            }
          }
          for (auto &node : shared_nodes) {
            if (node.second == face_node_count) {
              size_t offset                   = potential_offset[node.first];
              potential_faces[6 * offset + 0] = face.id_;
              potential_faces[6 * offset + 1] = face.connectivity_[0];
              potential_faces[6 * offset + 2] = face.connectivity_[1];
              potential_faces[6 * offset + 3] = face.connectivity_[2];
              potential_faces[6 * offset + 4] = face.connectivity_[3];
              potential_faces[6 * offset + 5] = face.element[0];
              assert(face.elementCount_ == 1);
              potential_offset[node.first]++;
            }
          }
        }
      }

      // Regenerate potential_offset since it was modified above...
      std::copy(potential_count.begin(), potential_count.end(), potential_offset.begin());
      generate_index(potential_offset);

      // Now need to send to the other processors...
      // For now, use all-to-all; optimization is just send to processors with
      // data...
      std::vector<INT> check_count(proc_count);
      MPI_Alltoall(TOPTR(potential_count), 1, Ioss::mpi_type((INT)0), TOPTR(check_count), 1,
                   Ioss::mpi_type((INT)0), region.get_database()->util().communicator());

      const int            values_per_face = 6; // id, 4-node-conn, element
      auto                 sum = std::accumulate(check_count.begin(), check_count.end(), 0);
      std::vector<int64_t> check_faces(values_per_face * sum);

      std::vector<INT> check_offset(check_count.begin(), check_count.end());
      generate_index(check_offset);

      // Need to adjust counts and offsets to account for sending 6 values per
      // face...
      for (size_t i = 0; i < proc_count; i++) {
        potential_count[i] *= values_per_face;
        potential_offset[i] *= values_per_face;
        check_count[i] *= values_per_face;
        check_offset[i] *= values_per_face;
      }

      Ioss::MY_Alltoallv(potential_faces, potential_count, potential_offset, check_faces,
                         check_count, check_offset, region.get_database()->util().communicator());

      // Now iterate the check_faces and see if any of them match one
      // of this processors faces...  If so, then mark as shared and
      // add the element...
      for (size_t i = 0; i < check_faces.size(); i += values_per_face) {
        size_t                id = check_faces[i + 0];
        std::array<size_t, 4> conn{};
        conn[0]            = check_faces[i + 1];
        conn[1]            = check_faces[i + 2];
        conn[2]            = check_faces[i + 3];
        conn[3]            = check_faces[i + 4];
        size_t     element = check_faces[i + 5];
        Ioss::Face face(id, conn);
        auto       face_iter = faces.find(face);
        if (face_iter != faces.end()) {
          // we have a match... This is a shared interior face
          (*face_iter).add_element(element);

          int proc = 0;
          for (size_t j = 0; j < check_count.size(); j++) {
            if (check_count[j] > 0 && check_offset[j] == (INT)i) {
              break;
            }
            proc++;
          }
          (*face_iter).sharedWithProc_ = proc;
        }
      }
    }
#endif
  }
} // namespace

namespace Ioss {
  FaceGenerator::FaceGenerator(Ioss::Region &region) : region_(region) {}

  template void FaceGenerator::generate_faces(int);
  template void FaceGenerator::generate_faces(int64_t);

  template <typename INT> void FaceGenerator::generate_faces(INT /*dummy*/)
  {
    Ioss::NodeBlock *nb = region_.get_node_blocks()[0];

    std::vector<INT> ids;
    nb->get_field_data("ids", ids);

    // Convert ids into hashed-ids
    auto                starth = std::chrono::high_resolution_clock::now();
    std::vector<size_t> hash_ids;
    hash_ids.reserve(ids.size());
    for (auto id : ids) {
      hash_ids.push_back(id_hash(id));
    }
    auto endh = std::chrono::high_resolution_clock::now();

    size_t numel = region_.get_property("element_count").get_int();

    faces_.reserve(3.3 * numel);

    const Ioss::ElementBlockContainer &ebs = region_.get_element_blocks();
    for (auto eb : ebs) {
      const Ioss::ElementTopology *topo = eb->topology();

      // Only handle continuum elements at this time...
      if (topo->parametric_dimension() != 3) {
        continue;
      }

      std::vector<INT> connectivity;
      eb->get_field_data("connectivity_raw", connectivity);

      std::vector<INT> elem_ids;
      eb->get_field_data("ids", elem_ids);

      int num_face_per_elem = topo->number_faces();
      assert(num_face_per_elem <= 6);
      std::array<Ioss::IntVector, 6> face_conn;
      std::array<int, 6>             face_count{};
      for (int face = 0; face < num_face_per_elem; face++) {
        face_conn[face]  = topo->face_connectivity(face + 1);
        face_count[face] = topo->face_type(face + 1)->number_corner_nodes();
      }

      int    num_node_per_elem = topo->number_nodes();
      size_t num_elem          = eb->entity_count();

      for (size_t elem = 0, offset = 0; elem < num_elem; elem++, offset += num_node_per_elem) {
        for (int face = 0; face < num_face_per_elem; face++) {
          size_t id = 0;
          assert(face_count[face] <= 4);
          std::array<size_t, 4> conn = {{0, 0, 0, 0}};
          for (int j = 0; j < face_count[face]; j++) {
            size_t fnode = offset + face_conn[face][j];
            size_t gnode = connectivity[fnode];
            conn[j]      = ids[gnode - 1];
            id += hash_ids[gnode - 1];
          }
          create_face(faces_, id, conn, elem_ids[elem]);
        }
      }
    }

    auto endf = std::chrono::high_resolution_clock::now();
    resolve_parallel_faces(region_, faces_, hash_ids, (INT)0);
    auto endp = std::chrono::high_resolution_clock::now();

    auto diffh = endh - starth;
    auto difff = endf - endh;
    std::cout << "Node ID hash time:   \t"
              << std::chrono::duration<double, std::milli>(diffh).count() << " ms\t"
              << hash_ids.size() / std::chrono::duration<double>(diffh).count()
              << " nodes/second\n";
    std::cout << "Face generation time:\t"
              << std::chrono::duration<double, std::milli>(difff).count() << " ms\t"
              << faces_.size() / std::chrono::duration<double>(difff).count() << " faces/second.\n";
#ifdef SEACAS_HAVE_MPI
    auto   diffp      = endp - endf;
    size_t proc_count = region_.get_database()->util().parallel_size();

    if (proc_count > 1) {
      std::cout << "Parallel time:       \t"
                << std::chrono::duration<double, std::milli>(diffp).count() << " ms\t"
                << faces_.size() / std::chrono::duration<double>(diffp).count()
                << " faces/second.\n";
    }
#endif
    std::cout << "Total time:          \t"
              << std::chrono::duration<double, std::milli>(endp - starth).count() << " ms\n\n";
  }
} // namespace Ioss

namespace {
#if defined(USE_MURMUR)
//-----------------------------------------------------------------------------
// MurmurHash2 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

// Note - This code makes a few assumptions about how your machine behaves -

// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4

// And it has a few limitations -

// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.

//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby

// The same caveats as 32-bit MurmurHash2 apply here - beware of alignment
// and endian-ness issues if used across multiple platforms.

// 64-bit hash for 64-bit platforms
#define BIG_CONSTANT(x) (x##LLU)
  uint64_t MurmurHash64A(const void *key, int len, uint64_t seed)
  {
    const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
    const int      r = 47;

    uint64_t h = seed ^ (len * m);

    const uint64_t *data = reinterpret_cast<const uint64_t *>(key);
    const uint64_t *end  = data + (len / 8);

    while (data != end) {
      uint64_t k = *data++;

      k *= m;
      k ^= k >> r;
      k *= m;

      h ^= k;
      h *= m;
    }

    const unsigned char *data2 = reinterpret_cast<const unsigned char *>(data);

    switch (len & 7) {
    case 7: h ^= uint64_t(data2[6]) << 48;
    case 6: h ^= uint64_t(data2[5]) << 40;
    case 5: h ^= uint64_t(data2[4]) << 32;
    case 4: h ^= uint64_t(data2[3]) << 24;
    case 3: h ^= uint64_t(data2[2]) << 16;
    case 2: h ^= uint64_t(data2[1]) << 8;
    case 1: h ^= uint64_t(data2[0]); h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
  }
#endif
} // namespace
