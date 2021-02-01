// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioss_FaceGenerator_h
#define IOSS_Ioss_FaceGenerator_h

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef> // for size_t

#define FG_USE_ROBIN
#if defined FG_USE_STD
#include <unordered_set>
#elif defined FG_USE_HOPSCOTCH
#include <hopscotch_set.h>
#elif defined FG_USE_ROBIN
#include <robin_set.h>
#endif

#include <utility>

namespace Ioss {
  class Region;

  class Face
  {
  public:
    Face() = default;
    Face(size_t id, std::array<size_t, 4> conn) : hashId_(id), connectivity_(conn) {}
    Face(std::array<size_t, 4> conn);

    void add_element(size_t element_id) const
    {
      if (elementCount_ < 2) {
        element[elementCount_++] = element_id;
      }
      else {
        face_element_error(element_id);
      }
    }

    void face_element_error(size_t element_id) const;

    size_t hashId_{0};

    // NOTE: Not used at all by `Face` or `FaceGenerator` class, but are used by
    // skinner to give a consistent element id in cases where there
    // is a hash collision (face.id).

    // NOTE: For interior faces, this will not be the same value for each
    // face where the `hashId_` *will* be consistent for interior faces.
    // Should only use this as an id if `elementCount_` is 1.

    // NOTE: This could be used to do parallel or block boundary
    // collision since it is calculated as 10*element_id + local_face,
    // you could recover element_id and local_face and then set up
    // parallel communication maps.  May need to save the proc it is
    // shared with also (which is available in git history)
    mutable std::array<size_t, 2> element{};
    mutable int                   elementCount_{0}; // Should be max of 2 solid elements...
    std::array<size_t, 4>         connectivity_{};
  };

  struct FaceHash
  {
    size_t operator()(const Face &face) const { return face.hashId_; }
  };

  struct FaceEqual
  {
    bool operator()(const Face &left, const Face &right) const
    {
      if (left.hashId_ != right.hashId_) {
        return false;
      }
      // Hash (hashId_) is equal
      // Check whether same vertices (can be in different order)
      // Most (All?) of the time, there are no hashId_ collisions, so this test will not
      // find a difference and the function will return 'true'
      // However, for some reason, removing this check does not change the execution time
      // appreiciably...
      for (auto lvert : left.connectivity_) {
        if (std::find(right.connectivity_.cbegin(), right.connectivity_.cend(), lvert) ==
            right.connectivity_.cend()) {
          // Not found, therefore not the same.
          return false;
        }
      }
      return true;
    }
  };

#if defined FG_USE_STD
  using FaceUnorderedSet = std::unordered_set<Face, FaceHash, FaceEqual>;
#elif defined FG_USE_HOPSCOTCH
  using FaceUnorderedSet = tsl::hopscotch_set<Face, FaceHash, FaceEqual>;
  // using FaceUnorderedSet = tsl::hopscotch_pg_set<Face, FaceHash, FaceEqual>;
#elif defined FG_USE_ROBIN
  using FaceUnorderedSet = tsl::robin_set<Face, FaceHash, FaceEqual>;
  // using FaceUnorderedSet = tsl::robin_pg_set<Face, FaceHash, FaceEqual>;
#endif
  class FaceGenerator
  {
  public:
    explicit FaceGenerator(Ioss::Region &region);
    ~FaceGenerator() = default;

    static size_t id_hash(size_t global_id);

    template <typename INT> void generate_faces(INT /*dummy*/, bool block_by_block = false);
    FaceUnorderedSet &           faces(const std::string &name = "ALL") { return faces_[name]; }

    //! Given a local node id (0-based), return the hashed value.
    size_t node_id_hash(size_t local_node_id) const { return hashIds_[local_node_id]; }

  private:
    template <typename INT> void hash_node_ids(const std::vector<INT> &node_ids);
    template <typename INT> void generate_block_faces(INT /*dummy*/);
    template <typename INT> void generate_model_faces(INT /*dummy*/);

    Ioss::Region &                          region_;
    std::map<std::string, FaceUnorderedSet> faces_;
    std::vector<size_t>                     hashIds_;
  };

} // namespace Ioss

#endif
