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

#ifndef IOSS_Ioss_FaceGenerator_h
#define IOSS_Ioss_FaceGenerator_h

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef> // for size_t
#include <unordered_set>
#include <utility>

namespace Ioss {
  class Region;

  class Face
  {
  public:
    Face() : id_(0), elementCount_(0), sharedWithProc_(-1) {}
    Face(size_t id, std::array<size_t, 4> conn)
        : id_(id), elementCount_(0), sharedWithProc_(-1), connectivity_(std::move(conn))
    {
    }

    void add_element(size_t element_id) const
    {
      assert(elementCount_ < 2);
      element[elementCount_++] = element_id;
    }

    size_t         id_;
    mutable size_t element[2]{};
    mutable int    elementCount_; // Should be max of 2 solid elements...
    mutable int    sharedWithProc_;
    std::array<size_t, 4> connectivity_{};
  };

  struct FaceHash
  {
    size_t operator()(const Face &face) const { return face.id_; }
  };

  struct FaceEqual
  {
    bool operator()(const Face &left, const Face &right) const
    {
      if (left.id_ != right.id_) {
        return false;
      }
      // Hash (id_) is equal
      // Check whether same vertices (can be in different order)
      for (auto lvert : left.connectivity_) {
        if (std::find(right.connectivity_.begin(), right.connectivity_.end(), lvert) ==
            right.connectivity_.end()) {
          // Not found, therefore not the same.
          return false;
        }
      }
      return true;
    }
  };

  using FaceUnorderedSet = std::unordered_set<Face, FaceHash, FaceEqual>;
  class FaceGenerator
  {
  public:
    explicit FaceGenerator(Ioss::Region &region);

    template <typename INT> void generate_faces(INT /*dummy*/);
    FaceUnorderedSet &           faces() { return faces_; }
  private:
    Ioss::Region &   region_;
    FaceUnorderedSet faces_;
  };
} // namespace Ioss

#endif
