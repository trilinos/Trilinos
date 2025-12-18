// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef DEVICEFIELDDATAMANAGERBASE_HPP
#define DEVICEFIELDDATAMANAGERBASE_HPP

#include "stk_mesh/base/Types.hpp"
#include <cstddef>
#include <any>

namespace stk::mesh {

class FieldDataBase;

class DeviceFieldDataManagerBase
{
public:

  struct BucketShift {
    BucketShift(int _oldIndex, int _newIndex)
      : oldIndex(_oldIndex),
        newIndex(_newIndex)
    {}
    int oldIndex;
    int newIndex;
  };

  struct SizeAndCapacity
  {
    size_t size = 0;
    size_t capacity = 0;
  };

  DeviceFieldDataManagerBase() = default;
  virtual ~DeviceFieldDataManagerBase() = default;

  virtual bool update_all_bucket_allocations() = 0;
  virtual void update_host_bucket_pointers(Ordinal fieldOrdinal) = 0;
  virtual void reorder_and_resize_buckets(EntityRank rank, const FieldVector& fields,
                                          const std::vector<SizeAndCapacity>& newBucketSizes,
                                          const std::vector<BucketShift>& bucketShifts) = 0;
  virtual void swap_field_data(Ordinal fieldOrdinal1, Ordinal fieldOrdinal2) = 0;
  virtual void clear_bucket_is_modified(Ordinal fieldOrdinal) = 0;
  virtual std::any get_device_bucket_is_modified(Ordinal fieldOrdinal, int& fieldRankedOrdinal) = 0;
  virtual size_t get_num_bytes_allocated_on_field(const FieldBase& field) const = 0;
  virtual bool has_unified_device_storage(Ordinal fieldOrdinal) const = 0;

  virtual void set_device_field_meta_data(FieldDataBase& fieldDataBase) = 0;

  virtual void add_new_bucket(EntityRank rank,
                              unsigned bucketSize,
                              unsigned bucketCapacity,
                              const PartVector& parts,
                              bool deviceMeshMod = false) = 0;

  virtual size_t synchronized_count() const = 0;

};

}

#endif // DEVICEFIELDDATAMANAGERBASE_HPP
