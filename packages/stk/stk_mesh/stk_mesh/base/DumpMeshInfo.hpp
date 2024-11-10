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
//
#ifndef DUMPMESHINFO_HPP
#define DUMPMESHINFO_HPP

#include <iostream>
#include <string>
#include "stk_mesh/base/Types.hpp"

namespace stk::mesh {
class BulkData;
class MetaData;
struct EntityKey;
class Bucket;
class Ghosting;
}

namespace stk::mesh::impl {

std::ostream & print_entity_key(std::ostream & os, const MetaData & meta_data, const EntityKey & key);
std::string print_entity_key(const MetaData & meta_data, const EntityKey & key);
bool print_comm_data_for_entity_in_ghosting(const BulkData & bulk, const Ghosting & ghosting,
                                            const EntityKey & entityKey, std::ostream & out);
void print_comm_data_for_entity(const BulkData & bulk, const EntityKey & entityKey, std::ostream & out);
void print_field_data_for_entity(const BulkData & bulk, const MeshIndex & meshIndex, std::ostream & out);
void print_entity_connectivity(const BulkData & bulk, const MeshIndex & meshIndex, std::ostream & out);
void print_bucket_parts(const BulkData & bulk, const Bucket * bucket, std::ostream & out);
void print_entity_offset_and_state(const BulkData & bulk, const MeshIndex & meshIndex, std::ostream & out);
void print_connectivity_of_rank(const BulkData & M, const Entity & entity, EntityRank connectedRank, std::ostream & out);

void dump_all_meta_info(const MetaData & meta, std::ostream & out);

void dump_mesh_bucket_info(const BulkData & bulk, std::ostream & out, Bucket * bucket);
void dump_all_mesh_info(const BulkData & bulk, std::ostream & out);
void dump_mesh_per_proc(const BulkData & bulk, const std::string & fileNamePrefix);

void dump_partition_summary(const BulkData & bulk, std::ostream & out);
void dump_partition_summary_per_proc(const BulkData & bulk, const std::string & fileNamePrefix);

std::vector<int> get_bucket_size_histogram(const BulkData & bulk, unsigned binWidth = 8);
void dump_bucket_size_histogram(const BulkData & bulk, std::ostream & out, unsigned binWidth = 8);
void dump_bucket_size_histogram_per_proc(const BulkData & bulk, const std::string& fileNamePrefix, unsigned binWidth = 8);
void dump_global_bucket_size_histogram(const BulkData & bulk, const std::string& fileName, unsigned binWidth = 8);

}

#endif
