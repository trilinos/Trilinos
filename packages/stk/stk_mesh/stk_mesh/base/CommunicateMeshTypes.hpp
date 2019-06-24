#ifndef STK_COMMUNICATE_MESH_TYPES_H
#define STK_COMMUNICATE_MESH_TYPES_H

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

#include <cassert>
#include <vector>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {

template<>
inline
CommBuffer & CommBuffer::pack<stk::mesh::PartVector>( const stk::mesh::PartVector & parts )
{
  unsigned num_parts = parts.size();
  std::vector<unsigned> part_ordinals;
  part_ordinals.reserve(num_parts);
  for(stk::mesh::PartVector::const_iterator it = parts.begin(); it != parts.end(); ++it)
  {
    const stk::mesh::Part * part = *it;
    part_ordinals.push_back(part->mesh_meta_data_ordinal());
  }

  pack(num_parts);
  pack(part_ordinals.data(), num_parts);

  return *this;
}

template<>
inline
CommBuffer & CommBuffer::pack<stk::mesh::Selector>( const stk::mesh::Selector & selector )
{
  assert(selector.is_all_unions());

  stk::mesh::PartVector parts;
  selector.get_parts(parts);
  pack(parts);

  return *this;
}

template<>
inline
CommBuffer & CommBuffer::unpack< std::pair<const stk::mesh::MetaData *, stk::mesh::PartVector *> >
( std::pair<const stk::mesh::MetaData *, stk::mesh::PartVector *> & pair )
{
  const stk::mesh::MetaData & meta = *pair.first;
  stk::mesh::PartVector & part_vector = *pair.second;

  unsigned num_parts;
  unpack(num_parts);
  std::vector<unsigned> part_ordinals(num_parts);
  unpack(part_ordinals.data(), num_parts);

  for(std::vector<unsigned>::const_iterator it = part_ordinals.begin(); it != part_ordinals.end(); ++it)
  {
    part_vector.push_back( &meta.get_part(*it) );
  }

  return *this;
}

template<>
inline
CommBuffer & CommBuffer::unpack< std::pair<const stk::mesh::MetaData *, stk::mesh::Selector *> >
( std::pair<const stk::mesh::MetaData *, stk::mesh::Selector *> & pair )
{
  const stk::mesh::MetaData * meta = pair.first;
  stk::mesh::Selector & selector = *pair.second;

  stk::mesh::PartVector parts;
  std::pair<const stk::mesh::MetaData *, stk::mesh::PartVector *> meta_parts_pair(meta, &parts);
  unpack(meta_parts_pair);

  selector = selectUnion(parts);

  return *this;
}

} // namespace stk

#endif
