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

#include <stk_mesh/base/Connectivity.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <algorithm>

namespace stk {
namespace mesh {

Connectivity::Connectivity()
 : m_connectivity(nullptr),
   m_offsetsAndOrdinals(nullptr)
{
  m_offsetsAndOrdinals = new ConnectivityOrdinal[m_capacityIndex + 1];
  for(unsigned i=0; i<m_capacityIndex; ++i) {
    m_offsetsAndOrdinals[i] = 0;
  }
  m_offsetsAndOrdinals[m_capacityIndex] = 0;
}

Connectivity::~Connectivity()
{
  delete [] m_connectivity;
  delete [] m_offsetsAndOrdinals;
}

int Connectivity::add_connectivity(EntityRank connectedRank,
                                   Entity connectedEntity,
                                   ConnectivityOrdinal ordinal,
                                   bool allowDuplicateOrdinals)
{
  ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
  int insertLocation = find_insert_location(connectedRank, connectedEntity,
                                            ordinal, allowDuplicateOrdinals);
  if (insertLocation >= 0) {
    insert(connectedRank, connectedEntity, ordinal, insertLocation);
  }
  return insertLocation;
}

int Connectivity::replace_or_add_connectivity(EntityRank connectedRank,
                                              Entity connectedEntity,
                                              ConnectivityOrdinal ordinal)
{
  ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
  int location = find_ordinal(connectedRank, ordinal);
  if (location >= 0) {
    bool different = replace(connectedRank, connectedEntity, ordinal, location);
    if (!different) {
      location = -1;
    }
  }
  else {
    location = find_insert_location(connectedRank, connectedEntity, ordinal, false);
    insert(connectedRank, connectedEntity, ordinal, location);
  }
  return location;
}

int Connectivity::remove_connectivity(EntityRank connectedRank,
                                      Entity connectedEntity,
                                      ConnectivityOrdinal ordinal)
{ 
  ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
  const int len = num_connectivity(connectedRank);
  const ConnectivityOrdinal* ordBegin = begin_ordinals(connectedRank);
  const ConnectivityOrdinal* ordEnd = end_ordinals(connectedRank);
  
  const ConnectivityOrdinal* lb = std::lower_bound(ordBegin, ordEnd, ordinal);
  int location = lb - ordBegin;
  const bool foundOrdinal = (lb != ordEnd) && (*lb == ordinal);
  if (foundOrdinal) {
    if (connectedEntity.local_offset() > 0) {
      const Entity* conn = begin_connectivity(connectedRank);
      while(location < len && ordinal == ordBegin[location]) {
        if (conn[location] == connectedEntity) {
          break;
        } 
        ++location;
      }
      if (location >= len || conn[location] != connectedEntity) {
        return -1;
      } 
    } 
  }
  else {
    return -1;
  }
      
  if (location >= 0) {
    remove(connectedRank, location);
  } 
      
  return location;
}

int Connectivity::find_insert_location(EntityRank connectedRank,
                                       Entity connectedEntity,
                                       ConnectivityOrdinal ordinal,
                                       bool allowDuplicateOrdinals)
{
  const int len = num_connectivity(connectedRank);
  const ConnectivityOrdinal* ordBegin = begin_ordinals(connectedRank);
  const ConnectivityOrdinal* ordEnd = end_ordinals(connectedRank);
  
  const ConnectivityOrdinal* lb = std::lower_bound(ordBegin, ordEnd, ordinal);
  int insertLocation = lb - ordBegin;
  const bool alreadyExists = (lb != ordEnd) && (*lb == ordinal);
  if (!allowDuplicateOrdinals) {
    if (alreadyExists) {
      return -1;
    }
  }
  else {
    if (alreadyExists) {
      const Entity* conn = begin_connectivity(connectedRank);
      while(insertLocation < len && ordinal == ordBegin[insertLocation] &&
            conn[insertLocation].local_offset() < connectedEntity.local_offset()) {
        ++insertLocation;
      }
      if (insertLocation < len && conn[insertLocation] == connectedEntity) {
        return -1;
      }
    }
  } 
  return insertLocation;
}

int Connectivity::find_ordinal(EntityRank connectedRank,
                               ConnectivityOrdinal ordinal)
{
  const ConnectivityOrdinal* ordBegin = begin_ordinals(connectedRank);
  const ConnectivityOrdinal* ordEnd = end_ordinals(connectedRank);
  const ConnectivityOrdinal* ordIter = std::find(ordBegin, ordEnd, ordinal);
  return ordIter==ordEnd ? -1 : ordIter - ordBegin;
}

void Connectivity::do_the_insert(unsigned insertOffset,
                                 Entity newEntity,
                                 ConnectivityOrdinal newOrdinal)
{
  bool needNewAlloc = m_offsetsAndOrdinals[MAX_NUM_RANKS] == get_capacity();
  if (needNewAlloc) {
    unsigned newCapacity = std::max(2u, static_cast<unsigned>(1.5*get_capacity()));
    Entity* newConnectivity = new Entity[newCapacity];
    ConnectivityOrdinal* newOrdinals = new ConnectivityOrdinal[m_firstOrdinalIndex + newCapacity];
    for(unsigned i=0; i<m_capacityIndex; ++i) {
      newOrdinals[i] = m_offsetsAndOrdinals[i];
    }
    newOrdinals[m_capacityIndex] = newCapacity;

    for(unsigned i=0; i<insertOffset; ++i) {
      newConnectivity[i] = m_connectivity[i];
      newOrdinals[m_firstOrdinalIndex + i] = m_offsetsAndOrdinals[m_firstOrdinalIndex + i];
    }
    newConnectivity[insertOffset] = newEntity;
    newOrdinals[m_firstOrdinalIndex + insertOffset] = newOrdinal;
    for(unsigned i=insertOffset; i<m_offsetsAndOrdinals[MAX_NUM_RANKS]; ++i) {
      newConnectivity[i+1] = m_connectivity[i];
      newOrdinals[m_firstOrdinalIndex + i+1] = m_offsetsAndOrdinals[m_firstOrdinalIndex + i];
    }
    delete [] m_connectivity;
    delete [] m_offsetsAndOrdinals;
    m_connectivity = newConnectivity;
    m_offsetsAndOrdinals = newOrdinals;
  }
  else {
    for(unsigned i=m_offsetsAndOrdinals[MAX_NUM_RANKS]; i>insertOffset; --i) {
      m_connectivity[i] = m_connectivity[i-1];
      m_offsetsAndOrdinals[m_firstOrdinalIndex + i] = m_offsetsAndOrdinals[m_firstOrdinalIndex + i-1];
    }
    m_connectivity[insertOffset] = newEntity;
    m_offsetsAndOrdinals[m_firstOrdinalIndex + insertOffset] = newOrdinal;
  }
}

void Connectivity::insert(EntityRank connectedRank,
                          Entity connectedEntity,
                          ConnectivityOrdinal ordinal,
                          unsigned insertLocation)
{     
  const unsigned insertOffset = m_offsetsAndOrdinals[connectedRank]+insertLocation;
  do_the_insert(insertOffset, connectedEntity, ordinal);
  for(unsigned idx = connectedRank+1; idx<MAX_NUM_RANKS+1; ++idx) {
    ++m_offsetsAndOrdinals[idx];
  }
}

bool Connectivity::replace(EntityRank connectedRank,
                           Entity connectedEntity,
                           ConnectivityOrdinal ordinal,
                           unsigned location)
{     
  const unsigned offset = m_offsetsAndOrdinals[connectedRank]+location;
  bool different = m_connectivity[offset] != connectedEntity;
  if (different) {
    m_connectivity[offset] = connectedEntity;
  }
  ThrowAssertMsg(m_offsetsAndOrdinals[m_firstOrdinalIndex + offset] == ordinal, "Connectivity::replace ordinal doesn't match.");
  return different;
}

void Connectivity::remove(EntityRank connectedRank, unsigned location)
{
  const int overallEnd = m_offsetsAndOrdinals[MAX_NUM_RANKS];
  int overallLocation = m_offsetsAndOrdinals[connectedRank] + location;
  for(int i=overallLocation; i<overallEnd-1; ++i) {
    m_connectivity[i] = m_connectivity[i+1];
    m_offsetsAndOrdinals[m_firstOrdinalIndex + i] = m_offsetsAndOrdinals[m_firstOrdinalIndex + i+1];
  }
  for(unsigned idx = connectedRank+1; idx<=MAX_NUM_RANKS; ++idx) {
    --m_offsetsAndOrdinals[idx];
  }
}

} // namespace mesh
} // namespace stk

