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

#ifndef stk_util_parallel_MPITagManager
#define stk_util_parallel_MPITagManager

#include <map>
#include <set>
#include <memory>
#include <iostream>
#include "stk_util/parallel/Parallel.hpp"   // for MPI


namespace stk {

class MPITag;

namespace Impl {
  class MPITagData;

  int delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);
}


// class for managing which MPI tags are in use on which communicator
// Note: all operations on this class are logically collective: they
//       must be called by all ranks on the given communicator in the
//       same order.
class MPITagManager
{
  public:
    MPITagManager();

    ~MPITagManager() {}

    MPITagManager(const MPITagManager&) = delete;

    MPITagManager& operator=(const MPITagManager&) = delete;

    using CommKey = int;

    // get an unused tag
    MPITag get_tag(MPI_Comm comm);

    // try to get tag given tag, or a nearby one if that specified one is unavailable
    MPITag get_tag(MPI_Comm comm, int tagHint);

    // sets the CommKey for the given comm if it hasn't been set already, returning the CommKey
    CommKey get_comm_key(MPI_Comm comm);

  private:

    struct CommRecord
    {
      CommRecord() = default;

      CommRecord(const CommRecord&) = delete;

      CommRecord& operator=(const CommRecord&) = delete;

      std::set<int> tags;
      std::map<int, std::weak_ptr<Impl::MPITagData>> tagPtrs;  // maps tag values to the MPITagData objects (pointers)
    };

    void free_tag(MPITag& tag);

    void free_tag(Impl::MPITagData& tag);

    const CommKey* generate_comm_key();

    void free_comm_key(CommKey* key);

    // checks that new_val is the same on all procs (in debug mode only)
    void debug_check(MPI_Comm comm, int newVal);

    std::map<CommKey, CommRecord> m_tags;
    std::set<CommKey> m_usedCommKeys;

    const int m_tagMin = 1024;
    int m_tagMax = 32767;
    int m_mpiAttrKey;

    friend Impl::MPITagData;
    friend MPITag;
    friend int Impl::delete_mpi_comm_key(MPI_Comm comm,int comm_keyval, void* attribute_val, void* extra_state);
};




namespace Impl {

class MPITagData
{
  public:

    MPITagData(MPITagManager* manager, MPITagManager::CommKey comm_key, int tag) :
      m_manager(manager),
      m_commKey(comm_key),
      m_tag(tag)
    {}

    ~MPITagData();

    int get_tag() { return m_tag;}

    MPITagManager::CommKey get_comm_key() { return m_commKey;}

    void set_free() { m_isFree = true;}

    bool is_free() { return m_isFree; }

  private:
    MPITagManager* m_manager;
    MPITagManager::CommKey m_commKey;
    int m_tag;
    bool m_isFree = false;
};
}

// an MPI tag that is currently in use.  Note that it is implicitly convertable
// to int, so it can be passed directly into MPI routines.
class MPITag
{
  public:
    explicit MPITag(std::shared_ptr<Impl::MPITagData> data) :
      m_data(data)
    {}

    operator int() const { return m_data->get_tag(); }

  private:
    void set_free()
    {
      m_data->set_free();
    }

    std::shared_ptr<Impl::MPITagData> m_data;

    friend MPITagManager;

    friend bool operator==(const MPITag& lhs, const MPITag& rhs);
    friend bool operator!=(const MPITag& lhs, const MPITag& rhs);
    friend std::ostream& operator<<(std::ostream& os, const MPITag& tag);
};

inline bool operator==(const MPITag& lhs, const MPITag& rhs)
{
  return static_cast<int>(lhs) == static_cast<int>(rhs) &&
         lhs.m_data->get_comm_key() == rhs.m_data->get_comm_key() &&
         !lhs.m_data->is_free() && !rhs.m_data->is_free();
}

inline bool operator!=(const MPITag& lhs, const MPITag& rhs)
{
  return !(lhs == rhs);
}

inline std::ostream&  operator<<(std::ostream& os, const MPITag& tag)
{
  os << "MPITag with value " << tag.m_data->get_tag() << " on Comm with key "
     << tag.m_data->get_comm_key() << ", is_free = " << tag.m_data->is_free();

  return os;
}




MPITagManager& get_mpi_tag_manager();

}

#endif