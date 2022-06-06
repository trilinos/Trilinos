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

#include "stk_util/stk_config.h"            // for STK_HAS_MPI
#ifdef STK_HAS_MPI

#include <map>
#include <set>
#include <forward_list>
#include <deque>
#include <memory>
#include <iostream>
#include <iterator>
#include "MPICommKey.hpp"
#include "stk_util/parallel/Parallel.hpp"   // for MPI
#include "stk_util/parallel/MPICommKey.hpp"  // for MPIKeyManager
#include "stk_util/parallel/MPITag.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, ThrowRequire
#include "stk_util/parallel/CommReplacer.hpp"
#include "stk_util/parallel/CommTagInUseList.hpp"


namespace stk {

// class for managing which MPI tags are in use on which communicator
// Note: all operations on this class are logically collective: they
//       must be called by all ranks on the given communicator in the
//       same order.
class MPITagManager
{
  protected:
    //TODO: writeup delayCount theory into manual
    MPITagManager(int deletionGroupSize, int delayCount);

  public:

    ~MPITagManager();

    MPITagManager(const MPITagManager&) = delete;

    MPITagManager& operator=(const MPITagManager&) = delete;

    MPITagManager(MPITagManager&&) = delete;

    MPITagManager& operator=(MPITagManager&&) = delete;

    // try to get tag given tag, or a nearby one if that specified one is unavailable
    // Note: this is a collective operation
    MPITag get_tag(MPI_Comm comm, int tagHint=MPI_ANY_TAG);

  private:
    //-------------------------------------------------------------------------
    // Tag generation

    int get_new_tag(impl::CommTagInUseList& commData, int tagHint);

    int get_any_tag(impl::CommTagInUseList& commData);

    int new_tag_search(impl::CommTagInUseList& commData, int tagHint);

    //-------------------------------------------------------------------------
    // Tag and Comm freeing

    void free_tag_local(impl::MPITagData& tag);

    void erase_comm(MPI_Comm comm);

    void check_same_value_on_all_procs_debug_only(MPI_Comm comm, int newVal);

    std::shared_ptr<impl::MPIKeyManager> m_keyManager;
    std::map<MPI_Comm, impl::CommTagInUseList, impl::CommCompare> m_commData;
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    impl::CommReplacer m_commReplacer;
#endif

    const unsigned int m_deletionGroupSize;  // number of tags to free as a group
    const unsigned int m_delayCount;         // wait on MPI_Ibarriers after this many
                                             // entries to get_tag() with the same MPI_Comm
    const int m_tagMin = 1024;
    int m_tagMax = 32767;
    impl::MPIKeyManager::CallerUID m_callbackUID = 0;

    friend impl::MPITagData;
    friend MPITagManager& get_mpi_tag_manager();
};


MPITagManager& get_mpi_tag_manager();
}


#endif
#endif
