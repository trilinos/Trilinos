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

#ifndef STK_MESH_FOREACHENTITY_HPP
#define STK_MESH_FOREACHENTITY_HPP

#include <stk_util/stk_config.h>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/baseImpl/ForEachEntityLoopAbstractions.hpp>

namespace stk {
namespace mesh {

class BulkData;

template <typename ALGORITHM_TO_RUN_PER_ENTITY>
void for_each_entity_run(const BulkData &mesh,
                         stk::topology::rank_t rank,
                         const Selector& selector,
                         const ALGORITHM_TO_RUN_PER_ENTITY &functor)
{
    impl::for_each_selected_entity_run(mesh, rank, selector, functor);
}

template <typename ALGORITHM_TO_RUN_PER_ENTITY>
void for_each_entity_run(const BulkData &mesh,
                         stk::topology::rank_t rank,
                         const ALGORITHM_TO_RUN_PER_ENTITY &functor)
{
    Selector selectAll = mesh.mesh_meta_data().universal_part();
    impl::for_each_selected_entity_run(mesh, rank, selectAll, functor);
}

template <typename ALGORITHM_TO_RUN_PER_ENTITY>
void for_each_entity_run_no_threads(const BulkData &mesh,
                         stk::topology::rank_t rank,
                         const Selector& selector,
                         const ALGORITHM_TO_RUN_PER_ENTITY &functor)
{
    impl::for_each_selected_entity_run_no_threads(mesh, rank, selector, functor);
}

template <typename ALGORITHM_TO_RUN_PER_ENTITY>
void for_each_entity_run_no_threads(const BulkData &mesh,
                         stk::topology::rank_t rank,
                         const ALGORITHM_TO_RUN_PER_ENTITY &functor)
{
    Selector selectAll = mesh.mesh_meta_data().universal_part();
    impl::for_each_selected_entity_run_no_threads(mesh, rank, selectAll, functor);
}

template <typename ALGORITHM_TO_RUN_PER_ENTITY>
void for_each_entity_run_with_nodes(const BulkData &mesh,
                                    stk::topology::rank_t rank,
                                    const Selector& selector,
                                    const ALGORITHM_TO_RUN_PER_ENTITY &functor)
{
    impl::for_each_selected_entity_run_with_nodes(mesh, rank, selector, functor);
}

} // namespace mesh
} // namespace stk

#endif

