// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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
//
#ifndef _UnitTestSkinMeshUseCaseUtils_hpp_
#define _UnitTestSkinMeshUseCaseUtils_hpp_

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ElemGraphTestUtils.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_unit_tests/stk_mesh/SetupKeyholeMesh.hpp>

#include <stk_unit_tests/stk_mesh_fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_unit_tests/stk_mesh_fixtures/heterogeneous_mesh.hpp>
#include <stk_unit_tests/stk_mesh_fixtures/degenerate_mesh.hpp>


namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;


inline stk::mesh::PartVector get_skin_parts(stk::mesh::MetaData &meta)
{
    stk::mesh::PartVector skin_parts;
    skin_parts.push_back(meta.get_part("skin"));
    skin_parts.push_back(meta.get_part("active"));
    return skin_parts;
}

inline void make_2_hex_mesh_with_element1_inactive(stk::mesh::BulkData& bulkData)
{
    stk::mesh::MetaData &meta = bulkData.mesh_meta_data();
    meta.declare_part_with_topology("skin", stk::topology::QUAD_4);
    stk::mesh::Part& active = meta.declare_part("active");
    stk::io::fill_mesh("generated:1x1x2", bulkData);

    stk::mesh::EntityVector entitiesToMakeActive;
    std::vector<stk::mesh::PartVector> add_parts;
    std::vector<stk::mesh::PartVector> rm_parts;

    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Entity element1 = bulkData.get_entity(rank, 1);

    if(bulkData.is_valid(element1) && bulkData.bucket(element1).owned())
    {
        entitiesToMakeActive.push_back(element1);
        add_parts.push_back(stk::mesh::PartVector(1, &active));
        rm_parts.push_back(stk::mesh::PartVector());
    }
    bulkData.batch_change_entity_parts(entitiesToMakeActive, add_parts, rm_parts);
}

} // namespace

#endif


