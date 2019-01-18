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

#ifndef UNITTEST_MESHFIXTURE_HPP
#define UNITTEST_MESHFIXTURE_HPP

#include "mpi.h"
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace unit_test_util
{

class MeshFixtureNoTest
{
protected:
    MeshFixtureNoTest()
    : communicator(MPI_COMM_WORLD), metaData(nullptr), bulkData(nullptr)
    {
        allocate_meta();
    }

    MeshFixtureNoTest(unsigned spatial_dim)
    : communicator(MPI_COMM_WORLD), metaData(nullptr), bulkData(nullptr)
    {
        allocate_meta(spatial_dim);
    }

    MeshFixtureNoTest(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
    : communicator(MPI_COMM_WORLD), metaData(nullptr), bulkData(nullptr)
    {
        allocate_meta(spatial_dim, entityRankNames);
    }

    virtual ~MeshFixtureNoTest()
    {
        delete bulkData;
        delete metaData;
        bulkData = nullptr;
        metaData = nullptr;
    }

    void setup_empty_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        allocate_bulk(auraOption);
    }

    virtual void setup_mesh(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        allocate_bulk(auraOption);
        stk::io::fill_mesh(meshSpecification, *bulkData);
    }

    void setup_mesh_with_cyclic_decomp(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        allocate_bulk(auraOption);
        stk::unit_test_util::generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(meshSpecification,*bulkData,"cyclic");
    }

    MPI_Comm get_comm() const
    {
        return communicator;
    }

    void reset_mesh()
    {
        delete bulkData;
        delete metaData;

        bulkData = nullptr;
        metaData = nullptr;
    }

    int get_parallel_rank() const
    {
        return stk::parallel_machine_rank(get_comm());
    }

    int get_parallel_size() const
    {
        return stk::parallel_machine_size(get_comm());
    }

    virtual stk::mesh::MetaData& get_meta()
    {
        ThrowRequireMsg(metaData!=nullptr, "Unit test error. Trying to get meta data before it has been initialized.");
        return *metaData;
    }

    stk::mesh::BulkData& get_bulk()
    {
        ThrowRequireMsg(bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulkData;
    }

    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        if(nullptr == metaData)
            allocate_meta();

        bulkData = new stk::mesh::BulkData(get_meta(), communicator, auraOption);
    }

    virtual void allocate_meta(unsigned spatialDim = 3, const std::vector<std::string>& entityRankNames = {})
    {
        ThrowRequireMsg(metaData==nullptr, "Unit test error. Trying to reset non NULL meta data.");
        metaData = new stk::mesh::MetaData(spatialDim, entityRankNames);
    }

    void set_bulk(stk::mesh::BulkData *inBulkData)
    {
        ThrowRequireMsg(bulkData==nullptr, "Unit test error. Trying to reset non NULL bulk data.");
        bulkData = inBulkData;
    }

    void delete_meta()
    {
        ThrowRequireMsg(bulkData==nullptr, "Unit test error. Trying to delete meta with non NULL bulk data.");
        delete metaData;
        metaData = nullptr;
    }

protected:
    MPI_Comm communicator;
    stk::mesh::MetaData *metaData = nullptr;
    stk::mesh::BulkData *bulkData = nullptr;
};

class MeshFixture : public MeshFixtureNoTest, public ::testing::Test {
 protected:
  MeshFixture(){}
  MeshFixture(unsigned spatial_dim) : MeshFixtureNoTest(spatial_dim) {}
  MeshFixture(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
  : MeshFixtureNoTest(spatial_dim,entityRankNames){}

};

class MeshTestFixture : public MeshFixture
{
protected:
    void run_test_on_num_procs(int numProcs, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        if(stk::parallel_machine_size(get_comm()) == numProcs)
        {
            run_test(auraOption);
        }
    }

    void run_test_on_num_procs_or_less(int numProcs, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        if(stk::parallel_machine_size(get_comm()) <= numProcs)
        {
            run_test(auraOption);
        }
    }

    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption) = 0;
};

}}

#endif

