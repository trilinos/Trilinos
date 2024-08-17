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

#ifndef UNITTEST_MESHFIXTURE_HPP
#define UNITTEST_MESHFIXTURE_HPP

#include "mpi.h"
#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_util/util/ReportHandler.hpp"          // for set_report_handler

namespace stk
{
namespace unit_test_util
{

class MeshFixtureNoTest
{
protected:
    MeshFixtureNoTest()
    : communicator(MPI_COMM_WORLD),
      m_spatialDim(3),
      m_entityRankNames()
    {
    }

    MeshFixtureNoTest(unsigned spatial_dim)
    : communicator(MPI_COMM_WORLD),
      m_spatialDim(spatial_dim),
      m_entityRankNames()
    {
    }

    MeshFixtureNoTest(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
    : communicator(MPI_COMM_WORLD),
      m_spatialDim(spatial_dim),
      m_entityRankNames(entityRankNames)
    {
    }

    MeshFixtureNoTest(unsigned spatial_dim, stk::mesh::BulkData::AutomaticAuraOption auraOption,
                      MPI_Comm comm = MPI_COMM_WORLD)
      : communicator(comm),
        m_spatialDim(spatial_dim),
        m_entityRankNames()
    {
      setup_empty_mesh(auraOption);
    }

    virtual ~MeshFixtureNoTest()
    {

    }

    void set_spatial_dimension(unsigned spatialDim)
    {
      m_spatialDim = spatialDim;
    }

    void setup_empty_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                          unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                          unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        allocate_bulk(auraOption, initialBucketCapacity, maximumBucketCapacity);
    }

    virtual void setup_mesh(const std::string &meshSpecification,
                            stk::mesh::BulkData::AutomaticAuraOption auraOption,
                            unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                            unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        allocate_bulk(auraOption, initialBucketCapacity, maximumBucketCapacity);
        stk::io::fill_mesh(meshSpecification, *bulkData);
    }

    void setup_mesh_with_cyclic_decomp(const std::string &meshSpecification,
                                       stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                       unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                                       unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        allocate_bulk(auraOption, initialBucketCapacity, maximumBucketCapacity);
        stk::unit_test_util::generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(meshSpecification,*bulkData,"cyclic");
    }

    MPI_Comm get_comm() const
    {
        return communicator;
    }

    void reset_mesh()
    {
        bulkData.reset();
        metaData.reset();
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
        STK_ThrowRequireMsg(metaData!=nullptr, "Unit test error. Trying to get meta data before it has been initialized.");
        return *metaData;
    }

    virtual stk::mesh::BulkData& get_bulk()
    {
        STK_ThrowRequireMsg(bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulkData;
    }

    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                               unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                               unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        stk::mesh::MeshBuilder builder(communicator);
        builder.set_spatial_dimension(m_spatialDim);
        builder.set_entity_rank_names(m_entityRankNames);
        builder.set_aura_option(auraOption);
        builder.set_initial_bucket_capacity(initialBucketCapacity);
        builder.set_maximum_bucket_capacity(maximumBucketCapacity);

        if(nullptr == metaData) {
          metaData = builder.create_meta_data();
        }

        if(nullptr == bulkData) {
          bulkData = builder.create(metaData);
          m_auraOption = auraOption;
          m_initialBucketCapacity = initialBucketCapacity;
          m_maximumBucketCapacity = maximumBucketCapacity;
        }

        STK_ThrowRequireMsg((auraOption == m_auraOption) &&
                            (initialBucketCapacity == m_initialBucketCapacity) &&
                            (maximumBucketCapacity == m_maximumBucketCapacity),
           "allocate_bulk() being called with different arguments from previous call:\n"
           "    auraOption = " << auraOption << " (previously: " << m_auraOption << ")\n"
           "    initialBucketCapacity = " << initialBucketCapacity << " (previously: " << m_initialBucketCapacity << ")\n"
           "    maximumBucketCapacity = " << maximumBucketCapacity << " (previously: " << m_maximumBucketCapacity << ")");
    }

    void set_meta(std::shared_ptr<stk::mesh::MetaData> inMetaData)
    {
        STK_ThrowRequireMsg(metaData==nullptr, "Unit test error. Trying to reset non NULL meta data.");
        metaData = inMetaData;
    }

    void set_bulk(std::shared_ptr<stk::mesh::BulkData> inBulkData)
    {
        STK_ThrowRequireMsg(bulkData==nullptr, "Unit test error. Trying to reset non NULL bulk data.");
        bulkData = inBulkData;

        STK_ThrowRequireMsg(metaData==nullptr || metaData==bulkData->mesh_meta_data_ptr(),
                        "Unit test error. Trying to reset non NULL meta data.");
    }

protected:
    MPI_Comm communicator;
    unsigned m_spatialDim;
    std::vector<std::string> m_entityRankNames;
    std::shared_ptr<stk::mesh::MetaData> metaData;
    std::shared_ptr<stk::mesh::BulkData> bulkData;

    stk::mesh::BulkData::AutomaticAuraOption m_auraOption{stk::mesh::BulkData::AUTO_AURA};
    unsigned m_initialBucketCapacity = 0;
    unsigned m_maximumBucketCapacity = 0;
};

class MeshFixture : public MeshFixtureNoTest, public ::ngp_testing::Test {
 protected:
  MeshFixture(){}
  MeshFixture(unsigned spatial_dim) : MeshFixtureNoTest(spatial_dim) {}
  MeshFixture(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
  : MeshFixtureNoTest(spatial_dim,entityRankNames){}

};

class MeshFixture2D : public MeshFixtureNoTest, public ::ngp_testing::Test {
 protected:
  MeshFixture2D() : MeshFixtureNoTest(2) {}
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

inline void delete_mesh(const std::string & baseFileName)
{
  const int pSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (pSize == 1) {
    unlink(baseFileName.c_str());
  }
  else {
    for (int proc = 0; proc < pSize; ++proc) {
      const std::string fileName = baseFileName + "." + std::to_string(pSize) + "." + std::to_string(proc);
      unlink(fileName.c_str());
    }
  }
}

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
MeshFixtureNoTest
{
protected:
    MeshFixtureNoTest()
    : communicator(MPI_COMM_WORLD),
      m_spatialDim(3),
      m_entityRankNames()
    {
    }

    MeshFixtureNoTest(unsigned spatial_dim)
    : communicator(MPI_COMM_WORLD),
      m_spatialDim(spatial_dim),
      m_entityRankNames()
    {
    }

    MeshFixtureNoTest(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
    : communicator(MPI_COMM_WORLD),
      m_spatialDim(spatial_dim),
      m_entityRankNames(entityRankNames)
    {
    }

    MeshFixtureNoTest(unsigned spatial_dim, stk::mesh::BulkData::AutomaticAuraOption auraOption,
                      MPI_Comm comm = MPI_COMM_WORLD)
      : communicator(comm),
        m_spatialDim(spatial_dim),
        m_entityRankNames()
    {
      setup_empty_mesh(auraOption);
    }

    virtual ~MeshFixtureNoTest()
    {

    }

    void set_spatial_dimension(unsigned spatialDim)
    {
      m_spatialDim = spatialDim;
    }

    void setup_empty_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                          unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                          unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        allocate_bulk(auraOption, initialBucketCapacity, maximumBucketCapacity);
    }

    virtual void setup_mesh(const std::string &meshSpecification,
                            stk::mesh::BulkData::AutomaticAuraOption auraOption,
                            unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                            unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        allocate_bulk(auraOption, initialBucketCapacity, maximumBucketCapacity);
        stk::io::fill_mesh(meshSpecification, *bulkData);
    }

    void setup_mesh_with_cyclic_decomp(const std::string &meshSpecification,
                                       stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                       unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                                       unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        allocate_bulk(auraOption, initialBucketCapacity, maximumBucketCapacity);
        stk::unit_test_util::generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(meshSpecification,*bulkData,"cyclic");
    }

    MPI_Comm get_comm() const
    {
        return communicator;
    }

    void reset_mesh()
    {
        bulkData.reset();
        metaData.reset();
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
        STK_ThrowRequireMsg(metaData!=nullptr, "Unit test error. Trying to get meta data before it has been initialized.");
        return *metaData;
    }

    virtual stk::mesh::BulkData& get_bulk()
    {
        STK_ThrowRequireMsg(bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulkData;
    }

    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                               unsigned initialBucketCapacity = mesh::get_default_initial_bucket_capacity(),
                               unsigned maximumBucketCapacity = mesh::get_default_maximum_bucket_capacity())
    {
        stk::mesh::MeshBuilder builder(communicator);
        builder.set_spatial_dimension(m_spatialDim);
        builder.set_entity_rank_names(m_entityRankNames);
        builder.set_aura_option(auraOption);
        builder.set_initial_bucket_capacity(initialBucketCapacity);
        builder.set_maximum_bucket_capacity(maximumBucketCapacity);

        if(nullptr == metaData) {
          metaData = builder.create_meta_data();
        }

        if(nullptr == bulkData) {
          bulkData = builder.create(metaData);
          m_auraOption = auraOption;
          m_initialBucketCapacity = initialBucketCapacity;
          m_maximumBucketCapacity = maximumBucketCapacity;
        }

        STK_ThrowRequireMsg((auraOption == m_auraOption) &&
                            (initialBucketCapacity == m_initialBucketCapacity) &&
                            (maximumBucketCapacity == m_maximumBucketCapacity),
           "allocate_bulk() being called with different arguments from previous call:\n"
           "    auraOption = " << auraOption << " (previously: " << m_auraOption << ")\n"
           "    initialBucketCapacity = " << initialBucketCapacity << " (previously: " << m_initialBucketCapacity << ")\n"
           "    maximumBucketCapacity = " << maximumBucketCapacity << " (previously: " << m_maximumBucketCapacity << ")");
    }

    void set_meta(std::shared_ptr<stk::mesh::MetaData> inMetaData)
    {
        STK_ThrowRequireMsg(metaData==nullptr, "Unit test error. Trying to reset non NULL meta data.");
        metaData = inMetaData;
    }

    void set_bulk(std::shared_ptr<stk::mesh::BulkData> inBulkData)
    {
        STK_ThrowRequireMsg(bulkData==nullptr, "Unit test error. Trying to reset non NULL bulk data.");
        bulkData = inBulkData;

        STK_ThrowRequireMsg(metaData==nullptr || metaData==bulkData->mesh_meta_data_ptr(),
                        "Unit test error. Trying to reset non NULL meta data.");
    }

protected:
    MPI_Comm communicator;
    unsigned m_spatialDim;
    std::vector<std::string> m_entityRankNames;
    std::shared_ptr<stk::mesh::MetaData> metaData;
    std::shared_ptr<stk::mesh::BulkData> bulkData;

    stk::mesh::BulkData::AutomaticAuraOption m_auraOption{stk::mesh::BulkData::AUTO_AURA};
    unsigned m_initialBucketCapacity = 0;
    unsigned m_maximumBucketCapacity = 0;
};

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
MeshFixture : public MeshFixtureNoTest, public ::ngp_testing::Test {
 protected:
  MeshFixture(){}
  MeshFixture(unsigned spatial_dim) : MeshFixtureNoTest(spatial_dim) {}
  MeshFixture(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
  : MeshFixtureNoTest(spatial_dim,entityRankNames){}

};

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
MeshFixture2D : public MeshFixtureNoTest, public ::ngp_testing::Test {
 protected:
  MeshFixture2D() : MeshFixtureNoTest(2) {}
};

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
MeshTestFixture : public MeshFixture
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

inline void delete_mesh(const std::string & baseFileName)
{
  const int pSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (pSize == 1) {
    unlink(baseFileName.c_str());
  }
  else {
    for (int proc = 0; proc < pSize; ++proc) {
      const std::string fileName = baseFileName + "." + std::to_string(pSize) + "." + std::to_string(proc);
      unlink(fileName.c_str());
    }
  }
}

}  // namespace simple_fields

}}

#endif

