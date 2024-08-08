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

#ifndef READWRITESIDESETTESTER_HPP_
#define READWRITESIDESETTESTER_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stddef.h>                                  // for size_t
#include <stk_io/InputFile.hpp>                      // for InputFile
#include <stk_io/StkMeshIoBroker.hpp>                // for StkMeshIoBroker, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/MetaData.hpp>                // for BulkData
#include <stk_mesh/base/Part.hpp>                // for BulkData
#include <stk_mesh/base/SideSetEntry.hpp>            // for SideSet, etc
#include <stk_util/parallel/Parallel.hpp>            // for ParallelMachine
#include <string>                                    // for string
#include <vector>                                    // for vector
#include "Ioss_Property.h"                           // for Property
#include "Ioss_Region.h"                             // for Region
#include "Ioss_SideBlock.h"                          // for SideBlock
#include "Ioss_SideSet.h"                            // for SideSet
#include "mpi.h"
#include "stk_io/IossBridge.hpp"                     // for include_entity
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_unit_test_utils/FaceTestingUtils.hpp"

namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk{ namespace unit_test_util{ namespace sideset{

struct IdAndSideSet {
    int id;
    stk::mesh::SideSet* sideSet = nullptr;

    ~IdAndSideSet() {
        delete sideSet;
    }
};
typedef std::vector<IdAndSideSet> SideSetData;

struct ElemIdSide {
    int elem_id;
    int side_ordinal;
};
typedef std::vector<ElemIdSide> ElemIdSideVector;

struct SideSetIdAndElemIdSides {
    int id;
    ElemIdSideVector sideSet;
};
typedef std::vector<SideSetIdAndElemIdSides> SideSetIdAndElemIdSidesVector;

struct ElemIdSideLess {
  inline bool operator()(const ElemIdSide &lhs, const ElemIdSide &rhs) const
  {
      if(lhs.elem_id < rhs.elem_id)
          return true;
      else if(lhs.elem_id == rhs.elem_id)
          return lhs.side_ordinal < rhs.side_ordinal;
      else
          return false;
  }
  inline ElemIdSideLess& operator=(const ElemIdSideLess& rhs);
};

enum ReadMode { READ_SERIAL_AND_DECOMPOSE, READ_ALREADY_DECOMPOSED };

class StkMeshIoBrokerTester : public stk::io::StkMeshIoBroker
{
public:
  StkMeshIoBrokerTester() {}

  virtual void write_output_mesh(size_t output_file_index)
  {
    m_outputFiles[output_file_index]->write_output_mesh(bulk_data(), attributeFieldOrderingByPartOrdinal);
  }
};

class BulkDataTester : public stk::mesh::BulkData
{
public:
  BulkDataTester(stk::mesh::MetaData & mesh_meta_data, stk::ParallelMachine parallel)
    : stk::mesh::BulkData(std::shared_ptr<stk::mesh::MetaData>(&mesh_meta_data, [](auto pointerWeWontDelete){}), parallel, stk::mesh::BulkData::AUTO_AURA)
  {
  }
};

stk::mesh::SideSet* get_stk_side_set(stk::mesh::BulkData &bulk, const ElemIdSideVector &ss);
SideSetData get_stk_side_set_data(stk::mesh::BulkData &bulk, const SideSetIdAndElemIdSidesVector &ssData);

void write_exo_file(BulkDataTester &bulkData, const std::string &filename);
void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode);

void load_mesh_and_fill_sideset_data(StkMeshIoBrokerTester &stkIo);
void setup_io_broker_for_read(stk::io::StkMeshIoBroker &stkIo, stk::mesh::BulkData &bulkData, std::string filename, ReadMode readMode);

void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& inputFileName, const std::string& outputFileName);

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& inputFileName, const std::string& outputFileName);

} // namespace simple_fields

void compare_sidesets(const std::string& inputFileName,
                      BulkDataTester &bulk1,
                      BulkDataTester &bulk2);

void compare_sidesets(const std::string& inputFileName,
                      stk::mesh::BulkData &bulk,
                      const SideSetIdAndElemIdSidesVector &expected);

void compare_sidesets(const std::string& input_file_name,
                      stk::mesh::BulkData &bulk,
                      const SideSetIdAndElemIdSidesVector &sideset,
                      const SideSetIdAndElemIdSidesVector &expected);


namespace simple_fields {

class StkMeshIoBrokerTester : public stk::io::StkMeshIoBroker
{
public:
  StkMeshIoBrokerTester() {}

  virtual void write_output_mesh(size_t output_file_index)
  {
    m_outputFiles[output_file_index]->write_output_mesh(bulk_data(), attributeFieldOrderingByPartOrdinal);
  }
};

class BulkDataTester : public stk::mesh::BulkData
{
public:
  BulkDataTester(stk::mesh::MetaData & mesh_meta_data, stk::ParallelMachine parallel)
    : stk::mesh::BulkData(std::shared_ptr<stk::mesh::MetaData>(&mesh_meta_data, [](auto pointerWeWontDelete){}), parallel, stk::mesh::BulkData::AUTO_AURA)
  {
  }
};

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::SideSet* get_stk_side_set(stk::mesh::BulkData &bulk, const ElemIdSideVector &ss);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
SideSetData get_stk_side_set_data(stk::mesh::BulkData &bulk, const SideSetIdAndElemIdSidesVector &ssData);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void write_exo_file(BulkDataTester &bulkData, const std::string &filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void load_mesh_and_fill_sideset_data(StkMeshIoBrokerTester &stkIo);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void setup_io_broker_for_read(stk::io::StkMeshIoBroker &stkIo, stk::mesh::BulkData &bulkData, std::string filename, ReadMode readMode);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& inputFileName, const std::string& outputFileName);

} // namespace simple_fields

}
}
}
#endif
