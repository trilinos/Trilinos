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

#ifndef READWRITESIDESETTESTER_HPP_
#define READWRITESIDESETTESTER_HPP_

#include <vector>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/InputFile.hpp>   // for InputFile for m_input_files
#include <stk_util/parallel/Parallel.hpp> // for stk::parallel_machine_size
#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc
#include "Ioss_SideBlock.h"             // for SideBlock
#include "Ioss_SideSet.h"               // for SideSet, SideBlockContainer
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"

namespace stk{ namespace unit_test_util{ namespace sideset{

struct IdAndSideSet {
    int id;
    stk::mesh::SideSet sideSet;
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

    virtual void populate_mesh(bool delay_field_data_allocation = true)
    {
        stk::io::StkMeshIoBroker::populate_mesh(delay_field_data_allocation);
        extract_sideset_data_from_io();
    }

    void write_output_mesh(size_t output_file_index)
    {
      m_output_files[output_file_index]->write_output_mesh(bulk_data(), attributeFieldOrderingByPartOrdinal);
    }

private:
    int convert_to_zero_based_ordinal(int one_based_ordinal) const
    {
        return one_based_ordinal-1;
    }

    stk::mesh::SideSetEntry add_elem_side_pair(stk::mesh::Entity elem, int one_based_ordinal)
    {
        int zero_based_side_ordinal = convert_to_zero_based_ordinal(one_based_ordinal);
        return stk::mesh::SideSetEntry{elem, zero_based_side_ordinal};
    }

    void convert_elem_sides_pairs_into_sideset(const stk::mesh::BulkData& bulk, const std::vector<int>& elem_side, stk::mesh::SideSet& sideset)
    {
        for(size_t is=0; is<elem_side.size() / 2; ++is)
        {
            stk::mesh::Entity const elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elem_side[is*2]);
            if (bulk.is_valid(elem))
                sideset.push_back(add_elem_side_pair(elem, elem_side[is*2+1]));
        }
    }

    void convert_block_to_sideset(const stk::mesh::BulkData& bulk, const Ioss::SideBlock* block, stk::mesh::SideSet& sideset)
    {
        if (stk::io::include_entity(block))
        {
            std::vector<int> elem_side ;
            block->get_field_data("element_side", elem_side);
            convert_elem_sides_pairs_into_sideset(bulk, elem_side, sideset);
        }
    }

    void convert_ioss_sideset_to_stk_sideset(const stk::mesh::BulkData& bulk, const Ioss::SideSet* sset, stk::mesh::SideSet& sideset)
    {
        if(stk::io::include_entity(sset))
        {
            for (size_t i=0; i < sset->block_count(); i++)
            {
                Ioss::SideBlock *block = sset->get_block(i);
                convert_block_to_sideset(bulk, block, sideset);
            }
        }
    }



private:
    void extract_sideset_data_from_io()
    {
        stk::mesh::BulkData &bulk = this->bulk_data();
        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        for ( const Ioss::SideSet * sset : region->get_sidesets() )
        {
            stk::mesh::SideSet &sideSet = bulk.create_sideset(sset->get_property("id").get_int());
            convert_ioss_sideset_to_stk_sideset(bulk, sset, sideSet);
        }
    }
};

class BulkDataTester : public stk::mesh::BulkData
{
public:
    BulkDataTester(stk::mesh::MetaData & mesh_meta_data
                   , stk::ParallelMachine parallel)
    :   stk::mesh::BulkData(mesh_meta_data, parallel)
    {
    }
};

stk::mesh::SideSet get_stk_side_set(stk::mesh::BulkData &bulk, const ElemIdSideVector &ss);
SideSetData get_stk_side_set_data(stk::mesh::BulkData &bulk, const SideSetIdAndElemIdSidesVector &ssData);

void write_exo_file(BulkDataTester &bulkData, const std::string &filename);
void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode);

void load_mesh_and_fill_sideset_data(StkMeshIoBrokerTester &stkIo);
void setup_io_broker_for_read(stk::io::StkMeshIoBroker &stkIo, stk::mesh::BulkData &bulkData, std::string filename, ReadMode readMode);

void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& inputFileName, const std::string& outputFileName);

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


}
}
}
#endif
