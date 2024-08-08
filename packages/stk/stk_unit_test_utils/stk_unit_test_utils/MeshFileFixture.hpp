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

#ifndef UNITTEST_MESHFILEFIXTURE_HPP
#define UNITTEST_MESHFILEFIXTURE_HPP

#include <gtest/gtest.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace unit_test_util
{

class MeshFileFixture : public MeshFixture
{
protected:
    MeshFileFixture()
    : stkIo(get_comm())
    {
        allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    }
    void read_mesh(const std::string &fileToRead)
    {
        stk::io::fill_mesh_preexisting(stkIo, fileToRead, get_bulk());
    }
    void NGPTearDown()
    {
        unlink(filename.c_str());
    }
    const std::string filename = "filename.exo";
    stk::io::StkMeshIoBroker stkIo;
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
MeshFileFixture : public MeshFixture
{
protected:
    MeshFileFixture()
    : stkIo(get_comm())
    {
        allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    }
    void read_mesh(const std::string &fileToRead)
    {
        stk::io::fill_mesh_preexisting(stkIo, fileToRead, get_bulk());
    }
    void NGPTearDown()
    {
        unlink(filename.c_str());
    }
    const std::string filename = "filename.exo";
    stk::io::StkMeshIoBroker stkIo;
};

} // namespace simple_fields

}}

#endif

