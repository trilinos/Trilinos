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

#ifndef TEXTMESHTOFILE_HPP
#define TEXTMESHTOFILE_HPP

#include <stddef.h>                             // for size_t
#include <stk_io/StkMeshIoBroker.hpp>           // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>           // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>           // for MetaData
#include <string>                               // for string
#include "stk_util/parallel/Parallel.hpp"       // for ParallelMachine
#include "stk_unit_test_utils/BuildMesh.hpp"

namespace stk
{
namespace unit_test_util
{

class TextMeshToFile
{
public:
  TextMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption);
  ~TextMeshToFile() = default;

  void setup_mesh(const std::string& meshDesc, const std::string& outputFileName);
  void write_mesh();

  stk::mesh::BulkData& get_bulk() { return m_bulk; }
  stk::mesh::MetaData& get_meta() { return m_meta; }

protected:
  std::shared_ptr<stk::mesh::BulkData> m_bulkPtr;
  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;
  stk::io::StkMeshIoBroker m_broker;
  size_t m_outputFileIndex = 0;
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
TextMeshToFile
{
public:
  TextMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption);
  ~TextMeshToFile() = default;

  void setup_mesh(const std::string& meshDesc, const std::string& outputFileName);
  void write_mesh();

  stk::mesh::BulkData& get_bulk() { return m_bulk; }
  stk::mesh::MetaData& get_meta() { return m_meta; }

protected:
  std::shared_ptr<stk::mesh::BulkData> m_bulkPtr;
  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;
  stk::io::StkMeshIoBroker m_broker;
  size_t m_outputFileIndex = 0;
};

} // namespace simple_fields

}
}
#endif // TEXTMESHTOFILE_HPP
