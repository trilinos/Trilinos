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

#include "TextMeshToFile.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

namespace stk
{
namespace unit_test_util
{

TextMeshToFile::TextMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  : m_bulkPtr(build_mesh(3, comm, auraOption)),
    m_bulk(*m_bulkPtr),
    m_meta(m_bulk.mesh_meta_data())
{
}

void
TextMeshToFile::write_mesh()
{
  m_broker.write_output_mesh(m_outputFileIndex);
}

void TextMeshToFile::setup_mesh(const std::string& meshDesc, const std::string& outputFileName)
{
  stk::unit_test_util::setup_text_mesh(m_bulk, meshDesc);

  m_broker.set_bulk_data(m_bulk);
  m_broker.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
  m_broker.property_add(Ioss::Property("INTEGER_SIZE_DB", 8));
  m_outputFileIndex = m_broker.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
}

namespace simple_fields {

TextMeshToFile::TextMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  : m_bulkPtr(build_mesh(3, comm, auraOption)),
    m_bulk(*m_bulkPtr),
    m_meta(m_bulk.mesh_meta_data())
{
}

void
TextMeshToFile::write_mesh()
{
  m_broker.write_output_mesh(m_outputFileIndex);
}

void TextMeshToFile::setup_mesh(const std::string& meshDesc, const std::string& outputFileName)
{
  stk::unit_test_util::setup_text_mesh(m_bulk, meshDesc);

  m_broker.set_bulk_data(m_bulk);
  m_broker.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
  m_broker.property_add(Ioss::Property("INTEGER_SIZE_DB", 8));
  m_outputFileIndex = m_broker.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
}

} // namespace simple_fields

}
}
