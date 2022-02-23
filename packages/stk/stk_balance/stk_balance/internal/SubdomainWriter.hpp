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

#ifndef SUBDOMAINWRITER_HPP
#define SUBDOMAINWRITER_HPP

#include "stk_balance/internal/OutputSerializerBulkData.hpp"
#include "stk_io/IOHelpers.hpp"

namespace Ioss { class Region; }

namespace stk {
namespace balance {

class SubdomainWriterBase {
public:
  SubdomainWriterBase(const stk::io::StkMeshIoBroker& inputBroker, OutputSerializerBulkData* bulk);
  virtual ~SubdomainWriterBase() = default;

  stk::mesh::BulkData& get_bulk_data() { return *m_bulk; }
  stk::mesh::MetaData& get_meta_data() { return m_bulk->mesh_meta_data(); }

  virtual void setup_output_file(const std::string& fileName, unsigned subdomain, unsigned numSubdomains,
                                 int globalNumNodes, int globalNumElems) = 0;
  virtual void write_mesh() = 0;
  virtual void write_transient_data(double timeStep) = 0;

protected:
  const stk::io::StkMeshIoBroker& m_inputBroker;
  OutputSerializerBulkData* m_bulk;
};


class SubdomainWriter : public SubdomainWriterBase {
public:
  SubdomainWriter(const stk::io::StkMeshIoBroker& inputBroker, OutputSerializerBulkData* bulk,
                  const stk::io::EntitySharingInfo& nodeSharingInfo);
  virtual ~SubdomainWriter() override;

  virtual void setup_output_file(const std::string & fileName, unsigned subdomain, unsigned numSubdomains,
                                 int globalNumNodes, int globalNumElems) override;
  virtual void write_mesh() override;
  virtual void write_transient_data(double time) override;

  void add_qa_records();
  void add_info_records();
  void add_global_variables();
  void write_global_variables(int step);

protected:
  Ioss::Region* m_outRegion;
  stk::io::EntitySharingInfo m_nodeSharingInfo;
};


class EmptySubdomainWriter : public SubdomainWriterBase {
public:
  EmptySubdomainWriter(const stk::io::StkMeshIoBroker& inputBroker, OutputSerializerBulkData* bulk);
  virtual ~EmptySubdomainWriter() override = default;

  virtual void setup_output_file(const std::string& fileName, unsigned subdomain, unsigned numSubdomains,
                                 int globalNumNodes, int globalNumElems) override;
  virtual void write_mesh() override;
  virtual void write_transient_data(double time) override;

  void parallel_consistent_global_variable_read();
};

}
}

#endif // SUBDOMAINWRITER_HPP
