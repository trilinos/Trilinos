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

#include "SubdomainWriter.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_util/registry/ProductRegistry.hpp"

namespace stk {
namespace balance {

SubdomainWriterBase::SubdomainWriterBase(const stk::io::StkMeshIoBroker& inputBroker, OutputSerializerBulkData* bulk)
  : m_inputBroker(inputBroker),
    m_bulk(bulk)
{}


SubdomainWriter::SubdomainWriter(const stk::io::StkMeshIoBroker& inputBroker, OutputSerializerBulkData* bulk,
                                 const stk::io::EntitySharingInfo& nodeSharingInfo)
  : SubdomainWriterBase(inputBroker, bulk),
    m_nodeSharingInfo(nodeSharingInfo)
{}

SubdomainWriter::~SubdomainWriter()
{
  stk::io::delete_selector_property(*m_outRegion);
  delete m_outRegion;
}

void
SubdomainWriter::setup_output_file(const std::string& fileName, unsigned subdomain, unsigned numSubdomains,
                                   int globalNumNodes, int globalNumElems)
{
  bool use64Bit = false;
  int dbIntSize = m_inputBroker.check_integer_size_requirements_serial();
  if (dbIntSize > 4) {
    use64Bit = true;
  }

  Ioss::DatabaseIO *dbo = stk::io::create_database_for_subdomain(fileName, subdomain, numSubdomains, use64Bit);

  const size_t inputFileIndex = m_inputBroker.get_active_mesh();
  const Ioss::DatabaseIO* inputDBIO = m_inputBroker.get_input_database(inputFileIndex);
  const int maxSymbolLength = inputDBIO->maximum_symbol_length();
  if (maxSymbolLength > 0) {
    dbo->set_maximum_symbol_length(maxSymbolLength);
  }

  m_outRegion = new Ioss::Region(dbo, fileName);
  stk::io::OutputParams params(*m_outRegion, *m_bulk);
  stk::io::add_properties_for_subdomain(params, subdomain, numSubdomains, globalNumNodes, globalNumElems);

  if (use64Bit) {
    m_outRegion->property_add(Ioss::Property("INTEGER_SIZE_API", dbIntSize));
    m_outRegion->property_add(Ioss::Property("INTEGER_SIZE_DB", dbIntSize));
  }
}

void
SubdomainWriter::add_qa_records()
{
  std::vector<stk::io::QaRecord> qaRecords = m_inputBroker.get_qa_records();

  for (const stk::io::QaRecord& qaRec : qaRecords) {
    m_outRegion->add_qa_record(qaRec.name, qaRec.version, qaRec.date, qaRec.time);
  }

  std::string codeName = "stk_balance";
  std::string codeVersion = stk::ProductRegistry::version();
  m_outRegion->property_add(Ioss::Property(std::string("code_name"), codeName));
  m_outRegion->property_add(Ioss::Property(std::string("code_version"), codeVersion));
}

void
SubdomainWriter::add_info_records()
{
  std::vector<std::string> infoRecords = m_inputBroker.get_info_records();

  m_outRegion->add_information_records(infoRecords);
}

void
SubdomainWriter::add_global_variables()
{
  std::vector<std::string> globalVariableNames;
  m_inputBroker.get_global_variable_names(globalVariableNames);
  
  std::shared_ptr<Ioss::Region> region(m_outRegion, [](auto /*pointerWeWontDelete*/){});

  for (const std::string& globalVariableName : globalVariableNames) {
    size_t length = m_inputBroker.get_global_variable_length(globalVariableName);
    stk::io::impl::add_global(region, globalVariableName, length, Ioss::Field::DOUBLE);
  }
  
}

void
SubdomainWriter::write_mesh()
{
  add_qa_records();
  add_info_records();
  stk::io::OutputParams params(*m_outRegion, *m_bulk);
  stk::io::write_file_for_subdomain(params, m_nodeSharingInfo);
  add_global_variables();
}

void
SubdomainWriter::write_global_variables(int step)
{
  m_outRegion->begin_mode(Ioss::STATE_TRANSIENT);
  m_outRegion->begin_state(step);

  std::vector<std::string> globalVariableNames;
  m_inputBroker.get_global_variable_names(globalVariableNames);

  std::shared_ptr<Ioss::Region> region(m_outRegion, [](auto /*pointerWeWontDelete*/){});

  std::vector<double> globalVariable;
  for (const std::string& globalVariableName : globalVariableNames) {
    m_inputBroker.get_global(globalVariableName, globalVariable);
    stk::io::impl::write_global(region, globalVariableName, globalVariable);
  }

  m_outRegion->end_state(step);
  m_outRegion->end_mode(Ioss::STATE_TRANSIENT);
}

void
SubdomainWriter::write_transient_data(double time)
{
  stk::io::OutputParams params(*m_outRegion, *m_bulk);
  const int step = stk::io::write_transient_data_for_subdomain(params, time);

  write_global_variables(step);
}


EmptySubdomainWriter::EmptySubdomainWriter(const stk::io::StkMeshIoBroker& inputBroker, OutputSerializerBulkData* bulk)
  : SubdomainWriterBase(inputBroker, bulk)
{}

void
EmptySubdomainWriter::setup_output_file(const std::string& , unsigned , unsigned , int , int )
{}

void
EmptySubdomainWriter::write_mesh()
{}

void
EmptySubdomainWriter::parallel_consistent_global_variable_read()
{
  std::vector<std::string> globalVariableNames;
  m_inputBroker.get_global_variable_names(globalVariableNames);

  std::vector<double> globalVariable;
  for (const std::string& globalVariableName : globalVariableNames) {
    m_inputBroker.get_global(globalVariableName, globalVariable);
  }
}

void
EmptySubdomainWriter::write_transient_data(double )
{
  parallel_consistent_global_variable_read();
}

}
}
