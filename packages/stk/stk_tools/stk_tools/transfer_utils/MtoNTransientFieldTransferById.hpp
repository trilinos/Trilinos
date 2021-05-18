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

#ifndef MTONTRANSIENTFIELDTRANSFERBYID_HPP
#define MTONTRANSIENTFIELDTRANSFERBYID_HPP

#include "TransientFieldTransferById.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_util/registry/ProductRegistry.hpp"
#include "stk_io/IOHelpers.hpp"
#include "Teuchos_RCP.hpp"                         // for RCP::RCP<T>, etc
#include "Teuchos_RCPDecl.hpp"                     // for RCP
#include <map>

namespace stk {
namespace transfer_utils {

class M2NOutputSerializerBulkData : public stk::mesh::BulkData
{
public:
  M2NOutputSerializerBulkData(stk::mesh::MetaData & mesh_meta_data, ParallelMachine parallel)
    : BulkData(mesh_meta_data, parallel, stk::mesh::BulkData::NO_AUTO_AURA, true)
  {}

  virtual ~M2NOutputSerializerBulkData() override = default;

  void switch_to_serial_mesh() {

    modification_begin();

    stk::mesh::EntityVector nodesToUnshare;
    stk::mesh::get_entities(*this, stk::topology::NODE_RANK, mesh_meta_data().globally_shared_part(), nodesToUnshare);

    for (const stk::mesh::Entity & node : nodesToUnshare) {
      internal_set_owner(node, 0);
      remove_entity_comm(node);
      entity_comm_map_clear(entity_key(node));
    }
    destroy_all_ghosting();

    stk::mesh::PartVector addParts{&mesh_meta_data().locally_owned_part()};
    stk::mesh::PartVector removeParts{&mesh_meta_data().globally_shared_part()};
    internal_verify_and_change_entity_parts(nodesToUnshare, addParts, removeParts);

    m_parallel = Parallel(MPI_COMM_SELF);

    modification_end();
  }
};

constexpr unsigned INVALID_SUBDOMAIN = std::numeric_limits<unsigned>::max();

inline bool is_valid_subdomain(unsigned subdomain) {
  return subdomain != INVALID_SUBDOMAIN;
}

class SubdomainWriterBase {
public:
  SubdomainWriterBase(const stk::io::StkMeshIoBroker & inputBroker, M2NOutputSerializerBulkData * bulk)
    : m_inputBroker(inputBroker),
      m_bulk(bulk)
  {}

  virtual ~SubdomainWriterBase() = default;

  stk::mesh::BulkData & get_bulk_data() { return *m_bulk; }
  stk::mesh::MetaData & get_meta_data() { return m_bulk->mesh_meta_data(); }

  virtual void setup_output_file(const std::string & fileName, unsigned subdomain, unsigned numSubdomains,
                                 int globalNumNodes, int globalNumElems) = 0;

  virtual void write_mesh() = 0;

  virtual void write_transient_data(double timeStep) = 0;

protected:
  const stk::io::StkMeshIoBroker & m_inputBroker;
  M2NOutputSerializerBulkData* m_bulk;
};


class SubdomainWriter : public SubdomainWriterBase {
public:
  SubdomainWriter(const stk::io::StkMeshIoBroker & inputBroker, M2NOutputSerializerBulkData * bulk,
                  const stk::io::EntitySharingInfo & nodeSharingInfo)
    : SubdomainWriterBase(inputBroker, bulk),
      m_nodeSharingInfo(nodeSharingInfo)
  {}

  virtual ~SubdomainWriter() override
  {
    stk::io::delete_selector_property(*m_outRegion);
    delete m_outRegion;
  }

  virtual void setup_output_file(const std::string & fileName, unsigned subdomain, unsigned numSubdomains,
                                 int globalNumNodes, int globalNumElems) override
  {
    Ioss::DatabaseIO *dbo = stk::io::create_database_for_subdomain(fileName, subdomain, numSubdomains);
    m_outRegion = new Ioss::Region(dbo, fileName);

    stk::io::add_properties_for_subdomain(*m_bulk, *m_outRegion, subdomain, numSubdomains, globalNumNodes, globalNumElems);

    int dbIntSize = m_inputBroker.check_integer_size_requirements_serial();
    if (dbIntSize > 4) {
      m_outRegion->property_add(Ioss::Property("INTEGER_SIZE_API", dbIntSize));
      m_outRegion->property_add(Ioss::Property("INTEGER_SIZE_DB", dbIntSize));
    }
  }

  void add_qa_records()
  {
    std::vector<stk::io::QaRecord> qaRecords = m_inputBroker.get_qa_records();

    for (const stk::io::QaRecord &qaRec : qaRecords) {
      m_outRegion->add_qa_record(qaRec.name, qaRec.version, qaRec.date, qaRec.time);
    }

    std::string codeName = "stk_balance_m2n";
    std::string codeVersion = stk::ProductRegistry::version();
    m_outRegion->property_add(Ioss::Property(std::string("code_name"), codeName));
    m_outRegion->property_add(Ioss::Property(std::string("code_version"), codeVersion));
  }

  void add_info_records()
  {
    std::vector<std::string> infoRecords = m_inputBroker.get_info_records();

    m_outRegion->add_information_records(infoRecords);
  }

  void add_global_variables()
  {
    std::vector<std::string> globalVariableNames;
    m_inputBroker.get_global_variable_names(globalVariableNames);
    Teuchos::RCP<Ioss::Region> rcpRegion;
    rcpRegion.reset(m_outRegion, false);

    for (const std::string& globalVariableName : globalVariableNames) {
      size_t length = m_inputBroker.get_global_variable_length(globalVariableName);
      stk::io::internal_add_global(rcpRegion, globalVariableName, length, Ioss::Field::DOUBLE);
    }
  }

  virtual void write_mesh() override
  {
    add_qa_records();
    add_info_records();
    stk::io::write_file_for_subdomain(*m_outRegion, *m_bulk, m_nodeSharingInfo);
    add_global_variables();
  }

  void write_global_variables(int step)
  {
    m_outRegion->begin_mode(Ioss::STATE_TRANSIENT);
    m_outRegion->begin_state(step);

    std::vector<std::string> globalVariableNames;
    m_inputBroker.get_global_variable_names(globalVariableNames);

    Teuchos::RCP<Ioss::Region> rcpRegion;
    rcpRegion.reset(m_outRegion, false);

    std::vector<double> globalVariable;
    for (const std::string& globalVariableName : globalVariableNames) {
      m_inputBroker.get_global(globalVariableName, globalVariable);
      stk::io::internal_write_global(rcpRegion, globalVariableName, globalVariable);
    }

    m_outRegion->end_state(step);
    m_outRegion->end_mode(Ioss::STATE_TRANSIENT);
  }

  virtual void write_transient_data(double time) override
  {
    const int step = stk::io::write_transient_data_for_subdomain(*m_outRegion, *m_bulk, time);

    write_global_variables(step);
  }

protected:
  Ioss::Region* m_outRegion;
  stk::io::EntitySharingInfo m_nodeSharingInfo;
};


class EmptySubdomainWriter : public SubdomainWriterBase {
public:
  EmptySubdomainWriter(const stk::io::StkMeshIoBroker & inputBroker, M2NOutputSerializerBulkData * bulk)
    : SubdomainWriterBase(inputBroker, bulk)
  {}

  virtual ~EmptySubdomainWriter() override = default;

  virtual void setup_output_file(const std::string & /*fileName*/, unsigned /*subdomain*/, unsigned /*numSubdomains*/,
                                 int /*globalNumNodes*/, int /*globalNumElems*/) override {}

  virtual void write_mesh() override {}

  void parallel_consistent_global_variable_read()
  {
    std::vector<std::string> globalVariableNames;
    m_inputBroker.get_global_variable_names(globalVariableNames);

    std::vector<double> globalVariable;
    for (const std::string& globalVariableName : globalVariableNames) {
      m_inputBroker.get_global(globalVariableName, globalVariable);
    }
  }

  virtual void write_transient_data(double /*time*/) override
  {
    parallel_consistent_global_variable_read();
  }
};

class MtoNTransientFieldTransferById
{
public:
    MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker &inputBroker, unsigned numSubDomain);

    MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker &inputBroker, unsigned numSubDomain,
                                   const std::vector<stk::mesh::EntityRank> &entityRanks);

    ~MtoNTransientFieldTransferById();

    void setup_subdomain(M2NOutputSerializerBulkData& bulk, const std::string &filename,
                         unsigned subdomain, const stk::io::EntitySharingInfo& nodeSharingInfo,
                         int global_num_nodes, int global_num_elems);

    void transfer_transient_data(unsigned subdomain);
    void transfer_and_write_transient_data(unsigned subdomain);

    SubdomainWriterBase & get_subdomain_writer(unsigned subdomain);

private:
    MtoNTransientFieldTransferById();

    void initialize_transfer(unsigned subDomainInfo);

    stk::io::StkMeshIoBroker &m_inputBroker;
    unsigned m_numSubDomain;
    std::map<unsigned, SubdomainWriterBase*> m_subdomainWriters;
    std::map<unsigned, std::vector<TransientTransferByIdForRank*>> m_subdomainTransfers;
    std::vector<stk::mesh::EntityRank> m_entityRanks;
};

}
}

#endif // MTONTRANSIENTFIELDTRANSFERBYID_HPP
