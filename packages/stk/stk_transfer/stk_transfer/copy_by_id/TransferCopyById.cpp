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
//

#include "TransferCopyById.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/diag/SlibDiagWriter.hpp>
#include <sstream>

namespace stk {
namespace transfer {

const DataTypeTranslator<unsigned> TransferCopyById::translateUnsigned;
const DataTypeTranslator<int64_t> TransferCopyById::translateInt64;
const DataTypeTranslator<uint64_t> TransferCopyById::translateUInt64;
const DataTypeTranslator<int> TransferCopyById::translateInt;
const DataTypeTranslator<long double> TransferCopyById::translateLongDouble;
const DataTypeTranslator<double> TransferCopyById::translateDouble;

void TransferCopyById::setup_translators()
{
  m_numFields = m_meshb.num_fields();

  ThrowRequireMsg(m_numFields == m_mesha.num_fields(), "Send and receive regions have different number of fields");

  m_sendFieldDataTypes.resize(m_numFields);
  m_recvFieldDataTypes.resize(m_numFields);
  m_fieldCompatibility.resize(m_numFields);

  for(unsigned i = 0; i < m_numFields; i++) {
    m_sendFieldDataTypes[i] = m_mesha.get_field_type(i);
    m_recvFieldDataTypes[i] = m_meshb.get_field_type(i);
    m_fieldCompatibility[i] = (m_sendFieldDataTypes[i] == m_recvFieldDataTypes[i]);
  }
  
  m_dataTranslators = { &translateUnsigned, 
                        &translateInt64,
                        &translateUInt64,
                        &translateInt,
                        &translateLongDouble,
                        &translateDouble };
}

void TransferCopyById::local_copy(const Mesh_ID key)
{
  for (unsigned f=0 ; f<m_numFields ; ++f) {
    const void * f_dataA = m_mesha.field_data(key,f);
    void * f_dataB = m_meshb.field_data(key,f);
    const unsigned this_field_size_a = m_mesha.field_data_size(key,f);
    const unsigned this_field_size_b = m_meshb.field_data_size(key,f);

    if(m_fieldCompatibility[f]) {
      unsigned fsize = std::min(this_field_size_a, this_field_size_b);
      if(fsize != 0) {
        ThrowRequireMsg(this_field_size_a == this_field_size_b, 
                        "field_size_a " << this_field_size_a << " field_size_b " << this_field_size_b);
      }
      std::memcpy(f_dataB, f_dataA, fsize);
    } else {
      DataTypeKey::data_t sentDataType = m_sendFieldDataTypes[f];
      m_dataTranslators[sentDataType]->translate(f_dataA, this_field_size_a, m_recvFieldDataTypes[f], f_dataB, this_field_size_b);
    }
  }
}

void TransferCopyById::pack_fields(const int target_proc, const Mesh_ID key, CommSparse& commSparse)
{
  commSparse.send_buffer(target_proc).pack<Mesh_ID>(key);

  for (unsigned f=0; f<m_numFields; ++f)  {
    const unsigned this_field_size = m_mesha.field_data_size(key,f);
    const uint8_t * f_data = reinterpret_cast<const uint8_t *>(m_mesha.field_data(key,f));

    unsigned sendFieldSize = 0;

    if(m_fieldCompatibility[f]) {
      sendFieldSize = this_field_size;
    } else {
      DataTypeKey dataTranslator(m_sendFieldDataTypes[f], this_field_size);
      sendFieldSize = dataTranslator.m_value;
    }
    commSparse.send_buffer(target_proc).pack<unsigned>(sendFieldSize);

    for (unsigned index=0 ; index<this_field_size ; ++index) {
      commSparse.send_buffer(target_proc).pack<uint8_t>(f_data[index]);
    }
  }
}

void TransferCopyById::send_fields(CommSparse& commSparse)
{
  const ParallelMachine comm = m_mesha.comm();
  const int my_proc = parallel_machine_rank(comm);

  for (int phase=0;phase<2;++phase) {
    KeyToTargetProcessor::const_iterator map_it = m_key_to_target_processor.begin();

    for ( ; map_it!=m_key_to_target_processor.end() ;  ++map_it) {
      const Mesh_ID key = map_it->first;
      const int target_proc = map_it->second;

      if (target_proc == my_proc) {
        if (0 == phase) {
          local_copy(key);
        }
        continue;
      }
      pack_fields(target_proc, key, commSparse);
    }

    if (phase == 0) {
      commSparse.allocate_buffers();
    } else {
      commSparse.communicate();
    }
  }
}

void TransferCopyById::unpack_and_copy_fields(const int recv_proc, CommSparse& commSparse, CopyTransferUnpackInfo& unpackInfo)
{
  unsigned fsize = std::min(unpackInfo.fieldSize, unpackInfo.sentDataTypeKey);
  unpackInfo.buffer.resize(unpackInfo.sentDataTypeKey);

  for (unsigned index = 0 ; index < unpackInfo.sentDataTypeKey; ++index) {
    commSparse.recv_buffer(recv_proc).unpack<uint8_t>(unpackInfo.buffer[index]);
  }
  std::memcpy(unpackInfo.fieldData, unpackInfo.buffer.data(), fsize);
}

void TransferCopyById::unpack_and_copy_fields_with_compatibility(const int recv_proc, CommSparse& commSparse, CopyTransferUnpackInfo& unpackInfo)
{
  DataTypeKey unpackedData(unpackInfo.sentDataTypeKey);
  DataTypeKey::data_t sentDataType = unpackedData.get_data_type();
  unsigned remoteFieldSize = unpackedData.get_data_length();

  unpackInfo.buffer.resize(remoteFieldSize);

  for(unsigned index = 0 ; index < remoteFieldSize; ++index) {
    commSparse.recv_buffer(recv_proc).unpack<uint8_t>(unpackInfo.buffer[index]);
  }
  m_dataTranslators[sentDataType]->translate(unpackInfo.buffer.data(), remoteFieldSize, m_recvFieldDataTypes[unpackInfo.fieldIndex], 
                                             unpackInfo.fieldData, unpackInfo.fieldSize);
}

void TransferCopyById::unpack_fields(MeshIDSet& remote_keys, const int recv_proc, CommSparse& commSparse)
{
  Mesh_ID key;
  std::vector<uint8_t> tmpBuffer;
  const int my_proc = parallel_machine_rank(m_mesha.comm());
  commSparse.recv_buffer(recv_proc).unpack<Mesh_ID>(key);

  if(remote_keys.find(key) == remote_keys.end()) {
    ++m_errorCount;
    m_errorMsg << "P" << my_proc << " Error, proc = " << recv_proc << " sent unrequested data for key = " << key << std::endl;
  } else {
    remote_keys.erase(key);
  }
  for(unsigned f=0 ; f<m_numFields ; ++f) {
    uint8_t * f_data = reinterpret_cast<uint8_t*>(m_meshb.field_data(key,f));
    unsigned sentDataTypeKey = 0;

    commSparse.recv_buffer(recv_proc).unpack<unsigned>(sentDataTypeKey);

    CopyTransferUnpackInfo unpackInfo(tmpBuffer, f_data, sentDataTypeKey, m_meshb.field_data_size(key,f), f);

    if(m_fieldCompatibility[f]) {
      unpack_and_copy_fields(recv_proc, commSparse, unpackInfo);
    } else {
      unpack_and_copy_fields_with_compatibility(recv_proc, commSparse, unpackInfo); 
    }
  }
}

void TransferCopyById::receive_fields(CommSparse& commSparse)
{
  const ParallelMachine comm = m_mesha.comm();
  const int my_proc = parallel_machine_rank(comm);
  const int num_proc = parallel_machine_size(comm);

  MeshIDSet remote_keys = m_search.get_remote_keys();

  for(int recv_proc=0;recv_proc<num_proc;++recv_proc) {
    if( my_proc != recv_proc ) {
      while(commSparse.recv_buffer(recv_proc).remaining()) {
        unpack_fields(remote_keys, recv_proc, commSparse);
      }
    }
  }
  if(!remote_keys.empty()) {
    m_errorMsg << "P" << my_proc << " Error, Did not receive all the requested keys from other processors, unfulfilled keys = [";
    for(MeshIDSet::const_iterator set_it=remote_keys.begin() ; set_it != remote_keys.end() ; ++set_it) {
      m_errorMsg << m_meshb.print_mesh_id(*set_it) << " ";
    }
    m_errorMsg << "]" << std::endl;
    ++m_errorCount;
  }
  all_reduce( comm , ReduceSum<1>( & m_errorCount ) );
  if(m_errorCount) {
    all_write_string( comm , std::cerr , m_errorMsg.str() );
    ThrowErrorMsg("Error in communication during CopyTransfer!\n");
  }
}

void TransferCopyById::do_transfer()
{
  stk::CommSparse commSparse(m_mesha.comm());
  m_errorCount = 0;
  m_errorMsg.str("");
  m_errorMsg.clear();

  send_fields(commSparse);
  receive_fields(commSparse); 
}

}  } // namespace transfer stk
