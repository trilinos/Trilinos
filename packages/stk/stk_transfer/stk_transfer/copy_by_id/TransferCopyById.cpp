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
#include <cstring>
#include <cstddef>

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

  STK_ThrowRequireMsg(m_numFields == m_mesha.num_fields(), "Send and receive regions have different number of fields");

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

void byte_copy(const int bytesPerScalar, const int scalarByteStride,
               const void * f_dataA, const int numSrcBytes,
               void * f_dataB, const int numRcvBytes)
{
  const std::byte* srcValueBytes = reinterpret_cast<const std::byte*>(f_dataA);
  std::byte* rcvValueBytes = reinterpret_cast<std::byte*>(f_dataB);

  int numBytes = std::min(numSrcBytes, numRcvBytes);
  while (numBytes) {
    for (int scalarByteIdx = 0; scalarByteIdx < bytesPerScalar; ++scalarByteIdx) {
      rcvValueBytes[scalarByteIdx] = srcValueBytes[scalarByteIdx];
    }
    numBytes -= bytesPerScalar;
    srcValueBytes += bytesPerScalar;
    rcvValueBytes += scalarByteStride;
  }
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
        STK_ThrowRequireMsg(this_field_size_a == this_field_size_b, 
                        "field_size_a " << this_field_size_a << " field_size_b " << this_field_size_b);
      }

      int numSrcBytes = m_mesha.num_bytes(key,f);
      int numRcvBytes = m_meshb.num_bytes(key,f);
      const int bytesPerScalar = m_meshb.bytes_per_scalar(key,f);
      const int scalarByteStride = m_meshb.scalar_byte_stride(key,f);

      byte_copy(bytesPerScalar, scalarByteStride, f_dataA, numSrcBytes, f_dataB, numRcvBytes);
    } else {
      DataTypeKey::data_t sentDataType = m_sendFieldDataTypes[f];

      int numSrcComponents = m_mesha.num_components(key,f);
      int srcComponentStride = m_mesha.component_stride(key,f);

      int numRcvComponents = m_meshb.num_components(key,f);
      int rcvComponentStride = m_meshb.component_stride(key,f);

      m_dataTranslators[sentDataType]->translate(f_dataA, numSrcComponents, srcComponentStride, m_recvFieldDataTypes[f],
                                                 f_dataB, numRcvComponents, rcvComponentStride);
    }
  }
}

void TransferCopyById::pack_fields(const int target_proc, const Mesh_ID key, CommSparse& commSparse)
{
  commSparse.send_buffer(target_proc).pack<Mesh_ID>(key);

  for (unsigned f=0; f<m_numFields; ++f)  {
    unsigned numBytes = m_mesha.num_bytes(key,f);
    uint8_t * f_data = reinterpret_cast<uint8_t *>(m_mesha.field_data(key,f));

    unsigned sendFieldSize = 0;

    if(m_fieldCompatibility[f]) {
      sendFieldSize = numBytes;
    } else {
      DataTypeKey dataTranslator(m_sendFieldDataTypes[f], numBytes);
      sendFieldSize = dataTranslator.m_value;
    }
    commSparse.send_buffer(target_proc).pack<unsigned>(sendFieldSize);

    const int bytesPerScalar = m_mesha.bytes_per_scalar(key,f);

    while (numBytes) {
      for (int scalarByteIdx = 0; scalarByteIdx < bytesPerScalar; ++scalarByteIdx) {
        commSparse.send_buffer(target_proc).pack<uint8_t>(f_data[scalarByteIdx]);
      }
      numBytes -= bytesPerScalar;
      f_data += bytesPerScalar;
    }
  }
}

void TransferCopyById::initialize_commsparse_buffers()
{
  m_mesha.acquire_field_data();
  m_meshb.acquire_field_data();

  pack_commsparse();
  m_commSparse.allocate_buffers();

  m_mesha.release_field_data();
  m_meshb.release_field_data();
}

void TransferCopyById::pack_commsparse()
{
  const int myProc = m_commSparse.parallel_rank();

  for (const auto & k2tEntry : m_key_to_target_processor) {
    const int targetProc = k2tEntry.second;
    if (targetProc != myProc) {
      const Mesh_ID key = k2tEntry.first;
      pack_fields(targetProc, key, m_commSparse);
    }
  }
}

void TransferCopyById::pack_send_fields(CommSparse& commSparse)
{
  const ParallelMachine comm = m_mesha.comm();
  const int my_proc = parallel_machine_rank(comm);

  for (const auto & k2tEntry : m_key_to_target_processor) {
    const Mesh_ID key = k2tEntry.first;
    const int target_proc = k2tEntry.second;

    if (target_proc == my_proc) {
      local_copy(key);
      continue;
    }
    pack_fields(target_proc, key, commSparse);
  }
}

void TransferCopyById::unpack_and_copy_fields(CommBuffer& recvBuffer, CopyTransferUnpackInfo& unpackInfo)
{
  unsigned f = unpackInfo.fieldIndex;
  Mesh_ID key = unpackInfo.key;

  unsigned fsize = std::min(m_meshb.field_data_size(key,f), unpackInfo.sentDataTypeKey);
  unpackInfo.buffer.resize(unpackInfo.sentDataTypeKey);

  for (unsigned index = 0 ; index < unpackInfo.sentDataTypeKey; ++index) {
    recvBuffer.unpack<uint8_t>(unpackInfo.buffer[index]);
  }

  const int bytesPerScalar = m_meshb.bytes_per_scalar(key,f);
  const int scalarByteStride = m_meshb.scalar_byte_stride(key,f);

  byte_copy(bytesPerScalar, scalarByteStride, unpackInfo.buffer.data(), unpackInfo.sentDataTypeKey, unpackInfo.fieldData, fsize);
}

void TransferCopyById::unpack_and_copy_fields_with_compatibility(CommBuffer& recvBuffer, CopyTransferUnpackInfo& unpackInfo)
{
  unsigned f = unpackInfo.fieldIndex;
  Mesh_ID key = unpackInfo.key;

  DataTypeKey unpackedData(unpackInfo.sentDataTypeKey);
  DataTypeKey::data_t sentDataType = unpackedData.get_data_type();
  unsigned remoteFieldSize = unpackedData.get_data_length();

  unpackInfo.buffer.resize(remoteFieldSize);

  for(unsigned index = 0 ; index < remoteFieldSize; ++index) {
    recvBuffer.unpack<uint8_t>(unpackInfo.buffer[index]);
  }

  int numSrcComponents = m_dataTranslators[sentDataType]->get_src_component_size(remoteFieldSize);
  int srcComponentStride = 1;

  int numRcvComponents = m_meshb.num_components(key,f);
  int rcvComponentStride = m_meshb.component_stride(key,f);

  m_dataTranslators[sentDataType]->translate(unpackInfo.buffer.data(), numSrcComponents, srcComponentStride, m_recvFieldDataTypes[unpackInfo.fieldIndex],
                                             unpackInfo.fieldData, numRcvComponents, rcvComponentStride);
}

void TransferCopyById::unpack_fields(MeshIDSet & remoteKeys, const int recv_proc, CommBuffer& recvBuffer)
{
  Mesh_ID key;
  std::vector<uint8_t> tmpBuffer;
  recvBuffer.unpack<Mesh_ID>(key);

  if(remoteKeys.find(key) == remoteKeys.end()) {
    ++m_errorCount;
    const int my_proc = parallel_machine_rank(m_mesha.comm());
    m_errorMsg << "P" << my_proc << " Error, proc = " << recv_proc << " sent unrequested data for key = " << key << std::endl;
  } else {
    remoteKeys.erase(key);
  }
  for(unsigned f=0 ; f<m_numFields ; ++f) {
    uint8_t * f_data = reinterpret_cast<uint8_t*>(m_meshb.field_data(key,f));
    unsigned sentDataTypeKey = 0;

    recvBuffer.unpack<unsigned>(sentDataTypeKey);

    CopyTransferUnpackInfo unpackInfo(tmpBuffer, f_data, sentDataTypeKey, key, f);

    if(m_fieldCompatibility[f]) {
      unpack_and_copy_fields(recvBuffer, unpackInfo);
    } else {
      unpack_and_copy_fields_with_compatibility(recvBuffer, unpackInfo); 
    }
  }
}

void TransferCopyById::receive_fields(MeshIDSet & remoteKeys)
{
  for(int fromProc=0; fromProc<m_commSparse.parallel_size(); ++fromProc) {
    stk::CommBuffer& recvBuffer = m_commSparse.recv_buffer(fromProc);
    while(recvBuffer.remaining()) {
      unpack_fields(remoteKeys, fromProc, recvBuffer);
    }
  }
}

void TransferCopyById::check_received_keys(MeshIDSet & remoteKeys)
{
  const ParallelMachine comm = m_mesha.comm();
  const int my_proc = parallel_machine_rank(comm);

  if(!remoteKeys.empty()) {
    m_errorMsg << "P" << my_proc << " Error, Did not receive all the requested keys from other processors, unfulfilled keys = [";
    for(MeshIDSet::const_iterator set_it=remoteKeys.begin() ; set_it != remoteKeys.end() ; ++set_it) {
      m_errorMsg << m_meshb.print_mesh_id(*set_it) << " ";
    }
    m_errorMsg << "]" << std::endl;
    ++m_errorCount;
  }
  all_reduce( comm , ReduceSum<1>( & m_errorCount ) );
  if(m_errorCount) {
    all_write_string( comm , std::cerr , m_errorMsg.str() );
    STK_ThrowErrorMsg("Error in communication during CopyTransfer!\n");
  }
}

void TransferCopyById::do_transfer()
{
  m_errorCount = 0;
  m_errorMsg.str("");
  m_errorMsg.clear();

  m_mesha.begin_transfer();
  m_meshb.begin_transfer();

  pack_send_fields(m_commSparse);
  const bool deallocateSendBuffers = false;
  m_commSparse.communicate(deallocateSendBuffers);

  MeshIDSet remoteKeys = m_search.get_remote_keys();

  receive_fields(remoteKeys);
  check_received_keys(remoteKeys); 

  m_mesha.end_transfer();
  m_meshb.end_transfer();

  m_commSparse.reset_buffers();
}

}  } // namespace transfer stk
