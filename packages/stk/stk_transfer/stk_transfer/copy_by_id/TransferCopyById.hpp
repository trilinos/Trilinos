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


#ifndef  STK_COPYTRANSFER_HPP
#define  STK_COPYTRANSFER_HPP

#include <stk_transfer/TransferBase.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include "SearchById.hpp"
#include "TransferCopyByIdMeshAdapter.hpp"

namespace stk { class CommSparse; }

namespace stk {
namespace transfer {

struct CopyTransferUnpackInfo {
  std::vector<uint8_t>& buffer;
  uint8_t* fieldData;
  unsigned sentDataTypeKey;
  unsigned fieldSize;
  unsigned fieldIndex;

  CopyTransferUnpackInfo(std::vector<uint8_t>& buffer_, uint8_t* fieldData_, unsigned sentDataTypeKey_, unsigned fieldSize_, unsigned fieldIndex_)
    : buffer(buffer_), fieldData(fieldData_), sentDataTypeKey(sentDataTypeKey_), fieldSize(fieldSize_), fieldIndex(fieldIndex_) {}
};

class TransferCopyById : public TransferBase {
public :

  typedef SearchById::Mesh_ID Mesh_ID;
  typedef SearchById::KeyToTargetProcessor KeyToTargetProcessor;
  typedef SearchById::MeshIDSet MeshIDSet;

  TransferCopyById(SearchById & search_in, TransferCopyByIdMeshAdapter & mesha_in, TransferCopyByIdMeshAdapter & meshb_in)
    :m_search(search_in)
    ,m_mesha(mesha_in)
    ,m_meshb(meshb_in)
    ,m_commSparse(mesha_in.comm())
  {
    setup_translators();
  }

  virtual ~TransferCopyById() {}
  virtual void coarse_search() override
  {
    m_search.do_search(m_mesha,m_meshb,m_key_to_target_processor);
    initialize_commsparse_buffers();
  }
  virtual void communication() override {}
  virtual void local_search() override {}
  virtual void apply() override { do_transfer(); }

  void setup_translators();

  static const stk::transfer::DataTypeTranslator<unsigned> translateUnsigned;
  static const stk::transfer::DataTypeTranslator<int64_t> translateInt64;
  static const stk::transfer::DataTypeTranslator<uint64_t> translateUInt64;
  static const stk::transfer::DataTypeTranslator<int> translateInt;
  static const stk::transfer::DataTypeTranslator<long double> translateLongDouble;
  static const stk::transfer::DataTypeTranslator<double> translateDouble;

private:
  void initialize_commsparse_buffers();
  void pack_commsparse();
  void do_transfer();
  void pack_send_fields(CommSparse& commSparse);
  void receive_fields(MeshIDSet & remoteKeys);
  void check_received_keys(MeshIDSet & remoteKeys);
  void local_copy(const Mesh_ID key);
  void pack_fields(const int target_proc, const Mesh_ID key, CommSparse& commSparse);
  void unpack_fields(MeshIDSet & remoteKeys, const int recv_proc, CommBuffer& recvBuffer);
  void unpack_and_copy_fields(CommBuffer& recvBuffer, CopyTransferUnpackInfo& unpackInfo);
  void unpack_and_copy_fields_with_compatibility(CommBuffer& recvBuffer, CopyTransferUnpackInfo& unpackInfo);

  SearchById & m_search;
  TransferCopyByIdMeshAdapter & m_mesha;
  TransferCopyByIdMeshAdapter & m_meshb;
  CommSparse m_commSparse;
  KeyToTargetProcessor m_key_to_target_processor;
  unsigned m_numFields;
  unsigned m_errorCount;
  std::ostringstream m_errorMsg;

  std::vector<DataTypeKey::data_t> m_sendFieldDataTypes;
  std::vector<DataTypeKey::data_t> m_recvFieldDataTypes;
  std::vector<int> m_fieldCompatibility;
  std::vector<const stk::transfer::TranslatorBase*> m_dataTranslators;
};

}  } // namespace transfer stk

#endif //  STK_COPYTRANSFER_HPP
