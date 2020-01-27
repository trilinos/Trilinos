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
#include "SearchById.hpp"
#include "TransferCopyByIdMeshAdapter.hpp"


namespace stk {
namespace transfer {


class TransferCopyById : public TransferBase {
public :

  typedef SearchById::Mesh_ID Mesh_ID;
  typedef SearchById::KeyToTargetProcessor KeyToTargetProcessor;
  typedef SearchById::MeshIDSet MeshIDSet;

  TransferCopyById(SearchById & search_in, TransferCopyByIdMeshAdapter & mesha_in, TransferCopyByIdMeshAdapter & meshb_in)
    :m_search(search_in)
    ,m_mesha(mesha_in)
    ,m_meshb(meshb_in)
  {}
  virtual ~TransferCopyById() {};
  virtual void coarse_search()
  {
    m_search.do_search(m_mesha,m_meshb,m_key_to_target_processor);
  }
  virtual void communication() {};
  virtual void local_search() {};
  virtual void apply()
  {
    do_transfer(m_key_to_target_processor,m_mesha,m_meshb);
  }

private:
  void do_transfer(const KeyToTargetProcessor & key_to_target_processor,
                   const TransferCopyByIdMeshAdapter & mesha,
                   TransferCopyByIdMeshAdapter & meshb);

  SearchById & m_search;
  TransferCopyByIdMeshAdapter & m_mesha;
  TransferCopyByIdMeshAdapter & m_meshb;
  KeyToTargetProcessor m_key_to_target_processor;
};

}  } // namespace transfer stk

#endif //  STK_COPYTRANSFER_HPP
