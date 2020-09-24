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
#ifndef _DisconnectBlocks_hpp_
#define _DisconnectBlocks_hpp_
#include "DisconnectBlocksImpl.hpp"
#include "DisconnectTypes.hpp"

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace tools {

enum DisconnectOption {
  DISCONNECT_GLOBAL,
  DISCONNECT_LOCAL
};

enum SnipOption {
  PRESERVE_INITIAL_HINGES,
  SNIP_ALL_HINGES
};

struct DisconnectBlocksOption {
  DisconnectBlocksOption()
   : disconnectOption(DISCONNECT_GLOBAL)
   , snipOption(PRESERVE_INITIAL_HINGES)
  {}

  DisconnectBlocksOption(DisconnectOption disconnectOption_, SnipOption snipOption_)
   : disconnectOption(disconnectOption_)
   , snipOption(snipOption_)
  {}

  DisconnectOption disconnectOption;
  SnipOption snipOption;
};

void disconnect_all_blocks(stk::mesh::BulkData& bulk, bool preserveOrphans = false);
void disconnect_all_blocks(stk::mesh::BulkData & bulk, impl::LinkInfo& info, bool preserveOrphans = false);
void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockPairVector& blockPairsToDisconnect,
                            DisconnectBlocksOption options = DisconnectBlocksOption());
void disconnect_user_blocks(stk::mesh::BulkData& bulk, const BlockNamePairVector& blockNamePairsToDisconnect,
                            DisconnectBlocksOption options = DisconnectBlocksOption());
}

}

#endif
