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

#ifndef STK_STK_UNIT_TEST_UTILS_TransferPingPong_hpp_
#define STK_STK_UNIT_TEST_UTILS_TransferPingPong_hpp_

#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/FieldEvaluator.hpp>
#include <stk_transfer_util/TransferMainSettings.hpp>
#include <stk_transfer_util/TransferMainBroker.hpp>
#include <stk_search_util/MasterElementProvider.hpp>

#include <string>
#include <memory>

namespace stk {
namespace unit_test_util {

class TransferPingPong {
public:
  TransferPingPong(MPI_Comm comm,
                   const std::string& sendMeshName,
                   const std::string& sendFieldName,
                   const FieldEvaluator& sendFieldEval,
                   const std::string& recvMeshName,
                   const std::string& recvFieldName,
                   const FieldEvaluator& recvFieldEval,
                   std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider);
                   
  void run_steps(int numSteps);

  double get_max_err() const;

private:
  transfer_util::TransferMainSettings m_pingSettings;
  transfer_util::TransferMainSettings m_pongSettings;
  std::shared_ptr<transfer_util::TransferMainBroker> m_pingBroker;
  std::shared_ptr<transfer_util::TransferMainBroker> m_pongBroker;
  const FieldEvaluator& m_sendFieldEval;
  std::string m_errFieldName;
};

double run_ping_pong_transfer(MPI_Comm comm,
                              int numSteps,
                              const FieldEvaluator& sendFieldEval,
                              std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
                              const std::string& sendMeshName = "generated:3x3x3|bbox:0,0,0,1,1,1",
                              const std::string& recvMeshName = "generated:5x5x5|bbox:0,0,0,1,1,1");

}} // namespace stk::unit_test_util

#endif /* STK_STK_UNIT_TEST_UTILS_TransferPingPong_hpp_ */
