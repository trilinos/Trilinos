// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test to exercise fix to Teuchos::Comm's createSubcommunicator, which
// had leaked memory.
// The fix added MPI_Comm_free to the opaqueWrappers.
// 8/2018  This test hangs of platform waterman.

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Comm, MessageTags, PacketType )
{
  auto comm = Teuchos::rcp_const_cast<Teuchos::Comm<int>>(Teuchos::DefaultComm<int>::getComm ());

  auto tag0 = comm->incrementTag();
  auto tag1 = comm->incrementTag();
  TEST_EQUALITY(tag0+1, tag1);

  std::vector<int> sub_comm_ranks = {0};
  auto sub_comm = comm->createSubcommunicator(sub_comm_ranks);

  if (!sub_comm.is_null()) {
    auto sub_comm_tag0 = sub_comm->incrementTag();
    auto sub_comm_tag1 = sub_comm->incrementTag();
    TEST_EQUALITY(sub_comm_tag0+1, sub_comm_tag1);
  }

  auto tag2 = comm->incrementTag();
  TEST_EQUALITY(tag0+2, tag2);

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Comm, MessageTags, int )
