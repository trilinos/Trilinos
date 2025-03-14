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

#include "stk_util/stk_config.h"            // for STK_HAS_MPI
#ifdef STK_HAS_MPI

#include "gtest/gtest.h"
#include "stk_util/parallel/MPICommKey.hpp"

namespace {

class KeyManagerTester : public ::testing::Test
{
  protected:
    KeyManagerTester()
    {
      comm1 = MPI_COMM_WORLD;
      MPI_Comm_dup(comm1, &comm2);
    }

    ~KeyManagerTester()
    {
      MPI_Comm_free(&comm2);
    }

    MPI_Comm comm1;
    MPI_Comm comm2;
    stk::impl::MPIKeyManager keyManager;
};

class CallbackChecker
{
  public:
    int callcount = 0;

    void set_called(MPI_Comm) { callcount++; }
};

}  // namespace

#ifndef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
namespace {
  MPI_Request req;
  bool is_complete = false;
  int comm_destructor_callback(MPI_Comm, int, void*, void*)
  {
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    is_complete = true;

    return MPI_SUCCESS;
  }
}

TEST(KeyManager, IntelMPIBug)
{
  MPI_Comm comm2;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm2);
  int keyval;
  MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, &comm_destructor_callback, &keyval, nullptr);
  MPI_Comm_set_attr(comm2, keyval, nullptr);

  is_complete = false;
  MPI_Ibarrier(comm2, &req);
  MPI_Comm_free(&comm2);
  EXPECT_TRUE(is_complete);
}
#endif


TEST_F(KeyManagerTester, KeyValuesUniqueness)
{
  MPI_Comm comm3, comm4;
  MPI_Comm_dup(comm1, &comm3);
  MPI_Comm_dup(comm1, &comm4);


  std::vector<stk::impl::MPIKeyManager::CommKey> keyVals;
  keyVals.push_back(keyManager.get_key(comm1));
  keyVals.push_back(keyManager.get_key(comm2));
  keyVals.push_back(keyManager.get_key(comm3));
  keyVals.push_back(keyManager.get_key(comm4));

  for (unsigned int i=0; i < keyVals.size(); ++i)
    for (unsigned int j=0; j < keyVals.size(); ++j)
      if (i == j)
        EXPECT_EQ(keyVals[i], keyVals[j]);
      else
        EXPECT_NE(keyVals[i], keyVals[j]);

  MPI_Comm_free(&comm3);
  MPI_Comm_free(&comm4);
}


TEST_F(KeyManagerTester, getTagConst)
{
  auto key1 = keyManager.get_key(comm1);
  auto key2 = keyManager.get_key(comm2);

  const stk::impl::MPIKeyManager& keyManagerConst = keyManager;

  EXPECT_EQ(key1, keyManagerConst.get_key(comm1));
  EXPECT_EQ(key2, keyManagerConst.get_key(comm2));

  MPI_Comm comm3;
  MPI_Comm_dup(comm1, &comm3);
  EXPECT_ANY_THROW(keyManagerConst.get_key(comm3));
  MPI_Comm_free(&comm3);
}


TEST_F(KeyManagerTester, KeyReuse)
{
  MPI_Comm comm3;
  MPI_Comm_dup(comm1, &comm3);

  auto key1 = keyManager.get_key(comm1);
  auto key2 = keyManager.get_key(comm2);
  keyManager.get_key(comm3);

  MPI_Comm_free(&comm3);
  MPI_Comm_dup(comm1, &comm3);

  auto key4 = keyManager.get_key(comm3);

  EXPECT_NE(key1, key4);
  EXPECT_NE(key2, key4);

  MPI_Comm_free(&comm3);
}


TEST_F(KeyManagerTester, SplitCommSubsetOfProcs)
{
  MPI_Comm comm3;
  int myrank;
  MPI_Comm_rank(comm1, &myrank);
  int color = myrank == 0 ? MPI_UNDEFINED : 1;
  MPI_Comm_split(comm1, color, 0, &comm3);

  auto key1 = keyManager.get_key(comm1);
  stk::impl::MPIKeyManager::CommKey key3;
  if (color != MPI_UNDEFINED)
    key3 = keyManager.get_key(comm3);
  auto key2 = keyManager.get_key(comm2);

  EXPECT_NE(key1, key2);
  if (color != MPI_UNDEFINED)
  {
    EXPECT_NE(key1, key3);
    EXPECT_NE(key2, key3);

    MPI_Comm_free(&comm3);
  }
}


TEST_F(KeyManagerTester, Callback)
{
  MPI_Comm comm3;
  MPI_Comm_dup(comm1, &comm3);

  CallbackChecker checker3;
  keyManager.get_key(comm3);

  EXPECT_EQ(checker3.callcount, 0);
  auto callerUID = keyManager.get_UID();
  keyManager.register_callback(comm3, callerUID,
    std::bind(&CallbackChecker::set_called, &checker3, std::placeholders::_1));
  MPI_Comm_free(&comm3);
  EXPECT_EQ(checker3.callcount, 1);  

}


TEST(KeyManager, CommNeverFreed)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  CallbackChecker checker;

  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm);
    auto callerUID = keyManager.get_UID();
    keyManager.register_callback(comm, callerUID,
      std::bind(&CallbackChecker::set_called, &checker, std::placeholders::_1));

   EXPECT_EQ(checker.callcount, 0);
  }

  EXPECT_EQ(checker.callcount, 1);
}  // comm is never freed, but MPIKeyManager shouldn't segfault


TEST(KeyManager, CommFreedAfterManager)
{
  MPI_Comm comm3;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm3);

  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm3);
  }  // keyManager destructor runs here, shouldn't segfault

  MPI_Comm_free(&comm3);
}


TEST(KeyManager, CommFreedAfterManagerWithCallback)
{
  MPI_Comm comm3;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm3);
  CallbackChecker checker3;

  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm3);
    auto callerUID = keyManager.get_UID();
    keyManager.register_callback(comm3, callerUID,
      std::bind(&CallbackChecker::set_called, &checker3, std::placeholders::_1));
  }

  EXPECT_EQ(checker3.callcount, 1);
  MPI_Comm_free(&comm3);
  EXPECT_EQ(checker3.callcount, 1);
}


TEST(KeyManager, CommFreedCallback)
{
  MPI_Comm comm3;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm3);
  CallbackChecker checker3;

  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm3);
    auto callerUID = keyManager.get_UID();
    keyManager.register_callback(comm3, callerUID,
      std::bind(&CallbackChecker::set_called, &checker3, std::placeholders::_1));

    EXPECT_EQ(checker3.callcount, 0);
    MPI_Comm_free(&comm3);
    EXPECT_EQ(checker3.callcount, 1);
  }

    
  EXPECT_EQ(checker3.callcount, 1);
}


TEST(KeyManager, CallbackImmediately)
{
  MPI_Comm comm3;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm3);
  CallbackChecker checker3;

  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm3);
    auto callerUID = keyManager.get_UID();
    keyManager.register_callback(comm3,
      callerUID, std::bind(&CallbackChecker::set_called, &checker3, std::placeholders::_1));
    keyManager.execute_callbacks_immediately(callerUID);

    EXPECT_EQ(checker3.callcount, 1);
  }

  EXPECT_EQ(checker3.callcount, 1);
  MPI_Comm_free(&comm3);
  EXPECT_EQ(checker3.callcount, 1);
}


TEST(KeyManager, UnregisterCallbacks)
{
  MPI_Comm comm3;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm3);
  CallbackChecker checker3;

  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm3);
    auto callerUID = keyManager.get_UID();
    keyManager.register_callback(comm3,
      callerUID, std::bind(&CallbackChecker::set_called, &checker3, std::placeholders::_1));
    keyManager.unregister_callbacks(callerUID);

    EXPECT_EQ(checker3.callcount, 0);
  }

  EXPECT_EQ(checker3.callcount, 0);
  MPI_Comm_free(&comm3);
  EXPECT_EQ(checker3.callcount, 0);
}


TEST(KeyManager, MultipleCallbacks)
{
  MPI_Comm comm3;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm3);
  CallbackChecker checker1, checker2, checker3;

  auto callback1 = std::bind(&CallbackChecker::set_called, &checker1, std::placeholders::_1);
  auto callback2 = std::bind(&CallbackChecker::set_called, &checker2, std::placeholders::_1);
  auto callback3 = std::bind(&CallbackChecker::set_called, &checker3, std::placeholders::_1);
  {
    stk::impl::MPIKeyManager keyManager;
    keyManager.get_key(comm3);
    auto callerUID1 = keyManager.get_UID();
    auto callerUID2 = keyManager.get_UID();
    auto callerUID3 = keyManager.get_UID();
    keyManager.register_callback(comm3, callerUID1, callback1);
    keyManager.register_callback(comm3, callerUID2, callback2);
    keyManager.register_callback(comm3, callerUID3, callback3);

    keyManager.execute_callbacks_immediately(callerUID1);
    EXPECT_EQ(checker1.callcount, 1);
    EXPECT_EQ(checker2.callcount, 0);
    EXPECT_EQ(checker3.callcount, 0);

    keyManager.unregister_callbacks(callerUID2);
    EXPECT_EQ(checker1.callcount, 1);
    EXPECT_EQ(checker2.callcount, 0);
    EXPECT_EQ(checker3.callcount, 0);
  }

  EXPECT_EQ(checker1.callcount, 1);
  EXPECT_EQ(checker2.callcount, 0);
  EXPECT_EQ(checker3.callcount, 1);

  MPI_Comm_free(&comm3);
  EXPECT_EQ(checker1.callcount, 1);
  EXPECT_EQ(checker2.callcount, 0);
  EXPECT_EQ(checker3.callcount, 1);
}


TEST(KeyManager, OrderedKeys)
{
  MPI_Comm comm1=MPI_COMM_WORLD, comm2, comm3, comm4;
  int myRank, commSize;
  MPI_Comm_rank(comm1, &myRank);
  MPI_Comm_size(comm1, &commSize);

  int color = myRank == commSize-1 ? MPI_UNDEFINED : 1;
  MPI_Comm_split(comm1, color, 0, &comm2);
  MPI_Comm_dup(comm1, &comm3);
  MPI_Comm_dup(comm1, &comm4);

  auto keyManager = std::make_shared<stk::impl::MPIKeyManager>();
  keyManager->get_key(comm1);
  if (color != MPI_UNDEFINED)
    keyManager->get_key(comm2);

  keyManager->get_key(comm3);
  if (color != MPI_UNDEFINED)
    MPI_Comm_free(&comm2);

  keyManager->get_key(comm4);

  std::vector<MPI_Comm> comms{comm1, comm3, comm4};
  std::sort(comms.begin(), comms.end(), stk::impl::CommCompare(keyManager));

  // if this doesn't hang, the comms are in the same order on all procs
  for (auto& comm : comms)
    MPI_Barrier(comm);

  MPI_Comm_free(&comm3);
  MPI_Comm_free(&comm4);
}

#endif