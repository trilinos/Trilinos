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

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for set_assert_handler, ThrowRequireMsg, set_erro...
#include "stk_util/ngp/NgpSpaces.hpp"
#include <iostream>                         // for basic_ostream::operator<<, operator<<, ostrin...
#include <stdexcept>                        // for logic_error, runtime_error, invalid_argument
#include <string>                           // for string


namespace {

bool test_assert_handler_called = false;
bool test_error_handler_called = false;
bool test_invarg_handler_called = false;

void
test_assert_handler(const char* expr,
                    const std::string& location,
                    std::ostringstream& message)
{
  test_assert_handler_called = true;
}

void
test_error_handler(const char* expr,
                   const std::string& location,
                   std::ostringstream& message)
{
  test_error_handler_called = true;
}

void
test_invarg_handler(const char* expr,
                    const std::string& location,
                    std::ostringstream& message)
{
  test_invarg_handler_called = true;
}

void force_throw_require_trigger(bool msg = true)
{
  if (msg) {
    STK_ThrowRequireMsg(false, "Will always fail, this is for testing");
  }
  else {
    STK_ThrowRequire(false);
  }
}

void force_throw_error_trigger(bool msg = true)
{
  if (msg) {
    STK_ThrowErrorMsgIf(true, "Will always fail, this is for testing");
  }
  else {
    STK_ThrowErrorIf(true);
  }
}

void force_throw_invarg_trigger(bool msg = true)
{
  if (msg) {
    STK_ThrowInvalidArgMsgIf(true, "Will always fail, this is for testing");
  }
  else {
    STK_ThrowInvalidArgIf(true);
  }
}

void check_interaction_with_if(bool msg = true)
{
  if (msg) {
    if (true)
      STK_ThrowRequireMsg(false, "Will always fail, this is for testing");
  }
  else {
    if (true)
      STK_ThrowRequire(false);
  }
}

void force_throw_assert()
{
  STK_ThrowAssert(false);
}

void test_no_expr_error()
{
  STK_ThrowErrorMsg("message");
}

} // namespace <empty>

TEST(UnitTestingOfThrowMacros, testUnit)
{
  // Setting assert handler to NULL should cause exception
  ASSERT_THROW(stk::set_assert_handler(0), std::runtime_error);

  // Check that Throw*Msg works
  ASSERT_THROW(force_throw_require_trigger(), std::logic_error);
  ASSERT_THROW(force_throw_error_trigger(), std::runtime_error);
  ASSERT_THROW(force_throw_invarg_trigger(), std::invalid_argument);

  // Check that Throw* works
  ASSERT_THROW(force_throw_require_trigger(false), std::logic_error);
  ASSERT_THROW(force_throw_error_trigger(false), std::runtime_error);
  ASSERT_THROW(force_throw_invarg_trigger(false), std::invalid_argument);

  // Check that macro interacts appropriately with if statements
  ASSERT_THROW(check_interaction_with_if(), std::logic_error);
  ASSERT_THROW(check_interaction_with_if(false), std::logic_error);

  // Check that usage of ThrowRequireMsg/ThrowAssertMsg does not change program
  // semantics. Code blocks that are not contained within braces seem to be
  // the most likely to be problematic.

  bool expected_execution_path = false;
  if (false)
    STK_ThrowRequireMsg(false, "test");
  else
    expected_execution_path = true;
  ASSERT_TRUE(expected_execution_path);

  expected_execution_path = false;
  if (false)
    STK_ThrowAssertMsg(false, "test");
  else
    expected_execution_path = true;
  ASSERT_TRUE(expected_execution_path);

  expected_execution_path = false;
  if (false)
    STK_ThrowErrorMsg("test");
  else
    expected_execution_path = true;
  ASSERT_TRUE(expected_execution_path);

  // These next four statements are to check compilation success

  if (false)
    STK_ThrowRequireMsg(false, "test");

  if (false)
    STK_ThrowAssertMsg(false, "test");

  if (false)
    STK_ThrowRequire(false);

  if (false)
    STK_ThrowAssert(false);

  // Check that do-while still works, again, we are mostly checking compilation
  // success here.

  do STK_ThrowRequireMsg(true, "test");
  while (false);

  do STK_ThrowAssertMsg(true, "test");
  while (false);

  // Check that message with put-tos compiles

  int temp = 0;
  STK_ThrowRequireMsg(true, "test: " << temp << " blah");
  STK_ThrowAssertMsg(true, "test: " << temp << " blah");

  // Check that assert behaves as expected (throws in debug, not in opt)
#ifdef NDEBUG
  force_throw_assert();
#else
  ASSERT_THROW(force_throw_assert(), std::logic_error);
#endif

  // Check that ThrowErrorMsg works

  ASSERT_THROW(test_no_expr_error(), std::runtime_error);
  
  // Check that setting handler for asserts works.

  stk::ErrorHandler orig = stk::set_assert_handler(test_assert_handler);

  STK_ThrowRequireMsg(false, "test");

  ASSERT_TRUE(test_assert_handler_called);

  stk::set_assert_handler(orig);

  ASSERT_THROW(force_throw_require_trigger(), std::logic_error);

  // Check that setting handler for errors works.

  orig = stk::set_error_handler(test_error_handler);

  STK_ThrowErrorMsgIf(true, "test");

  ASSERT_TRUE(test_error_handler_called);

  stk::set_error_handler(orig);

  ASSERT_THROW(force_throw_error_trigger(), std::runtime_error);

  // Check that setting handler for invalid args works.

  orig = stk::set_invalid_arg_handler(test_invarg_handler);

  STK_ThrowInvalidArgMsgIf(true, "test");

  ASSERT_TRUE(test_invarg_handler_called);

  stk::set_invalid_arg_handler(orig);

  ASSERT_THROW(force_throw_invarg_trigger(), std::invalid_argument);
}

void testNGPThrowRequireMsg()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    bool test = false;
    STK_NGP_ThrowRequireMsg(test == true, "Error testing whatever");
  });
}

#if (defined(__GNUC__) && (__GNUC__ > 4))
#define NEW_ENOUGH_GCC
#endif

TEST(UnitTestingOfThrowMacros, NGP_ThrowRequireMsg)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  // Also, OpenMP seems to produce an abort (in adddition to a throw?).
  //
  // testNGPThrowRequireMsg();
  std::cout<<"NGP_ThrowRequireMsg: #if cuda, openmp or hip"<<std::endl;
#else
#ifdef NEW_ENOUGH_GCC
  std::cout<<"NGP_ThrowRequireMsg: #ifdef new-enough-gcc"<<std::endl;
//For now, only test this on gcc compilers more recent than major version 4.
//A Trilinos pre-push test platform, 4.8.4 seems to produce an abort instead
//of a throw for this test.
  try {
    testNGPThrowRequireMsg();
  }
  catch (std::exception & ex) {
    std::cerr<<"ex.what(): "<<ex.what()<<std::endl;
    std::string expectedMsg1 = "Requirement( test == true ) FAILED\n"
                               "Error occurred at: stk_unit_tests/stk_util/util/UnitTestThrowMacros.cpp:";
    std::string expectedMsg2 = "Error: Error testing whatever\n";
    std::string message = ex.what();
  std::cout<<"NGP_ThrowRequireMsg: caught, comparing strings"<<std::endl;
    EXPECT_NE(message.find(expectedMsg1), std::string::npos);
    EXPECT_NE(message.find(expectedMsg2), std::string::npos);
  }
  catch (...) {
    std::cerr<<"Unexpected exception-type thrown."<<std::endl;
  }
#endif
#endif
}

void testNGPThrowRequire()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    bool test = false;
    STK_NGP_ThrowRequire(test == true);
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowRequire)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  //
  // testNGPThrowRequire();
  std::cout<<"NGP_ThrowRequire: #if cuda, openmp or hip"<<std::endl;
#else
#ifdef NEW_ENOUGH_GCC
  std::cout<<"NGP_ThrowRequire: #ifdef new-enough-gcc"<<std::endl;
  try {
    testNGPThrowRequire();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( test == true ) FAILED\n"
                               "Error occurred at: stk_unit_tests/stk_util/util/UnitTestThrowMacros.cpp:";
    std::string message = ex.what();
  std::cout<<"NGP_ThrowRequire: caught, comparing strings"<<std::endl;
    EXPECT_NE(message.find(expectedMsg), std::string::npos);
  }
#endif
#endif
}

// Debug testing
#ifndef NDEBUG
void testNGPThrowAssertMsg()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    bool test = false;
    STK_NGP_ThrowAssertMsg(test == true, "Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowAssertMsg_debug)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  //
  // testNGPThrowAssertMsg();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowAssertMsg();
  }
  catch (std::exception & ex) {
    const char * expectedMsg1 = "Requirement( test == true ) FAILED\n"
                               "Error occurred at: stk_unit_tests/stk_util/util/UnitTestThrowMacros.cpp:";
    const char * expectedMsg2 = "Error: Error testing whatever\n";
    std::string message = ex.what();
    EXPECT_NE(message.find(expectedMsg1), std::string::npos);
    EXPECT_NE(message.find(expectedMsg2), std::string::npos);
  }
#endif
#endif
}
#endif

// Release testing
#ifdef NDEBUG
void testNGPThrowAssertMsg()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    STK_NGP_ThrowAssertMsg(false, "Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowAssertMsg_release)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  //
  // testNGPThrowAssertMsg();
#else
#ifdef NEW_ENOUGH_GCC
  EXPECT_NO_THROW(testNGPThrowAssertMsg());
#endif
#endif
}
#endif

void testNGPThrowErrorMsgIf()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    bool test = true;
    STK_NGP_ThrowErrorMsgIf(test == true, "Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowErrorMsgIf)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  //
  // testNGPThrowErrorMsgIf();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowErrorMsgIf();
  }
  catch (std::exception & ex) {
    const char * expectedMsg1 = "Requirement( !(test == true) ) FAILED\n"
                               "Error occurred at: stk_unit_tests/stk_util/util/UnitTestThrowMacros.cpp:";
    const char * expectedMsg2 = "Error: Error testing whatever\n";
    std::string message = ex.what();
    EXPECT_NE(message.find(expectedMsg1), std::string::npos);
    EXPECT_NE(message.find(expectedMsg2), std::string::npos);
  }
#endif
#endif
}

void testNGPThrowErrorIf()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    bool test = true;
    STK_NGP_ThrowErrorIf(test == true);
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowErrorIf)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  //
  // testNGPThrowErrorIf();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowErrorIf();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( !(test == true) ) FAILED\n"
                               "Error occurred at: stk_unit_tests/stk_util/util/UnitTestThrowMacros.cpp:";
    std::string message = ex.what();
    EXPECT_NE(message.find(expectedMsg), std::string::npos);
  }
#endif
#endif
}

void testNGPThrowErrorMsg()
{
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int & i){
    STK_NGP_ThrowErrorMsg("Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowErrorMsg)
{
#if defined(STK_ENABLE_GPU) || defined(_OPENMP)
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize().
  //
  // testNGPThrowErrorMsg();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowErrorMsg();
  }
  catch (std::exception & ex) {
    const char * expectedMsg1 = "Error occurred at: stk_unit_tests/stk_util/util/UnitTestThrowMacros.cpp:";
    const char * expectedMsg2 = "Error: Error testing whatever\n";
    std::string message = ex.what();
    EXPECT_NE(message.find(expectedMsg1), std::string::npos);
    EXPECT_NE(message.find(expectedMsg2), std::string::npos);
  }
#endif
#endif
}

#if !defined(NDEBUG)
// clang-format off
#define TEST_THROW_MACROS_IN_BLOCKS_ARG(testThrow, cond) \
{                                                        \
  for (int i = 0; i < 1; ++i) {                          \
    testThrow(++i == -1 || cond);                        \
  }                                                      \
  for (int i = 0; i < 1; ++i)                            \
    testThrow(++i == -1 || cond);                        \
  if (int i = 1) {                                       \
    testThrow(++i == -1 || cond);                        \
  }                                                      \
  if (int i = 1)                                         \
    testThrow(++i == -1 || cond);                        \
}

#define TEST_THROW_MSG_MACROS_IN_BLOCKS_ARG(testThrowMsg, cond) \
{                                                               \
  for (int i = 0; i < 1; ++i) {                                 \
    testThrowMsg(++i == -1 || cond, "");                        \
  }                                                             \
  for (int i = 0; i < 1; ++i)                                   \
    testThrowMsg(++i == -1 || cond, "");                        \
  if (int i = 1) {                                              \
    testThrowMsg(++i == -1 || cond, "");                        \
  }                                                             \
  if (int i = 1)                                                \
    testThrowMsg(++i == -1 || cond, "");                        \
}
// clang-format on

#define TEST_THROW_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MACROS_IN_BLOCKS_ARG(testThrow, false)
#define TEST_THROW_IF_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MACROS_IN_BLOCKS_ARG(testThrow, true)
#define TEST_THROW_MSG_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MSG_MACROS_IN_BLOCKS_ARG(testThrow, false)
#define TEST_THROW_MSG_IF_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MSG_MACROS_IN_BLOCKS_ARG(testThrow, true)

#define TEST_NGP_THROW_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MACROS_IN_BLOCKS_ARG(testThrow, true)
#define TEST_NGP_THROW_IF_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MACROS_IN_BLOCKS_ARG(testThrow, false)
#define TEST_NGP_THROW_MSG_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MSG_MACROS_IN_BLOCKS_ARG(testThrow, true)
#define TEST_NGP_THROW_MSG_IF_MACROS_IN_BLOCKS(testThrow)  TEST_THROW_MSG_MACROS_IN_BLOCKS_ARG(testThrow, false)

TEST(UnitTestingOfThrowMacros, STK_ThrowsInSingleStatementBlock)
{
  EXPECT_ANY_THROW(TEST_THROW_MACROS_IN_BLOCKS(STK_ThrowRequireWithSierraHelpMsg));
  EXPECT_ANY_THROW(TEST_THROW_MACROS_IN_BLOCKS(STK_ThrowRequire));
  EXPECT_ANY_THROW(TEST_THROW_IF_MACROS_IN_BLOCKS(STK_ThrowErrorIf));
  EXPECT_ANY_THROW(TEST_THROW_IF_MACROS_IN_BLOCKS(STK_ThrowInvalidArgIf));

  EXPECT_ANY_THROW(TEST_THROW_MSG_MACROS_IN_BLOCKS(STK_ThrowRequireMsg));
  EXPECT_ANY_THROW(TEST_THROW_MSG_IF_MACROS_IN_BLOCKS(STK_ThrowErrorMsgIf));
  EXPECT_ANY_THROW(TEST_THROW_MSG_IF_MACROS_IN_BLOCKS(STK_ThrowInvalidArgMsgIf));

#ifdef NDEBUG
  EXPECT_NO_THROW(TEST_THROW_MACROS_IN_BLOCKS(STK_ThrowAssert));
  EXPECT_NO_THROW(TEST_THROW_MSG_MACROS_IN_BLOCKS(STK_ThrowAssertMsg));
#else
  EXPECT_ANY_THROW(TEST_THROW_MACROS_IN_BLOCKS(STK_ThrowAssert));
  EXPECT_ANY_THROW(TEST_THROW_MSG_MACROS_IN_BLOCKS(STK_ThrowAssertMsg));
#endif
}

template <typename Lambda>
void run_ngp_macros_test(Lambda lambda) {
  EXPECT_NO_THROW( Kokkos::parallel_for(1, lambda) );
}

template <typename... Lambdas>
void run_ngp_macros_in_single_statement_block(Lambdas... lambdas) {
  (run_ngp_macros_test(lambdas), ...);
}

void test_stk_ngp_macros_in_single_statement_block() {
  auto ngpThrowRequire = KOKKOS_LAMBDA(const int) { TEST_NGP_THROW_MACROS_IN_BLOCKS(STK_NGP_ThrowRequire); };
  auto ngpThrowAssert = KOKKOS_LAMBDA(const int) { TEST_NGP_THROW_MACROS_IN_BLOCKS(STK_NGP_ThrowAssert); };
  auto ngpThrowErrorIf = KOKKOS_LAMBDA(const int) { TEST_NGP_THROW_IF_MACROS_IN_BLOCKS(STK_NGP_ThrowErrorIf); };

  auto ngpThrowRequireMsg = KOKKOS_LAMBDA(const int) { TEST_NGP_THROW_MSG_MACROS_IN_BLOCKS(STK_NGP_ThrowRequireMsg); };
  auto ngpThrowAssertMsg = KOKKOS_LAMBDA(const int) { TEST_NGP_THROW_MSG_MACROS_IN_BLOCKS(STK_NGP_ThrowAssertMsg); };
  auto ngpThrowErrorMsgIf = KOKKOS_LAMBDA(const int) { TEST_NGP_THROW_MSG_IF_MACROS_IN_BLOCKS(STK_NGP_ThrowErrorMsgIf); };

  run_ngp_macros_in_single_statement_block(ngpThrowRequire, ngpThrowAssert, ngpThrowErrorIf, ngpThrowRequireMsg,
                                           ngpThrowRequire, ngpThrowAssertMsg,  ngpThrowErrorMsgIf);
}

TEST(UnitTestingOfThrowMacros, STK_NGPThrowsInSingleStatementBlock)
{
  test_stk_ngp_macros_in_single_statement_block();
}

#endif
