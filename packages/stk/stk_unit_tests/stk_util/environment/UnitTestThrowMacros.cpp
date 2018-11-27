// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <iostream>                     // for ostringstream, etc
#include <stdexcept>                    // for logic_error, runtime_error, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg, etc
#include <gtest/gtest.h>
#include <string>                       // for string


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
    ThrowRequireMsg(false, "Will always fail, this is for testing");
  }
  else {
    ThrowRequire(false);
  }
}

void force_throw_error_trigger(bool msg = true)
{
  if (msg) {
    ThrowErrorMsgIf(true, "Will always fail, this is for testing");
  }
  else {
    ThrowErrorIf(true);
  }
}

void force_throw_invarg_trigger(bool msg = true)
{
  if (msg) {
    ThrowInvalidArgMsgIf(true, "Will always fail, this is for testing");
  }
  else {
    ThrowInvalidArgIf(true);
  }
}

void check_interaction_with_if(bool msg = true)
{
  if (msg) {
    if (true)
      ThrowRequireMsg(false, "Will always fail, this is for testing");
  }
  else {
    if (true)
      ThrowRequire(false);
  }
}

void force_throw_assert()
{
  ThrowAssert(false);
}

void test_no_expr_error()
{
  ThrowErrorMsg("message");
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
    ThrowRequireMsg(false, "test");
  else
    expected_execution_path = true;
  ASSERT_TRUE(expected_execution_path);

  expected_execution_path = false;
  if (false)
    ThrowAssertMsg(false, "test");
  else
    expected_execution_path = true;
  ASSERT_TRUE(expected_execution_path);

  expected_execution_path = false;
  if (false)
    ThrowErrorMsg("test");
  else
    expected_execution_path = true;
  ASSERT_TRUE(expected_execution_path);

  // These next four statements are to check compilation success

  if (false)
    ThrowRequireMsg(false, "test");

  if (false)
    ThrowAssertMsg(false, "test");

  if (false)
    ThrowRequire(false);

  if (false)
    ThrowAssert(false);

  // Check that do-while still works, again, we are mostly checking compilation
  // success here.

  do ThrowRequireMsg(true, "test");
  while (false);

  do ThrowAssertMsg(true, "test");
  while (false);

  // Check that message with put-tos compiles

  int temp = 0;
  ThrowRequireMsg(true, "test: " << temp << " blah");
  ThrowAssertMsg(true, "test: " << temp << " blah");

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

  ThrowRequireMsg(false, "test");

  ASSERT_TRUE(test_assert_handler_called);

  stk::set_assert_handler(orig);

  ASSERT_THROW(force_throw_require_trigger(), std::logic_error);

  // Check that setting handler for errors works.

  orig = stk::set_error_handler(test_error_handler);

  ThrowErrorMsgIf(true, "test");

  ASSERT_TRUE(test_error_handler_called);

  stk::set_error_handler(orig);

  ASSERT_THROW(force_throw_error_trigger(), std::runtime_error);

  // Check that setting handler for invalid args works.

  orig = stk::set_invalid_arg_handler(test_invarg_handler);

  ThrowInvalidArgMsgIf(true, "test");

  ASSERT_TRUE(test_invarg_handler_called);

  stk::set_invalid_arg_handler(orig);

  ASSERT_THROW(force_throw_invarg_trigger(), std::invalid_argument);
}

void testNGPThrowRequireMsg()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    bool test = false;
    NGP_ThrowRequireMsg(test == true, "Error testing whatever");
  });
}

#if (defined(__GNUC__) && (__GNUC__ > 4))
#define NEW_ENOUGH_GCC
#endif

TEST(UnitTestingOfThrowMacros, NGP_ThrowRequireMsg)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
  //
  // testNGPThrowRequireMsg();
#else
#ifdef NEW_ENOUGH_GCC
//For now, only test this on gcc compilers more recent than major version 4.
//A Trilinos pre-push test platform, 4.8.4 seems to produce an abort instead
//of a throw for this test.
  try {
    testNGPThrowRequireMsg();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( test == true ) FAILED\n"
                               "Error occured at: stk_unit_tests/stk_util/environment/UnitTestThrowMacros.cpp:250\n"
                               "Error: Error testing whatever\n";
    EXPECT_STREQ(ex.what(), expectedMsg);
  }
#endif
#endif
}

void testNGPThrowRequire()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    bool test = false;
    NGP_ThrowRequire(test == true);
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowRequire)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
  //
  // testNGPThrowRequire();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowRequire();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( test == true ) FAILED\n"
                               "Error occured at: stk_unit_tests/stk_util/environment/UnitTestThrowMacros.cpp:287\n";
    EXPECT_STREQ(ex.what(), expectedMsg);
  }
#endif
#endif
}

// Debug testing
#ifndef NDEBUG
void testNGPThrowAssertMsg()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    bool test = false;
    NGP_ThrowAssertMsg(test == true, "Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowAssertMsg_debug)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
  //
  // testNGPThrowAssertMsg();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowAssertMsg();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( test == true ) FAILED\n"
                               "Error occured at: stk_unit_tests/stk_util/environment/UnitTestThrowMacros.cpp:318\n"
                               "Error: Error testing whatever\n";
    EXPECT_STREQ(ex.what(), expectedMsg);
  }
#endif
#endif
}
#endif

// Release testing
#ifdef NDEBUG
void testNGPThrowAssertMsg()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    NGP_ThrowAssertMsg(false, "Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowAssertMsg_release)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
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
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    bool test = true;
    NGP_ThrowErrorMsgIf(test == true, "Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowErrorMsgIf)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
  //
  // testNGPThrowErrorMsgIf();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowErrorMsgIf();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( !(test == true) ) FAILED\n"
                               "Error occured at: stk_unit_tests/stk_util/environment/UnitTestThrowMacros.cpp:373\n"
                               "Error: Error testing whatever\n";
    EXPECT_STREQ(ex.what(), expectedMsg);
  }
#endif
#endif
}

void testNGPThrowErrorIf()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    bool test = true;
    NGP_ThrowErrorIf(test == true);
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowErrorIf)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
  //
  // testNGPThrowErrorIf();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowErrorIf();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Requirement( !(test == true) ) FAILED\n"
                               "Error occured at: stk_unit_tests/stk_util/environment/UnitTestThrowMacros.cpp:403\n";
    EXPECT_STREQ(ex.what(), expectedMsg);
  }
#endif
#endif
}

void testNGPThrowErrorMsg()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & i){
    NGP_ThrowErrorMsg("Error testing whatever");
  });
}

TEST(UnitTestingOfThrowMacros, NGP_ThrowErrorMsg)
{
#ifdef KOKKOS_HAVE_CUDA
  // Unable to test a device-side abort, as this eventually results in a throw
  // inside Kokkos::finalize_all().
  //
  // testNGPThrowErrorMsg();
#else
#ifdef NEW_ENOUGH_GCC
  try {
    testNGPThrowErrorMsg();
  }
  catch (std::exception & ex) {
    const char * expectedMsg = "Error occured at: stk_unit_tests/stk_util/environment/UnitTestThrowMacros.cpp:431\n"
                               "Error: Error testing whatever\n";
    EXPECT_STREQ(ex.what(), expectedMsg);
  }
#endif
#endif
}

