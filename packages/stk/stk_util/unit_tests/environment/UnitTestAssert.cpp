/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <iostream>
#include <stdexcept>

#include <stk_util/environment/ReportHandler.hpp>

namespace {

bool test_assert_handler_called = false;

void
test_assert_handler(const char* expr,
                    const std::string& location,
                    std::ostringstream& message)
{
  test_assert_handler_called = true;
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

} // namespace <empty>



STKUNIT_UNIT_TEST(UnitTestingOfThrowRequireMsg, testUnit)
{
  // Setting assert handler to NULL should cause exception
  STKUNIT_ASSERT_THROW(stk::set_assert_handler(0), std::runtime_error);

  // Check that ThrowRequireMsg works
  STKUNIT_ASSERT_THROW(force_throw_require_trigger(), std::logic_error);

  STKUNIT_ASSERT_THROW(check_interaction_with_if(), std::logic_error);

  // Check that ThrowRequire works
  STKUNIT_ASSERT_THROW(force_throw_require_trigger(false), std::logic_error);

  STKUNIT_ASSERT_THROW(check_interaction_with_if(false), std::logic_error);

  // Check that usage of ThrowRequireMsg/ThrowAssertMsg does not change program
  // semantics. Code blocks that are not contained within braces seem to be
  // the most likely to be problematic.

  bool expected_execution_path = false;
  if (false)
    ThrowRequireMsg(false, "test");
  else
    expected_execution_path = true;
  STKUNIT_ASSERT(expected_execution_path);

  expected_execution_path = false;
  if (false)
    ThrowAssertMsg(false, "test");
  else
    expected_execution_path = true;
  STKUNIT_ASSERT(expected_execution_path);

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
  STKUNIT_ASSERT_THROW(force_throw_assert(), std::logic_error);
#endif

  // Check that setting handler works.

  stk::ErrorHandler orig = stk::set_assert_handler(test_assert_handler);

  ThrowRequireMsg(false, "test");

  STKUNIT_ASSERT(test_assert_handler_called);

  stk::set_assert_handler(orig);

  STKUNIT_ASSERT_THROW(force_throw_require_trigger(), std::logic_error);
}
