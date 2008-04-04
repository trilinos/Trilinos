#include "../CppUnitLite/TestHarness.h"

#include "Rythmos_DataStore_UnitTest.hpp"
#include "Rythmos_InterpolationBuffer_UnitTest.hpp"
#include "Rythmos_TimeRange_UnitTest.hpp"

int main()
{
  TestResult tr;
  TestRegistry::runAllTests(tr);
  return 0;
}
