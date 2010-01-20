#include "Teuchos_TimeMonitor.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR  )
{
  TEUCHOS_FUNC_TIME_MONITOR("FUNC_TIME_MONITOR1");
  std::ostringstream oss;
  TimeMonitor::summarize(oss);
  out << oss.str() << "\n";
  const size_t substr_i = oss.str().find("FUNC_TIME_MONITOR1");
  TEST_INEQUALITY(substr_i, std::string::npos);
}


TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR_tested  )
{
  TEUCHOS_FUNC_TIME_MONITOR("FUNC_TIME_MONITOR2");
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("FUNC_TIME_MONITOR2_inner", inner);
  }
  std::ostringstream oss;
  TimeMonitor::summarize(oss);
  out << oss.str() << "\n";
  const size_t substr_i = oss.str().find("FUNC_TIME_MONITOR2");
  TEST_INEQUALITY(substr_i, std::string::npos);
  const size_t substr_inner_i = oss.str().find("FUNC_TIME_MONITOR2_inner");
  TEST_INEQUALITY(substr_inner_i, std::string::npos);
}


} // namespace Teuchos
