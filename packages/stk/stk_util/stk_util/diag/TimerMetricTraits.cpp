#include <sstream>

#include <stk_util/diag/TimerMetricTraits.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/FormatTime.hpp>

namespace stk {
namespace diag {

namespace {

int s_timeFormat = TIMEFORMAT_HMS | TIMEFORMAT_MILLIS;

} // namespace <empty>


int
getTimerTimeFormat() 
{
  return s_timeFormat;
}

void
setTimerTimeFormat(
  int           time_format)
{
  s_timeFormat = time_format;
}


MetricTraits<LapCount>::Type
MetricTraits<LapCount>::value_now()
{
  return 1;
}

MetricTraits<CPUTime>::Type
MetricTraits<CPUTime>::value_now()
{
  return stk::cpu_time();
}

MetricTraits<WallTime>::Type
MetricTraits<WallTime>::value_now()
{
  return stk::wall_time();
}

MetricTraits<MPICount>::Type
MetricTraits<MPICount>::value_now()
{
  return 0;
}

MetricTraits<MPIByteCount>::Type
MetricTraits<MPIByteCount>::value_now()
{
  return 0;
}

std::string
MetricTraits<LapCount>::table_header() {
  return "Count";
}

std::string
MetricTraits<CPUTime>::table_header() {
  return "CPU Time";
}

std::string
MetricTraits<WallTime>::table_header() {
  return "Wall Time";
}

std::string
MetricTraits<MPICount>::table_header() {
  return "MPI Count";
}

std::string
MetricTraits<MPIByteCount>::table_header() {
  return "MPI Byte Count";
}


std::string
MetricTraits<CPUTime>::format(
  MetricTraits<CPUTime>::Type   time)
{
  return formatTime(time, getTimerTimeFormat());
}


std::string
MetricTraits<WallTime>::format(
  MetricTraits<WallTime>::Type   time)
{
  return formatTime(time, getTimerTimeFormat());
}


std::string
MetricTraits<MPICount>::format(
  MetricTraits<MPICount>::Type   count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}


std::string
MetricTraits<MPIByteCount>::format(
  MetricTraits<MPIByteCount>::Type      count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}

} // namespace diag
} // namespace stk
