#ifndef STK_NGP_TEST_GLOBAL_REPORTER_HPP
#define STK_NGP_TEST_GLOBAL_REPORTER_HPP

#include "NgpTestDeviceMacros.hpp"

namespace ngp_testing {

class ReporterBase;

void initialize_reporters();
void finalize_reporters();

NGP_TEST_FUNCTION ReporterBase* get_reporter();
ReporterBase* get_device_reporter();

}

#include "GlobalReporter.cpp"

#endif
