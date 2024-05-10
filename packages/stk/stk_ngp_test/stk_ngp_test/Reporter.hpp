#ifndef STK_NGP_TEST_REPORTER_HPP
#define STK_NGP_TEST_REPORTER_HPP

#include "NgpTestDeviceMacros.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_ErrorReporter.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>

namespace ngp_testing {

class TruncatedString {
 public:
  static constexpr int maxNumChar = 63;
  NGP_TEST_FUNCTION TruncatedString() {}
  NGP_TEST_FUNCTION TruncatedString(const char* str);

  operator const char*() const { return string; }

 private:
  char string[maxNumChar + 1] = {0};

  NGP_TEST_FUNCTION int get_string_length(const char* str) const;
  NGP_TEST_FUNCTION bool should_truncate(const int strLen) const;
  NGP_TEST_FUNCTION int get_truncation_offset(const int strLen) const;
  NGP_TEST_FUNCTION void prepend_truncation_indicator(char* dest);
};

struct Report {
  NGP_TEST_FUNCTION Report()
   : condition(), location() {}
  NGP_TEST_FUNCTION Report(const char* cond, const char* loc)
   : condition(cond), location(loc) {}


  TruncatedString condition;
  TruncatedString location;
};

inline std::ostream& operator<<(std::ostream& out, const Report& r) {
  out << "At location: " << r.location << ", failed condition: " << r.condition;
  return out;
}

template<class Device>
class Reporter {
 public:
  Reporter(int numReports) : reporter(numReports) {}
  void clear() { reporter.clear(); }
  void resize(int n) { reporter.resize(n); }
  int get_capacity() const { return reporter.getCapacity(); }

  int report_failures(std::ostream& out = std::cout) {
    int numAttemptedReports = reporter.getNumReportAttempts();
    int numReports = reporter.getNumReports();

    if (numAttemptedReports != 0) {
      out << "In test " << get_gtest_name() << ", " << numAttemptedReports << " failures occurred." << std::endl
          << "    Limited summary of failures: " << std::endl;

      std::vector<ReportT> failures;
      std::vector<int> failureWorkIndices;
      reporter.getReports(failureWorkIndices, failures);

      for (int i = 0; i < numReports; ++i) {
        auto&& failure = failures[i];
        out << "      " << failure << std::endl;
      }
    }

    return numAttemptedReports;
  }

  NGP_TEST_INLINE
  void add_failure(const char* condition,
                   const char* location) const {
    static const int dummyWorkIndex = -1;
    reporter.add_report(dummyWorkIndex, ReportT(condition, location));
  }

 private:
  using ReportT = Report;
  Kokkos::Experimental::ErrorReporter<ReportT, Device> reporter;

  inline std::string get_gtest_name() const {
    std::stringstream ss;
    const ::testing::TestInfo* const info = ::testing::UnitTest::GetInstance()->current_test_info();
    ss << info->test_case_name() << "." << info->name();
    return ss.str();
  }
};

}

namespace ngp_testing {

NGP_TEST_INLINE void bounded_strcpy(const char* src, const int maxChar, char* dest) {
  int idx = 0;
  if(src) {
    while(src[idx] != '\0' && idx < maxChar) {
      dest[idx] = src[idx];
      ++idx;
    }
  }
  dest[idx] = '\0';
}

NGP_TEST_INLINE
TruncatedString::TruncatedString(const char* src) {
  const int srcLen = get_string_length(src);
  if(should_truncate(srcLen)) {
    src += get_truncation_offset(srcLen);
  }
  bounded_strcpy(src, maxNumChar, string);
  if(should_truncate(srcLen)) {
    prepend_truncation_indicator(string);
  }
}

NGP_TEST_INLINE
int TruncatedString::get_string_length(const char* str) const {
  int len = 0;
  if(str) {
    while (str[len] != '\0') ++len;
  }
  return len;
}

NGP_TEST_INLINE
bool TruncatedString::should_truncate(const int strLen) const {
  return strLen > maxNumChar;
}

NGP_TEST_INLINE
int TruncatedString::get_truncation_offset(const int strLen) const {
  return strLen - maxNumChar + 1;
}

NGP_TEST_INLINE
void TruncatedString::prepend_truncation_indicator(char* dest) {
  dest[0] = '.';
  dest[1] = '.';
  dest[2] = '.';
}

}

#endif
