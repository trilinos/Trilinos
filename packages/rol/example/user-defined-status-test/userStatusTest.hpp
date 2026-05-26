
#include "ROL_StatusTest.hpp"
#include <fstream>

template<class Real>
class FileStatusTest : public ROL::StatusTest<Real> {
private:
  std::string filename_;

public:
  FileStatusTest(std::string filename = "terminate.txt")
    : filename_(filename) {
    std::ifstream file(filename_);
    ROL_TEST_FOR_EXCEPTION(file.is_open(),std::invalid_argument,
      ">>> UserStatusTest: Terminate file already exists!  Please remove and try again.");
  }

  bool check(ROL::AlgorithmState<Real> &state) {
    std::ifstream file(filename_);
    bool flag = file.is_open();
    if (flag) state.statusFlag = ROL::EXITSTATUS_USERDEFINED;
    return !flag;
  }
};

#include <chrono>

template<class Real>
class TimeOutStatusTest : public ROL::StatusTest<Real> {
private:
  std::chrono::steady_clock::time_point startTime_;
  Real timeLimit_; // Time limit in seconds

public:
  TimeOutStatusTest(std::chrono::steady_clock::time_point startTime, Real timeLimit = 10.0)
    : startTime_(startTime), timeLimit_(timeLimit) {}

  bool check(ROL::AlgorithmState<Real> &state) {
    auto currentTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::duration<Real>>(currentTime-startTime_);
    bool flag = elapsedTime.count() >= timeLimit_;
    if (flag) state.statusFlag = ROL::EXITSTATUS_USERDEFINED;
    return !flag;
  }
};
