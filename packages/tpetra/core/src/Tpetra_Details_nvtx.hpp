#pragma once

#if defined(TPETRA_HAVE_NVTX)
#include <nvToolsExt.h>
#define NVTX_RANGE_PUSH(s) nvtxRangePush(s)
#define NVTX_RANGE_POP() nvtxRangePop()
#define NVTX_MARK(s) nvtxMark(s)
#else
#define NVTX_RANGE_PUSH(s) ((void)(s))
#define NVTX_RANGE_POP()
#define NVTX_MARK(s) ((void)(s))
#endif

namespace Tpetra {
namespace Details {

class Range {
public:
  Range() : owns_(false) {}
  Range(const Range &other) = delete;
  Range(Range &&other) = delete;
  Range(const char *s) : owns_(true) {
    NVTX_RANGE_PUSH(s);
  }
  Range(const std::string &s) : Range(s.c_str()) {}

  ~Range() {
    if (owns_) {
      NVTX_RANGE_POP();
    }
  }
private:
  bool owns_;
};

inline void mark(const char *s) {
  NVTX_MARK(s);
}

}
}
