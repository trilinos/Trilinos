#pragma once

#if defined(TPETRA_HAVE_NVTX) && TPETRA_HAVE_NVTX == 1
#include <nvToolsExt.h>
#define TPETRA_NVTX_RANGE_PUSH(x) nvtxRangePush(x)
#define TPETRA_NVTX_RANGE_POP() nvtxRangePop()
#define TPETRA_NVTX_MARK(x) nvtxMark(x)
#else
#define TPETRA_NVTX_RANGE_PUSH(x) ((void)(x))
#define TPETRA_NVTX_RANGE_POP()
#define TPETRA_NVTX_MARK(x) ((void)(x))
#endif

struct NvtxRange {
  NvtxRange(const char *msg) {
    TPETRA_NVTX_RANGE_PUSH(msg);
  }
  ~NvtxRange() {
    TPETRA_NVTX_RANGE_POP();
  }
  NvtxRange(const NvtxRange &rhs) = delete;
};
