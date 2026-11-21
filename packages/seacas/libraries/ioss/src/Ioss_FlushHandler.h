// Copyright(C) 1999-2020, 2022, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "Ioss_CodeTypes.h"
#include "Ioss_ParallelUtils.h"
#include "ioss_export.h"
#include <ctime>
#include <functional>

namespace Ioss {

  // Flush the files buffer to disk...
  // If:
  //  flushInterval == -1 (default) -- flush if there is more
  // than 10 seconds since the last flush to avoid
  // the flush eating up cpu time for small fast jobs...
  //
  //  flushInterval == 0 -- do not flush until file is closed.
  //
  //  flushInterval == 1 -- flush every step
  //
  //  flushInterval > 1 -- flush if step % flushInterval == 0
  //
  //  if time between begin_state and end_state is > 10 seconds,

  class IOSS_EXPORT FlushHandler
  {
  private:
    int                                               flushInterval;
    bool                                              isParallel;
    bool                                              isFirstOutput;
    bool                                              flushOnFirstOutput;
    unsigned int                                      timeLastFlushInterval;
    unsigned int                                      timeStepBeginFlushInterval;
    time_t                                            timeLastFlush;
    time_t                                            timeStepBegin;
    std::reference_wrapper<const Ioss::ParallelUtils> util_;

  public:
    FlushHandler() = delete;

    explicit FlushHandler(const Ioss::ParallelUtils &util)
        : flushInterval(-1), isParallel(true), isFirstOutput(true), flushOnFirstOutput(false),
          timeLastFlushInterval(10), timeStepBeginFlushInterval(10), timeLastFlush(time(nullptr)),
          timeStepBegin(time(nullptr)), util_(util)
    {
    }

    const Ioss::ParallelUtils &util() const { return util_; }

    int getFlushInterval() { return flushInterval; }

    void setFlushInterval(int interval) { flushInterval = interval; }

    bool getIsParallel() const { return isParallel; }

    void setIsParallel(bool parallel) { isParallel = parallel; }

    bool getFlushOnFirstOutput() const { return flushOnFirstOutput; }

    void setFlushOnFirstOutput(bool flush) { flushOnFirstOutput = flush; }

    unsigned int getTimeLastFlushInterval() { return timeLastFlushInterval; }

    void setTimeLastFlushInterval(unsigned int timeInterval)
    {
      timeLastFlushInterval = timeInterval;
    }

    unsigned int getTimeStepBeginFlushInterval() { return timeStepBeginFlushInterval; }

    void setTimeStepBeginFlushInterval(unsigned int timeInterval)
    {
      timeStepBeginFlushInterval = timeInterval;
    }

    void resetTimeStepBegin() { timeStepBegin = time(nullptr); }

    bool doFlush(int state);

    bool allReduce(bool do_flush);

    bool isGreaterThanFlushInterval(time_t lastFlushTime, unsigned int interval);
  };

} // namespace Ioss