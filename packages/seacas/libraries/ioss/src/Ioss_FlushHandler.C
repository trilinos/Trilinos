// Copyright(C) 1999-2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_FlushHandler.h"
#include "Ioss_SerializeIO.h"

namespace Ioss {
  bool FlushHandler::doFlush(int state)
  {
    bool do_flush = false;
    if (flushOnFirstOutput && isFirstOutput) {
      do_flush = true;
    }
    else if (flushInterval == 1) {
      do_flush = true;
    }
    else if (flushInterval < 0) {
      do_flush = isGreaterThanFlushInterval(timeLastFlush, timeLastFlushInterval);
    }
    else if (flushInterval > 1) {
      if (state % flushInterval == 0) {
        do_flush = true;
      }
    }

    if (flushInterval != 0 && !do_flush) {
      // One last check -- if output took more than timeStepBeginFlushInterval (seconds)
      // then flush since the relative flush cost is outweighted by the time
      // it took to do the output (Basically, we have a lot of data being output...)
      do_flush = isGreaterThanFlushInterval(timeStepBegin, timeStepBeginFlushInterval);
    }

    if (isFirstOutput) {
      isFirstOutput = false;
    }
    return allReduce(do_flush);
  }

  bool FlushHandler::isGreaterThanFlushInterval(time_t lastFlushTime, unsigned int interval)
  {
    time_t cur_time = time(nullptr);
    bool   do_flush = false;

    if (std::difftime(cur_time, lastFlushTime) >= interval) {
      timeLastFlush = cur_time;
      do_flush      = true;
    }
    return do_flush;
  }

  bool FlushHandler::allReduce(bool do_flush)
  {
    int iflush = do_flush ? 1 : 0;
#ifdef SEACAS_HAVE_MPI
    if (isParallel) {
      if (Ioss::SerializeIO::isEnabled()) {
        util().broadcast(iflush);
      }
      else {
        iflush = util().global_minmax(iflush, Ioss::ParallelUtils::DO_MAX);
      }
    }
#endif
    return iflush == 1;
  }

} // namespace Ioss