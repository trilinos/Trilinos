// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   PANZER_MEMUTILS_HPP
#define   PANZER_MEMUTILS_HPP

#include <Teuchos_Comm.hpp>
#include <iostream>

namespace panzer
{
  //! The memory usage information.
  struct MemUsage
  {
    size_t currMin, currMax, currTot, peakMin, peakMax, peakTot;
    MemUsage(size_t inCurrMin = 0, size_t inCurrMax = 0, size_t inCurrTot = 0,
      size_t inPeakMin = 0, size_t inPeakMax = 0, size_t inPeakTot = 0)
      :
      currMin(inCurrMin), currMax(inCurrMax), currTot(inCurrTot), peakMin(inPeakMin),
      peakMax(inPeakMax), peakTot(inPeakTot)
    {
    }
    MemUsage operator+(const MemUsage& that) const
    {
      return MemUsage(currMin + that.currMin, currMax + that.currMax,
        currTot + that.currTot, peakMin + that.peakMin, peakMax + that.peakMax,
        peakTot + that.peakTot);
    }
    MemUsage operator-(const MemUsage& that) const
    {
      return MemUsage(currMin - that.currMin, currMax - that.currMax,
        currTot - that.currTot, peakMin - that.peakMin, peakMax - that.peakMax,
        peakTot - that.peakTot);
    }
    MemUsage& operator*=(const size_t& that)
    {
      currMin *= that; peakMin *= that;
      currMax *= that; peakMax *= that;
      currTot *= that; peakTot *= that;
      return *this;
    }
    MemUsage& operator/=(const size_t& that)
    {
      currMin /= that; peakMin /= that;
      currMax /= that; peakMax /= that;
      currTot /= that; peakTot /= that;
      return *this;
    }
  };
  enum MemUsageType {MEM_USAGE_CURRENT, MEM_USAGE_PEAK};

  //! Print memory usage to stream.
  void printMemoryUsage(std::ostream& os, const Teuchos::Comm<int>& comm);
  void printMemoryUsage(std::ostream& os, const Teuchos::Comm<int>& comm,
    const MemUsage& mem);
  void pretty(std::ostream& s, size_t num);

  //! Get memory usage in B.
  MemUsage getMemoryUsage(const Teuchos::Comm<int>& comm);

  //! Returns the peak (maximum so far) resident set size (physical memory use)
  //! measured in bytes, or zero if the value cannot be determined on this OS.
  MemUsage getPeakRSS(const Teuchos::Comm<int>& comm);

  //! Returns the current resident set size (physical memory use) measured in
  //! bytes, or zero if the value cannot be determined on this OS.
  MemUsage getCurrentRSS(const Teuchos::Comm<int>& comm);

  //! Reduce the memory usage over all the processors.
  MemUsage reduceMemUsage(size_t& mem, const Teuchos::Comm<int>& comm,
    const MemUsageType& type);
} // end namespace panzer

#endif // PANZER_MEMUTILS_HPP
