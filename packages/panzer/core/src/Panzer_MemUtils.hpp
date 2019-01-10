// @HEADER                                                                                                                                    
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
