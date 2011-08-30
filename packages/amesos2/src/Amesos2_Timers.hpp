// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_Timers.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Sat Jan 16 09:11:37 2010
  
  \brief  Container class for Timers used with the Amesos2::Solver
          class.
*/

#ifndef AMESOS2_TIMERS_HPP
#define AMESOS2_TIMERS_HPP

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Hashtable.hpp>
using Teuchos::TimeMonitor;
using Teuchos::Time;

namespace Amesos2 {


struct Timers {

  Timers()
    : mtxRedistTime_(*(TimeMonitor::getNewTimer("Time to redistribute data structures")))
    , mtxConvTime_(*(TimeMonitor::getNewTimer("Time to convert matrix to solver format")))
    , vecRedistTime_(*(TimeMonitor::getNewTimer("Time to redistribute vectors")))
    , vecConvTime_(*(TimeMonitor::getNewTimer("Time to convert vectors to solver format")))
    , preOrderTime_(*(TimeMonitor::getNewTimer("Time for matrix pre-order")))
    , symFactTime_(*(TimeMonitor::getNewTimer("Time for symbolic factorization")))
    , numFactTime_(*(TimeMonitor::getNewTimer("Time for numeric factorization")))
    , solveTime_(*(TimeMonitor::getNewTimer("Time for solve")))
    , totalTime_(*(TimeMonitor::getNewTimer("Total Time in Amesos2 interface")))
    {}

  Time mtxRedistTime_;
  Time mtxConvTime_;
  Time vecRedistTime_;
  Time vecConvTime_;
  Time preOrderTime_;
  Time symFactTime_;
  Time numFactTime_;
  Time solveTime_;
  Time totalTime_;
};


} // end namespace Amesos2

#endif  // AMESOS2_TIMERS_HPP
