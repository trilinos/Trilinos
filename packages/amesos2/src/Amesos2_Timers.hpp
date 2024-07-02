// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

namespace Amesos2 {


struct Timers {

  Timers()
    : mtxRedistTime_(*(Teuchos::TimeMonitor::getNewTimer("Time to redistribute data structures")))
    , mtxConvTime_(*(Teuchos::TimeMonitor::getNewTimer("Time to convert matrix to solver format")))
    , vecRedistTime_(*(Teuchos::TimeMonitor::getNewTimer("Time to redistribute vectors")))
    , vecConvTime_(*(Teuchos::TimeMonitor::getNewTimer("Time to convert vectors to solver format")))
    , preOrderTime_(*(Teuchos::TimeMonitor::getNewTimer("Time for matrix pre-order")))
    , symFactTime_(*(Teuchos::TimeMonitor::getNewTimer("Time for symbolic factorization")))
    , numFactTime_(*(Teuchos::TimeMonitor::getNewTimer("Time for numeric factorization")))
    , solveTime_(*(Teuchos::TimeMonitor::getNewTimer("Time for solve")))
    , totalTime_(*(Teuchos::TimeMonitor::getNewTimer("Total Time in Amesos2 interface")))
    {}

  Teuchos::Time mtxRedistTime_;
  Teuchos::Time mtxConvTime_;
  Teuchos::Time vecRedistTime_;
  Teuchos::Time vecConvTime_;
  Teuchos::Time preOrderTime_;
  Teuchos::Time symFactTime_;
  Teuchos::Time numFactTime_;
  Teuchos::Time solveTime_;
  Teuchos::Time totalTime_;
};


} // end namespace Amesos2

#endif  // AMESOS2_TIMERS_HPP
