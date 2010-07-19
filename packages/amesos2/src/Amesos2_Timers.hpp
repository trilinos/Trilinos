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

namespace Amesos{


struct Timers {

  Timers()
    : mtxRedistTime_(*(TimeMonitor::getNewTimer("Time to redistribute data structures")))
    , mtxConvTime_(*(TimeMonitor::getNewTimer("Time to convert matrix to solver format")))
    , vecRedistTime_(*(TimeMonitor::getNewTimer("Time to redistribute vectors")))
    , vecConvTime_(*(TimeMonitor::getNewTimer("Time to convert vectors to solver format")))
    , symFactTime_(*(TimeMonitor::getNewTimer("Time for symbolic factorization")))
    , numFactTime_(*(TimeMonitor::getNewTimer("Time for numeric factorization")))
    , solveTime_(*(TimeMonitor::getNewTimer("Time for solve")))
    , totalTime_(*(TimeMonitor::getNewTimer("Total Time in Amesos2 interface")))
    {}

  Time mtxRedistTime_;
  Time mtxConvTime_;
  Time vecRedistTime_;
  Time vecConvTime_;
  Time symFactTime_;
  Time numFactTime_;
  Time solveTime_;
  Time totalTime_;
};


} // end namespace Amesos

#endif  // AMESOS2_TIMERS_HPP
