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
using Teuchos::TimeMonitor;

namespace Amesos{


struct Timers {
  // Timers(){
  //   // initialize Timers
  //   mtxRedistTime_ =
  //     *(TimeMonitor::getNewTimer("Time to redistribute matrix"));
  //   mtxConvTime_ =
  //     *(TimeMonitor::getNewTimer("Time to convert matrix to solver format"));
  //   vecRedistTime_ =
  //     *(TimeMonitor::getNewTimer("Time to redistribute vectors"));
  //   symFactTime_ =
  //     *(TimeMonitor::getNewTimer("Time for symbolic factorization"));
  //   numFactTime_ =
  //     *(TimeMonitor::getNewTimer("Time for numeric factorization"));
  //   solveTime_ =
  //     *(TimeMonitor::getNewTimer("Time for solve"));
  //   overheadTime_ =
  //     *(TimeMonitor::getNewTimer("Amesos2 overhead Time"));
  // }

  Timers()
    : mtxRedistTime_(*(TimeMonitor::getNewTimer("Time to redistribute matrix")))
    , mtxConvTime_(*(TimeMonitor::getNewTimer("Time to convert matrix to solver format")))
    , vecRedistTime_(*(TimeMonitor::getNewTimer("Time to redistribute vectors")))
    , symFactTime_(*(TimeMonitor::getNewTimer("Time for symbolic factorization")))
    , numFactTime_(*(TimeMonitor::getNewTimer("Time for numeric factorization")))
    , solveTime_(*(TimeMonitor::getNewTimer("Time for solve")))
    , overheadTime_(*(TimeMonitor::getNewTimer("Amesos2 overhead Time")))
    {}

  Teuchos::Time mtxRedistTime_;
  Teuchos::Time mtxConvTime_;
  Teuchos::Time vecRedistTime_;
  Teuchos::Time symFactTime_;
  Teuchos::Time numFactTime_;
  Teuchos::Time solveTime_;
  Teuchos::Time overheadTime_;
};				// end struct Amesos2::Timers


} // end namespace Amesos

#endif  // AMESOS2_TIMERS_HPP
