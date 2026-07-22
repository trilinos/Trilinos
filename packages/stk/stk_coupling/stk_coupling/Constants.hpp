/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_CONSTANTS_HPP
#define STK_COUPLING_CONSTANTS_HPP

#include <string>

namespace stk
{
namespace coupling
{

enum SyncMode {
  Minimum = 0,
  Receive,
  Send,
  Any
};

//BEGINCouplingReservedNames
static const std::string AppName = "Application Name";
static const std::string TimeSyncMode = "Time Sync Mode";
static const std::string InitialTime = "Initial Time";
static const std::string CurrentTime = "Current Time";
static const std::string TimeStep = "Time Step";
static const std::string FinalTime = "Final Time";
static const std::string IsFinished = "Is Finished";
static const std::string SuccessFlag = "Is Successful";
static const std::string GlobalVars = "Global Vars";
static const std::string CouplingVersion = "CouplingVersion";
static const std::string ConvergenceStatus = "iteration_convergence_status";
static const std::string StepContinuationStatus = "solve_step_continuation_status";
//ENDCouplingReservedNames

}
}

#endif /* STK_COUPLING_CONSTANTS_HPP */
