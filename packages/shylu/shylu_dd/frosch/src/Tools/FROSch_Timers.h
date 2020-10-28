//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_TIMERS_HPP
#define _FROSCH_TIMERS_HPP

#include <ShyLU_DDFROSch_config.h>

// The level of detail of the timers used in FROSch is determined by FROSCH_TIMER_DETAILS
// Currently, we support no timers (0), selected timers (1), and all timers (2)
#ifndef FROSCH_TIMER_DETAILS
    #define FROSCH_TIMER_DETAILS 1
#endif

#if FROSCH_TIMER_DETAILS > 1
    #ifndef FROSCH_TIMER_START
        #define FROSCH_TIMER_START(A,S) RCP<TimeMonitor> A = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S))));
    #endif

    #ifndef FROSCH_TIMER_START_LEVELID
        #define FROSCH_TIMER_START_LEVELID(A,S) RCP<TimeMonitor> A = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (Level " + std::to_string(this->LevelID_) + std::string(")"))));
    #endif

    #ifndef FROSCH_TIMER_STOP
        #define FROSCH_TIMER_STOP(A) A.reset();
    #endif

    #ifndef FROSCH_DETAILTIMER_START
        #define FROSCH_DETAILTIMER_START(A,S) FROSCH_DETAILTIMER_START(A,S);
    #endif

    #ifndef FROSCH_DETAILTIMER_START_LEVELID
        #define FROSCH_DETAILTIMER_START_LEVELID(A,S) FROSCH_TIMER_START_LEVELID(A,S);
    #endif

    #ifndef FROSCH_DETAILTIMER_STOP
        #define FROSCH_DETAILTIMER_STOP(A,S) FROSCH_TIMER_STOP(A);
    #endif

    #ifndef FROSCH_TIMER_START_SUBDOMAINSOLVER
        #define FROSCH_TIMER_START_SUBDOMAINSOLVER(A,S) RCP<TimeMonitor> A = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (" + this->Description_ + std::string(")"))));
    #endif
#elif FROSCH_TIMER_DETAILS == 1
    #ifndef FROSCH_TIMER_START
        #define FROSCH_TIMER_START(A,S) RCP<TimeMonitor> A = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S))));
    #endif

    #ifndef FROSCH_TIMER_START_LEVELID
        #define FROSCH_TIMER_START_LEVELID(A,S) RCP<TimeMonitor> A = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (Level " + std::to_string(this->LevelID_) + std::string(")"))));
    #endif

    #ifndef FROSCH_TIMER_STOP
        #define FROSCH_TIMER_STOP(A) A.reset();
    #endif

    #ifndef FROSCH_DETAILTIMER_START
        #define FROSCH_DETAILTIMER_START(A,S)
    #endif

    #ifndef FROSCH_DETAILTIMER_START_LEVELID
        #define FROSCH_DETAILTIMER_START_LEVELID(A,S)
    #endif

    #ifndef FROSCH_DETAILTIMER_STOP
        #define FROSCH_DETAILTIMER_STOP(A)
    #endif

    #ifndef FROSCH_TIMER_START_SUBDOMAINSOLVER
        #define FROSCH_TIMER_START_SUBDOMAINSOLVER(A,S)
    #endif
#else
    #ifndef FROSCH_TIMER_START
        #define FROSCH_TIMER_START(A,S)
    #endif

    #ifndef FROSCH_TIMER_START_LEVELID
        #define FROSCH_TIMER_START_LEVELID(A,S)
    #endif

    #ifndef FROSCH_TIMER_STOP
        #define FROSCH_TIMER_STOP(A)
    #endif

    #ifndef FROSCH_DETAILTIMER_START
        #define FROSCH_DETAILTIMER_START(A,S)
    #endif

    #ifndef FROSCH_DETAILTIMER_START_LEVELID
        #define FROSCH_DETAILTIMER_START_LEVELID(A,S)
    #endif

    #ifndef FROSCH_DETAILTIMER_STOP
        #define FROSCH_DETAILTIMER_STOP(A)
    #endif

    #ifndef FROSCH_TIMER_START_SUBDOMAINSOLVER
        #define FROSCH_TIMER_START_SUBDOMAINSOLVER(A,S)
    #endif
#endif

#endif
