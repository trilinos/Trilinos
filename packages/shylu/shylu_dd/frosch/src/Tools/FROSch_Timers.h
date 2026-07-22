// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
        #define FROSCH_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S))));
    #endif

    #ifndef FROSCH_TIMER_START_LEVELID
        #define FROSCH_TIMER_START_LEVELID(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (Level " + std::to_string(this->LevelID_) + std::string(")"))));
    #endif

    #ifndef FROSCH_TIMER_STOP
        #define FROSCH_TIMER_STOP(A) A.reset();
    #endif

    #ifndef FROSCH_DETAILTIMER_START
        #define FROSCH_DETAILTIMER_START(A,S) FROSCH_TIMER_START(A,S);
    #endif

    #ifndef FROSCH_DETAILTIMER_START_LEVELID
        #define FROSCH_DETAILTIMER_START_LEVELID(A,S) FROSCH_TIMER_START_LEVELID(A,S);
    #endif

    #ifndef FROSCH_DETAILTIMER_STOP
        #define FROSCH_DETAILTIMER_STOP(A) FROSCH_TIMER_STOP(A);
    #endif

    #ifndef FROSCH_TIMER_START_SOLVER
        #define FROSCH_TIMER_START_SOLVER(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (" + this->Description_ + std::string(")"))));
    #endif

    #ifndef FROSCH_TIMER_START_SUBDOMAINSOLVER
        #define FROSCH_TIMER_START_SUBDOMAINSOLVER(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (" + this->Description_ + std::string(")"))));
    #endif
#elif FROSCH_TIMER_DETAILS == 1
    #ifndef FROSCH_TIMER_START
        #define FROSCH_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S))));
    #endif

    #ifndef FROSCH_TIMER_START_LEVELID
        #define FROSCH_TIMER_START_LEVELID(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FROSch: ") + std::string(S) + " (Level " + std::to_string(this->LevelID_) + std::string(")"))));
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

    #ifndef FROSCH_TIMER_START_SOLVER
        #define FROSCH_TIMER_START_SOLVER(A,S)
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

    #ifndef FROSCH_TIMER_START_SOLVER
        #define FROSCH_TIMER_START_SOLVER(A,S)
    #endif

    #ifndef FROSCH_TIMER_START_SUBDOMAINSOLVER
        #define FROSCH_TIMER_START_SUBDOMAINSOLVER(A,S)
    #endif
#endif

#endif
