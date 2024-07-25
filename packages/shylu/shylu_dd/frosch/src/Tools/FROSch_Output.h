// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_OUTPUT_HPP
#define _FROSCH_OUTPUT_HPP

#include <ShyLU_DDFROSch_config.h>

// Width of indent before FROSch output
#ifndef FROSCH_OUTPUT_INDENT
    #define FROSCH_OUTPUT_INDENT 5
#endif

#ifndef FROSCH_ASSERT
    #define FROSCH_ASSERT(COND,MSG) \
    { \
        const bool throw_exception = !(COND); \
        if(throw_exception) { \
            Teuchos::TestForException_incrThrowNumber(); \
            const int throwNumber = Teuchos::TestForException_getThrowNumber(); \
            std::ostringstream omsg; \
            omsg \
                << std::setw(FROSCH_OUTPUT_INDENT) << " " << __FILE__ << ":" << __LINE__ << ":\n\n" \
                << "Throw number = " << throwNumber \
                << "\n\n" \
                << std::setw(FROSCH_OUTPUT_INDENT) << " " << "Throw test that evaluated to true: "#COND \
                << "\n\n" \
                << std::setw(FROSCH_OUTPUT_INDENT) << " " << "[ERROR] " << MSG; \
            const std::string &omsgstr = omsg.str(); \
            TEUCHOS_STORE_STACKTRACE(); \
            Teuchos::TestForException_break(omsgstr, throwNumber); \
            throw std::logic_error(omsgstr); \
        } \
    }
#endif

#ifndef FROSCH_WARNING
    #define FROSCH_WARNING(CLASS,VERBOSE,OUTPUT) if (VERBOSE) std::cerr << std::setw(FROSCH_OUTPUT_INDENT) << " " << "[WARNING] " << CLASS << ": " << OUTPUT << std::endl;
#endif

#ifndef FROSCH_NOTIFICATION
    #define FROSCH_NOTIFICATION(CLASS,VERBOSE,OUTPUT) if (VERBOSE) std::cout << std::setw(FROSCH_OUTPUT_INDENT) << " " << "[NOTIFICATION] " << CLASS << ": " << OUTPUT << std::endl;
#endif

#ifndef FROSCH_TEST_OUTPUT
    #define FROSCH_TEST_OUTPUT(COMM,VERBOSE,OUTPUT) COMM->barrier(); COMM->barrier(); COMM->barrier(); if (VERBOSE) std::cout << std::setw(FROSCH_OUTPUT_INDENT) << " " << OUTPUT << std::endl;
#endif

#endif
