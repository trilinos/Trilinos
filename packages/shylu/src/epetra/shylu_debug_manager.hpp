
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/** \file shylu_debug_manager.hpp

    \brief Debug Manager for ShyLU

    \author Siva Rajamanickam
*/
#ifndef SHYLU_DEBUG_MANAGER_HPP
#define SHYLU_DEBUG_MANAGER_HPP

// Using the #define for SHYLU_DEBUG makes the print() much faster when there
// is nothing to print. However, it would be nice to allow print() even in the
// release build. We can revisit this if the "if (debug level)" check 
// becomes too expensive, especially while printing something in the innermost
// loop.
//#define SHYLU_DEBUG

#include "Isorropia_config.h" // Just for HAVE_MPI

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

class DebugManager
{
    public:

    DebugManager (int debugLevel = 0,
     const Teuchos::RCP< std::ostream > &os = Teuchos::rcp(&std::cout,false)
     );

    virtual ~DebugManager() {};

    inline void setOStream(const Teuchos::RCP<std::ostream> &os)
    {
        myOS_ = os;
    };

    inline void setDebugLevel(int debugLevel) { debugLevel_ = debugLevel; };

    // TODO: Do we need this ?
    inline std::ostream& stream()
    {
        if ( debugLevel_ && iPrint_ )
            return *myOS_;
        else
            return myBHS_;
    }

    inline Teuchos::RCP<std::ostream> getOStream() { return myOS_; };

    inline void print(int debugLevel, const std::string &output);

    inline void printInAllTasks(int debugLevel, const std::string &output);

    // The const char * versions of print functions are needed to avoid the
    // expensive conversion in code like
    //          print(5, "I am here");
    inline void print(int debugLevel, const char *output);

    inline void printInAllTasks(int debugLevel, const char *output);

    inline void error(const std::string &output);

    private:

    int debugLevel_;
    Teuchos::RCP<std::ostream> myOS_;
    Teuchos::oblackholestream myBHS_;
    bool iPrint_;
    int myPID_;
};


DebugManager::DebugManager(int debugLevel, const Teuchos::RCP<std::ostream> &os)
    :
    debugLevel_(debugLevel),
    myOS_(os)
{
#ifdef HAVE_MPI
        int mpiStarted = 0;
        MPI_Initialized(&mpiStarted);
        if (mpiStarted) MPI_Comm_rank(MPI_COMM_WORLD, &myPID_);
        else myPID_=0;
#else
        myPID_ = 0;
#endif
        iPrint_ = (myPID_ == 0);
}

inline void DebugManager::print(int debugLevel, const std::string &output)
{
//#ifdef SHYLU_DEBUG
    if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output;
//#endif
}

inline void DebugManager::print(int debugLevel, const char *output)
{
//#ifdef SHYLU_DEBUG
    if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output;
//#endif
}

inline void DebugManager::printInAllTasks(int debugLevel,
            const std::string &output)
{
//#ifdef SHYLU_DEBUG
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output;
//#endif
}

inline void DebugManager::printInAllTasks(int debugLevel,
            const char *output)
{
//#ifdef SHYLU_DEBUG
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output;
//#endif
}

// Errors show up for all PIDs, even in release builds
inline void DebugManager::error(const std::string &output)
{
    *myOS_ << "PID =" << myPID_ << " " << output;
}

#endif // SHYLU_DEBUG_MANAGER_HPP
