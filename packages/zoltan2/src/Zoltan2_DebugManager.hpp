/** \file Zoltan2_DebugManager.hpp

    \brief Debug Manager for Zoltan2

    \author Siva Rajamanickam
*/
#ifndef _ZOLTAN2_DEBUG_MANAGER_HPP_
#define _ZOLTAN2_DEBUG_MANAGER_HPP_

// Using the #define for ZOLTAN2_DEBUG makes the print() much faster when there
// is nothing to print. However, it would be nice to allow print() even in the
// release build. We can revisit this if the "if (debug level)" check 
// becomes too expensive, especially while printing something in the innermost
// loop.
//#define ZOLTAN2_DEBUG

#include "Zoltan2_config.h" // Just for HAVE_MPI

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Z2
{

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
//#ifdef ZOLTAN2_DEBUG
    if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output;
//#endif
}

inline void DebugManager::print(int debugLevel, const char *output)
{
//#ifdef ZOLTAN2_DEBUG
    if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output;
//#endif
};

inline void DebugManager::printInAllTasks(int debugLevel,
            const std::string &output)
{
//#ifdef ZOLTAN2_DEBUG
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output;
//#endif
}

inline void DebugManager::printInAllTasks(int debugLevel,
            const char *output)
{
//#ifdef ZOLTAN2_DEBUG
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output;
//#endif
}

// Errors show up for all PIDs, even in release builds
inline void DebugManager::error(const std::string &output)
{
    *myOS_ << "PID =" << myPID_ << " " << output;
}

#endif // _ZOLTAN2_DEBUG_MANAGER_HPP_
