/** \file Zoltan2_DebugManager.cpp

    \brief Debug Manager for Zoltan2

    \author Siva Rajamanickam
*/
// Using the #define for ZOLTAN2_DEBUG makes the print() much faster when there
// is nothing to print. However, it would be nice to allow print() even in the
// release build. We can revisit this if the "if (debug level)" check 
// becomes too expensive, especially while printing something in the innermost
// loop.
//#define ZOLTAN2_DEBUG

#include <Zoltan2_DebugManager.hpp>

namespace Zoltan2
{

DebugManager::DebugManager(Teuchos::RCP<const Teuchos::Comm<int> > comm,
    int debugLevel, std::ostream *os)
    :
    comm_(comm),
    debugLevel_(debugLevel),
    myOS_(os)
{
    myPID_ = (*comm_).getRank();
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
}

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

} //namespace Zoltan2
