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

#include <string>
#include <iostream>
#include "Teuchos_oblackholestream.hpp"

namespace Zoltan2
{
/*! Zoltan2::DebugManager
    \brief Methods to display debugging statements.

*/

class DebugManager
{
    public:

    DebugManager ( int rank, bool doPrinting, std::ostream &debugOs, int debugLevel){
      myPID_ = rank;
      iPrint_ = doPrinting;
      myOS_ = &debugOs;
      debugLevel_ = debugLevel;
    }

    virtual ~DebugManager() {};

    inline void setRank(int p)
    {
      myPID_ = p;
    }

    inline void setIPrint(bool p)
    {
      iPrint_ = p;
    }

    inline void setOStream(std::ostream &os)
    {
      myOS_ = &os;
    };

    inline void setDebugLevel(int debugLevel) { debugLevel_ = debugLevel; };

    inline std::ostream *getOStream() const { return myOS_; };

    inline int getDebugLevel() const { return debugLevel_; };

    inline void print(int debugLevel, const std::string &output){
//#ifdef ZOLTAN2_DEBUG
      if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output;
//#endif
    }

    inline void printInAllTasks(int debugLevel, const std::string &output){
//#ifdef ZOLTAN2_DEBUG
      if (debugLevel <= debugLevel_)
        *myOS_ << output;
//#endif
    }

    // The const char * versions of print functions are needed to avoid the
    // expensive conversion in code like
    //          print(5, "I am here");
    inline void print(int debugLevel, const char *output){
//#ifdef ZOLTAN2_DEBUG
        if (debugLevel <= debugLevel_ && iPrint_)
            *myOS_ << output;
//#endif
    }

    inline void printInAllTasks(int debugLevel, const char *output) {
//#ifdef ZOLTAN2_DEBUG
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output;
//#endif
    }

    inline void error(const std::string &output) {
      *myOS_ << "PID =" << myPID_ << " " << output;
    }

    private:

    int myPID_;
    int debugLevel_;
    std::ostream *myOS_;
    Teuchos::oblackholestream myBHS_;
    bool iPrint_;
};

} //namespace Zoltan2

#endif
