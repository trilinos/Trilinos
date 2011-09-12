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

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"

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
      myOS_ = Teuchos::rcp(&debugOs);
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
      myOS_ = Teuchos::rcp(&os);
    };

    inline void setDebugLevel(int debugLevel) { debugLevel_ = debugLevel; };

    inline Teuchos::RCP<std::ostream> getOStream() const { return myOS_; };

    inline int getDebugLevel() const { return debugLevel_; };

    inline void print(int debugLevel, const std::string &output);

    inline void printInAllTasks(int debugLevel, const std::string &output);

    // The const char * versions of print functions are needed to avoid the
    // expensive conversion in code like
    //          print(5, "I am here");
    inline void print(int debugLevel, const char *output);

    inline void printInAllTasks(int debugLevel, const char *output);

    inline void error(const std::string &output);

    private:

    int myPID_;
    int debugLevel_;
    Teuchos::RCP<std::ostream> myOS_;
    Teuchos::oblackholestream myBHS_;
    bool iPrint_;
};

} //namespace Zoltan2

#endif // _ZOLTAN2_DEBUG_MANAGER_HPP_
