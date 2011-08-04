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

namespace Z2
{
/*! Zoltan2::DebugManager
    \brief Methods to display debugging statements.

*/

class DebugManager
{
    public:

    DebugManager ( 
     Teuchos::RCP<const Teuchos::Comm<int> > comm =
                 Teuchos::DefaultComm<int>::getComm(),
     int debugLevel = 0,
     std::ostream *os = &std::cout
     );

    virtual ~DebugManager() {};

    inline void setComm(const Teuchos::RCP<Teuchos::Comm<int> > &comm)
    {
      comm_ = comm;
    }

    inline void setOStream(std::ostream *os)
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

    inline Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
        return comm_;
    };

    inline std::ostream *getOStream() const { return myOS_; };

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

    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    int debugLevel_;
    std::ostream *myOS_;
    Teuchos::oblackholestream myBHS_;
    bool iPrint_;
    int myPID_;
};

} //namespace Z2

#endif // _ZOLTAN2_DEBUG_MANAGER_HPP_
