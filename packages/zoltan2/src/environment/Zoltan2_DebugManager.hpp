// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_DebugManager.hpp
    \brief Debug output manager for Zoltan2
    \author Siva Rajamanickam
*/
#ifndef ZOLTAN2_DEBUGMANAGER_HPP
#define ZOLTAN2_DEBUGMANAGER_HPP

#include <string>
#include <iostream>

namespace Zoltan2
{
/*! \brief DebugManager contains the methods that perform output of
                      debug and status messages.

   An Environment has a DebugManager.

   Parameters governing debug/status output:

     \li \c debug_level
     \li \c debug_procs
     \li \c debug_output_stream
     \li \c debug_output_file

   For more information see at their definitions in 
   createAllParameters() in Zoltan2_Parameters.cpp.

   If Zoltan2 is compiled with \b Z2_OMIT_ALL_STATUS_MESSAGES, no status
   messages will be displayed and status message code is ifdef'd out.

   \todo For nightly testing, add a build for -DZ2_OMIT_ALL_STATUS_MESSAGES.
*/

class DebugManager
{
    public:

    /*! \brief Constructor
     *   \param rank  the MPI rank of this process.
     *   \param doPrinting  true if this process is one that outputs messages.
     *   \param debugOs      the output stream for debug messages.
     *   \param debugLevel   the highest level of message to print, messages
     *                      that are below this level will be ignored.
     */
    DebugManager ( int rank, bool doPrinting, std::ostream &debugOs, 
      int debugLevel){
      myPID_ = rank;
      iPrint_ = doPrinting;
      myOS_ = &debugOs;
      debugLevel_ = debugLevel;
    }

    /*! \brief Destructor
     */
    virtual ~DebugManager() {};

    /*! \brief Return the output stream for debug/status messages.
     */
    inline std::ostream *getOStream() const { return myOS_; };

    /*! \brief Return the highest level of message that will be printed.
     */
    inline int getDebugLevel() const { return debugLevel_; };

    /*! \brief Print a debug or status message, if this process 
     *         is one of those that is supposed to be doing output.
     *
     *  \param  debugLevel  The level of this message.  If it is higher than
     *                 the debugLevel stated in the constructor, it will be
     *                 ignored.
     *  \param  output      The message.
     */
    inline void print(int debugLevel, const std::string &output){
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
      if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output << std::endl;
#endif
    }

    /*! \brief Print a debug or status message regardless of 
     *   whether this process is one of those that is supposed to be 
     *    doing output.
     *
     *  \param  debugLevel  The level of this message.  If it is higher than
     *                 the debugLevel stated in the constructor, it will be
     *                 ignored.
     *  \param  output      The message.
     */
    inline void printInAllTasks(int debugLevel, const std::string &output){
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
      if (debugLevel <= debugLevel_)
        *myOS_ << output << std::endl;
#endif
    }

    /*! \brief The const char * versions of print functions are needed to 
     *      avoid the expensive conversion to string.
     *
     *  \param  debugLevel  The level of this message.  If it is higher than
     *                 the debugLevel stated in the constructor, it will be
     *                 ignored.
     *  \param  output      The message.
     */
    inline void print(int debugLevel, const char *output){
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
      if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << output << std::endl;
#endif
    }

    /*! \brief The const char * versions of print functions are needed to 
     *      avoid the expensive conversion to string.
     *
     *  \param  debugLevel  The level of this message.  If it is higher than
     *                 the debugLevel stated in the constructor, it will be
     *                 ignored.
     *  \param  output      The message.
     */
    inline void printInAllTasks(int debugLevel, const char *output) {
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output << std::endl;
#endif
    }

    private:

    int myPID_;
    int debugLevel_;
    std::ostream *myOS_;
    bool iPrint_;
};

} //namespace Zoltan2

#endif
