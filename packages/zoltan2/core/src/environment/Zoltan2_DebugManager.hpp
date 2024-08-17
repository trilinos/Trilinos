// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_DebugManager.hpp
    \brief Debug output manager for Zoltan2
    \author Siva Rajamanickam
*/
#ifndef ZOLTAN2_DEBUGMANAGER_HPP
#define ZOLTAN2_DEBUGMANAGER_HPP

#include <Zoltan2_Parameters.hpp>
#include <string>
#include <iostream>
#include <fstream>

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

    /*! \brief Constructor for output to an ofstream.
     *   \param rank  the MPI rank of this process.
     *   \param doPrinting  true if this process is one that outputs messages.
     *   \param debugOs      the output stream for debug messages.
     *   \param debugLevel   the highest level of message to print, messages
     *                      that are below this level will be ignored.
     *
     * Different constructor for file output so we can close the file
     * in the destructor.
     */
    DebugManager ( int rank, bool doPrinting, std::ofstream &debugOs, 
      MessageOutputLevel debugLevel) : myPID_(rank), debugLevel_(debugLevel),
        myOS_(static_cast<std::ostream *>(&debugOs)), fileOS_(&debugOs), 
        iPrint_(doPrinting) {}

    /*! \brief Constructor for output to an iostream.
     *   \param rank  the MPI rank of this process.
     *   \param doPrinting  true if this process is one that outputs messages.
     *   \param debugOs      the output stream for debug messages.
     *   \param debugLevel   the highest level of message to print, messages
     *                      that are below this level will be ignored.
     */
    DebugManager ( int rank, bool doPrinting, std::ostream &debugOs, 
      MessageOutputLevel debugLevel) : myPID_(rank), debugLevel_(debugLevel),
        myOS_(&debugOs), fileOS_(NULL), iPrint_(doPrinting) {}

    /*! \brief Destructor
     */
    virtual ~DebugManager()
    {
      myOS_->flush();
      if (fileOS_)
        fileOS_->close();
    }

    /*! \brief Return the output stream for debug/status messages.
     */
    inline std::ostream *getOStream() const { return myOS_; };

    /*! \brief Return the highest level of message that will be printed.
     */
    inline MessageOutputLevel getDebugLevel() const { return debugLevel_; };

    /*! \brief Print a debug or status message, if this process 
     *         is one of those that is supposed to be doing output.
     *
     *  \param  debugLevel  The level of this message.  If it is higher than
     *                 the debugLevel stated in the constructor, it will be
     *                 ignored.
     *  \param  output      The message.
     */
    inline void print(MessageOutputLevel debugLevel, const std::string &output){
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
      if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << myPID_ << ": " << output << std::endl;
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
    inline void printInAllTasks(MessageOutputLevel debugLevel, const std::string &output){
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
      if (debugLevel <= debugLevel_)
        *myOS_ << myPID_ << ": " << output << std::endl;
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
    inline void print(MessageOutputLevel debugLevel, const char *output){
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
      if (debugLevel <= debugLevel_ && iPrint_)
        *myOS_ << myPID_ << ": " << output << std::endl;
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
    inline void printInAllTasks(MessageOutputLevel debugLevel, const char *output) {
#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
    if (debugLevel <= debugLevel_)
        *myOS_ << "PID =" << myPID_ << " " << output << std::endl;
#endif
    }

    private:

    int myPID_;
    MessageOutputLevel debugLevel_;
    std::ostream *myOS_;
    std::ofstream *fileOS_;
    bool iPrint_;
};

} //namespace Zoltan2

#endif
