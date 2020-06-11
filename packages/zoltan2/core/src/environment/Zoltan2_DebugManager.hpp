// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
