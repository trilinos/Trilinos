// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_VERBOSITYLEVEL_HPP
#define MUELU_VERBOSITYLEVEL_HPP

#include <Teuchos_VerbosityLevel.hpp>

namespace MueLu {

  enum MsgType
    {
      Errors          = 0x0000001, //!< Errors

      Warnings0       = 0x0000010, //!< Important warning messages (one line)
      Warnings00      = 0x0000020, //!< Important warning messages (more verbose)
      Warnings1       = 0x0000040, //!< Additional warnings
      PerfWarnings    = 0x0000080, //!< Performance warnings

      Runtime0        = 0x0000100, //!< One-liner description of what is happening
      Runtime1        = 0x0000200, //!< Description of what is happening (more verbose)
      RuntimeTimings  = 0x0000400, //!< Timers that are enabled (using Timings0/Timings1) will be printed during the execution
      NoTimeReport    = 0x0000800, //!< By default, enabled timers appears in the teuchos time monitor summary. Use this option if you do not want to record timing information.

      Parameters0     = 0x0001000, //!< Print class parameters
      Parameters1     = 0x0002000, //!< Print class parameters (more parameters, more verbose)

      Statistics0     = 0x0010000, //!< Print statistics that do not involve significant additional computation
      Statistics1     = 0x0020000, //!< Print more statistics

      Timings0        = 0x0100000, //!< High level timing information (use Teuchos::TimeMonitor::summarize() to print)
      Timings1        = 0x0200000, //!< Detailed timing information   (use Teuchos::TimeMonitor::summarize() to print)
      TimingsByLevel  = 0x0400000, //!< Record timing information level by level. Must be used in combinaison with Timings0/Timings1

      External        = 0x1000000, //!< Print external lib objects
      Debug           = 0x2000000, //!< Print additional debugging information

      // Predefined combinations of MsgType
      // Can be used in user code or examples. Do not used as input parameters of IsPrint() or GetOStream().
      Warnings        = Warnings0 | Warnings00 | Warnings1 | PerfWarnings, //!< Print all warning messages
      Runtime         = Runtime0 | Runtime1,                               //!< Print description of what is going on
      Parameters      = Parameters0 | Parameters1,                         //!< Print parameters
      Statistics      = Statistics0 | Statistics1,                         //!< Print all statistics
      Timings         = Timings0 | Timings1 | TimingsByLevel,              //!< Print all timing information

      //
      None    = 0,
      Low     = Errors | Statistics0 | Timings0,
      Medium  = Errors | Warnings0 | Runtime0 | Parameters0 | Statistics0 | Timings0,
      High    = Errors | Warnings  | Runtime  | Parameters  | Statistics  | Timings,
#ifdef HAVE_MUELU_DEBUG
      Extreme = Errors | Warnings  | Runtime  | Parameters  | Statistics  | Timings | External | Debug,
#else
      Extreme = Errors | Warnings  | Runtime  | Parameters  | Statistics  | Timings | External,
#endif
      Default = High, // This is the default of print() methods. For VerboseObject, another default is set by VerboseObject::globalVerbLevel_ // TODO: move it to the VerboseObject class

      NotSpecified = -1
    };

  //!
  typedef int VerbLevel;

  //!
  VerbLevel toMueLuVerbLevel(const Teuchos::EVerbosityLevel verbLevel);

} // namespace MueLu

#endif
