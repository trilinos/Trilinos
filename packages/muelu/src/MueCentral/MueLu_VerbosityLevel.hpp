#ifndef MUELU_VERBOSITYLEVEL_HPP
#define MUELU_VERBOSITYLEVEL_HPP

#include <Teuchos_VerbosityLevel.hpp>

namespace MueLu {

  enum MsgType 
    {
      Errors       = 0x0000001, //!< Errors

      Warnings0    = 0x0000010, //!< Important warning messages (one line)
      Warnings00   = 0x0000020, //!< Important warning messages (more verbose)
      Warnings1    = 0x0000040, //!< Additional warnings
      PerfWarnings = 0x0000080, //!< Performance warnings
   
      Runtime0     = 0x0000100, //!< One-liner description of what is happening
      Runtime1     = 0x0000200, //!< description of what is happening (more verbose)

      Parameters0  = 0x0001000, //!< Print class parameters
      Parameters1  = 0x0002000, //!< Print class parameters (more parameters, more verbose)

      Statistics0  = 0x0010000, //!< Print statistics that do not involve significant additional computation
      Statistics1  = 0x0020000, //!< Print more statistics

      Timings0     = 0x0100000, //!< Print high level timing information
      Timings1     = 0x0200000, //!< Print detailed timing information

      External     = 0x1000000, //!< Print external lib objects
      Debug        = 0x2000000, //!< Print additional debugging information

      // Predefined combinaisons of MsgType
      // Can be used in user code or examples. Do not used as input parameters of IsPrint() or GetOStream().
      Warnings     = Warnings0 + Warnings00 + Warnings1 + PerfWarnings, //!< Print all warning messages
      Runtime      = Runtime0 + Runtime1,                               //!< Print description of what is going on
      Parameters   = Parameters0 + Parameters1,                         //!< Print parameters
      Statistics   = Statistics0 + Statistics1,                         //!< Print all statistics
      Timings      = Timings0 + Timings1,                               //!< Print all timing information

      //
      None    = 0,
      Low     = Errors + Statistics0 + Timings0,
      Medium  = Errors + Warnings0 + Runtime0 + Parameters0 + Statistics0 + Timings0,
      High    = Errors + Warnings  + Runtime  + Parameters  + Statistics  + Timings,
#ifdef HAVE_MUELU_DEBUG
      Extreme = Errors + Warnings  + Runtime  + Parameters  + Statistics  + Timings + External + Debug,
#else
      Extreme = Errors + Warnings  + Runtime  + Parameters  + Statistics  + Timings + External,
#endif
      Default = Low
    };

  //!
  typedef int VerbLevel;    

  //!
  VerbLevel toMueLuVerbLevel(const Teuchos::EVerbosityLevel verbLevel);

} // namespace MueLu

#endif
