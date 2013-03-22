/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_LogManager_hpp_
#define _fei_LogManager_hpp_

#include <fei_fwd.hpp>

#include <string>
#include <vector>

namespace fei {

/** Singleton class to manage attributes controlling the type and
   amount of data that should be written to the fei log file.
*/
class LogManager {
 public:
  /** destructor */
  virtual ~LogManager();

  /** Accessor for the one-and-only instance of LogManager.
      Constructs a LogManager instance on the first call, returns
      that same instance on the first and all subsequent calls.
  */
  static LogManager& getLogManager();

  /** Query output-level. Result is an enumeration. The enumeration is
   defined in fei_fwd.hpp. */
  OutputLevel getOutputLevel();

  /** Set output-level, using an enumeration. The enumeration is
   defined in fei_fwd.hpp. */
  void setOutputLevel(OutputLevel olevel);

  /** Set output-level, using a string. Valid values are strings that
   match the names of the enumeration values. e.g., "MATRIX_FILES", etc.
   */
  void setOutputLevel(const char* olevel);

  /** Specify path where debug-log files should be written. */
  void setOutputPath(const std::string& opath);

  /** Query for string specifying path to where debug-log files should
      be written. */
  const std::string& getOutputPath();

  /** Set numProcs and localProc (which will be used in the log-file-name).
  */
  void setNumProcs(int nprocs, int localproc);

 private:
  /** constructor */
  LogManager();

  OutputLevel output_level_;
  std::string output_path_;
  int numProcs_;
  int localProc_;
}; //class LogManager
}//namespace fei
#endif

