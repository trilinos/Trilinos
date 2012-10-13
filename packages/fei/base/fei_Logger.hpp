/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Logger_hpp_
#define _fei_Logger_hpp_

#include <fei_fwd.hpp>
#include <fei_iosfwd.hpp>
#include <set>

namespace fei {
/** Class to be inherited by fei classes that wish to write
    to the fei debug-log file. */
class Logger {
 public:
  /** constructor */
  Logger();
  /** destructor */
  virtual ~Logger();

  /** set specified output-level. */
  void setOutputLevel(OutputLevel olevel);

  void addLogID(int ID);
  void addLogEqn(int eqn);

  bool isLogID(int ID);
  bool isLogEqn(int eqn);

  std::set<int>& getLogIDs();
  std::set<int>& getLogEqns();

 protected:
  /** output level
    Note that the OutputLevel enum is defined in fei_fwd.hpp.
  */
  OutputLevel output_level_;
  /** output stream */
  FEI_OSTREAM* output_stream_;

  std::set<int> logIDs_;
  std::set<int> logEqns_;

 private:
  Logger(const Logger& src);
  Logger& operator=(const Logger& src);
};//class Logger
}//namespace fei
#endif

