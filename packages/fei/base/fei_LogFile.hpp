/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_LogFile_hpp_
#define _fei_LogFile_hpp_

#include "fei_iosfwd.hpp"

namespace fei {

/** Singleton class to manage (open, close, etc.) the one-and-only
   fei log file.
*/
class LogFile {
 public:
  /** destructor */
  virtual ~LogFile();

  /** Open a log-file ostream. If one is already open, it is closed
    before the new one is opened.

    The file name is 'fei_log.counter.nprocs.localproc', where
    counter is the number of times this function has been called,
    and nprocs and localproc are specified in the arguments.

    @param path Path, not including file-name, to log-file.
    @param nprocs Number of processors.
    @param localproc Rank of local processor.
  */
  void openOutputStream(const char* path=NULL,
                        int nprocs=1,
                        int localproc=0);

  /** Query for the log-file ostream. */
  FEI_OSTREAM* getOutputStream();

  /** Destroy the log-file ostream (closes the file).
  */
  void closeOutputStream();

  /** Accessor for the one-and-only instance of LogFile.
      Constructs a LogFile instance on the first call, returns
      that same instance on the first and all subsequent calls.
  */
  static LogFile& getLogFile();

 private:
  /** constructor */
  LogFile();

  FEI_OSTREAM* output_stream_;
  unsigned counter_;
}; //class LogFile
}//namespace fei
#endif
