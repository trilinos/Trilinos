/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _tester_hpp_
#define _tester_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>

#include <string>

/**
   A simple, general test "harness". This interface is intended to be specialized
   by more specific test harnesses that run tests on particular classes.

*/
class tester {
 public:
  tester(MPI_Comm comm);

  virtual ~tester();

  /** A name describing this test.
  */
  virtual const char* getName() = 0;

  virtual int runtests() = 0;

  void setPath(const std::string& path);

 protected:
  MPI_Comm comm_;
  int numProcs_, localProc_;
  std::string path_;
};

#endif // _tester_hpp_
