/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _feitester_hpp_
#define _feitester_hpp_

#include <string>

/**
   A test "harness" for fei implementations. This interface is intended to be
   specialized by more specific test harnesses that exercise various specific
   fei implementations.

   All methods in this interface (except getName()) return int error-codes, and
   it is expected that an error-return of 0 indicates success while any
   non-zero error-return is an indication of failure.

   The methods on this interface should be called in the following
   order:
<ol>
<li>testInitialization()
<li>testLoading()
<li>testSolve()
<li>testCheckResult()
</ol>
*/
class feitester {
 public:
  feitester() : path_() {}
  virtual ~feitester(){}

  /** Method to obtain a name describing this test.
  */
  virtual const char* getName() = 0;

  virtual int testInitialization() = 0;

  virtual int testLoading() = 0;

  virtual int testSolve() = 0;

  virtual int testCheckResult() = 0;

  virtual void dumpMatrixFiles() = 0;

  virtual void setParameter(const char* param) = 0;

  void setPath(const std::string& path)
  { path_ = path; }

  void setPath(const char* path)
  { path_ = path; }

 protected:
  std::string path_;
};

#endif // _feitester_hpp_
