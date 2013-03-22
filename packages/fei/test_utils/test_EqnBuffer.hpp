/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_EqnBuffer_h_
#define _test_EqnBuffer_h_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

class test_EqnBuffer : public tester {
 public:
  test_EqnBuffer(MPI_Comm comm);
  virtual ~test_EqnBuffer();

  const char* getName()
    {
      static const char name[] = "EqnBuffer";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();
};

#endif // _test_EqnBuffer_h_
