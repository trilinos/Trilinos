/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_FillableVec_hpp_
#define _test_FillableVec_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

/** Unit-tester for fei::FillableVec
*/
class test_FillableVec : public tester{
 public:
  test_FillableVec(MPI_Comm comm);
  virtual ~test_FillableVec();

  const char* getName()
    {
      static const char name[] = "fei::FillableVec";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();

  int test5();

  int test6();

  int test7();

  int test8();

};

#endif // _test_FillableVec_hpp_
