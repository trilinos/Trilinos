/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_Reducer_hpp_
#define _test_Reducer_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_fwd.hpp>

/** Unit-tester for fei::Reducer.
*/
class test_Reducer : public tester {
 public:
  test_Reducer(MPI_Comm comm);
  virtual ~test_Reducer();

  const char* getName()
    {
      static const char name[] = "fei::Reducer";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();

  int serialtest1();

  int serialtest2();

  int serialtest3();
};


#endif // _test_Reducer_hpp_
