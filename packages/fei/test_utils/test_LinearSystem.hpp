/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_LinearSystem_h_
#define _test_LinearSystem_h_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_fwd.hpp>

/** Unit-tester for snl_fei::LinearSystem.
*/
class test_LinearSystem : public tester {
 public:
  test_LinearSystem(MPI_Comm comm);
  virtual ~test_LinearSystem();

  const char* getName()
    {
      static const char name[] = "snl_fei::LinearSystem";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();

  int test5();

};


#endif // _test_LinearSystem_h_
