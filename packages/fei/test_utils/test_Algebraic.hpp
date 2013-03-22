/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_Algebraic_h_
#define _test_Algebraic_h_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

/** Unit-tester for using new fei classes with purely algebraic data.
*/
class test_Algebraic : public tester {
 public:
  test_Algebraic(MPI_Comm comm);
  virtual ~test_Algebraic();

  const char* getName()
    {
      static const char name[] = "snl_fei:: classes, purely algebraic data";
      return((const char*)name);
    }

  int runtests();

  int serialtest1();

  int serialtest2();

  int test1();

  int test2();

  int test3();

  int test4();

};

#endif // _test_Algebraic_h_
