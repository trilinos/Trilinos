/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    OF the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_misc_hpp_
#define _test_misc_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

/** Unit-tester for miscellaneous fei:: classes and methods.
*/
class test_misc : public tester {
 public:
  test_misc(MPI_Comm comm);
  virtual ~test_misc();

  const char* getName()
    {
      static const char name[] = "misc fei::";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();

  int serialtest0();

  int serialtest1();

  int serialtest2();

  int serialtest3();
};


#endif // _test_misc_hpp_
