/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    OF the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_Utils_h_
#define _test_Utils_h_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

/** Unit-tester for snl_fei_Utils methods.
*/
class test_Utils : public tester {
 public:
  test_Utils(MPI_Comm comm);
  virtual ~test_Utils();

  const char* getName()
    {
      static const char name[] = "snl_fei_Utils";
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


#endif // _test_Utils_h_
