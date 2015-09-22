/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_FEI_Impl_hpp_
#define _test_FEI_Impl_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

/** Unit-tester for fei::FEI_Impl.
*/
class test_FEI_Impl : public tester {
 public:
  test_FEI_Impl(MPI_Comm comm);
  virtual ~test_FEI_Impl();

  const char* getName()
    {
      static const char name[] = "fei::FEI_Impl";
      return((const char*)name);
    }

  int runtests();

  int serialtest1();

  int test1();

  int test2();

  int test3();

  int test4();

  int compareCoefs(int n,
		   const double*const* coefs1,
		   const double*const* coefs2);
};


#endif // _test_FEI_Impl_hpp_

