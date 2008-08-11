/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_SNL_FEI_Structure_h_
#define _test_SNL_FEI_Structure_h_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

/** Unit-tester for SNL_FEI_Structure.
*/
class test_SNL_FEI_Structure : public tester {
 public:
  test_SNL_FEI_Structure(MPI_Comm comm);
  virtual ~test_SNL_FEI_Structure();

  const char* getName()
    {
      static const char name[] = "SNL_FEI_Structure";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();
};


#endif // _test_SNL_FEI_Structure_h_
