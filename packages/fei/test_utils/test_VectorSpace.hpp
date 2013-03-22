/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_VectorSpace_h_
#define _test_VectorSpace_h_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_fwd.hpp>

/** Unit-tester for fei::VectorSpace.
*/
class test_VectorSpace : public tester {
 public:
  test_VectorSpace(MPI_Comm comm);
  virtual ~test_VectorSpace();

  const char* getName()
    {
      static const char name[] = "fei::VectorSpace";
      return((const char*)name);
    }

  int runtests();

  int test0();

  int test1();

  int test2();

  int test3();

  int test4();

  static fei::SharedPtr<fei::VectorSpace> create_VectorSpace(MPI_Comm comm);

  static fei::SharedPtr<fei::VectorSpace> create_VectorSpace(MPI_Comm comm,
				testData* testdata,
				int localProc, int numProcs,
				bool defineBothFields, bool initSolnBothFields,
				const char* name,
				fei::SharedPtr<fei::Factory> factory,
				bool turnOnDebugOutput=false);
};

#endif // _test_VectorSpace_h_
