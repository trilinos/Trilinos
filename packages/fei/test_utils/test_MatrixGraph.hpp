/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_MatrixGraph_h_
#define _test_MatrixGraph_h_

#include <fei_macros.hpp>

#include <fei_defs.h>
#include <fei_SharedPtr.hpp>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_fwd.hpp>

/** Unit-tester for fei::MatrixGraph.
*/
class test_MatrixGraph : public tester {
 public:
  test_MatrixGraph(MPI_Comm comm);
  virtual ~test_MatrixGraph();

  const char* getName()
    {
      static const char name[] = "fei::MatrixGraph";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int serialtest1();

  int test2();

  int test3();

  int test4();

  int test5();

  static fei::SharedPtr<fei::MatrixGraph> create_MatrixGraph(testData* testdata,
				int localProc, int numProcs,
				bool bothFields, bool nonSymmetric,
				const char* name,
				fei::SharedPtr<fei::VectorSpace> vectorSpacePtr,
				fei::SharedPtr<fei::Factory> factory,
                                const std::string& path,
				bool turnOnDebugOutput=false);
};

#endif // _test_MatrixGraph_h_
