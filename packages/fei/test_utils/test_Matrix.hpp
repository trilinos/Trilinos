/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_Matrix_hpp_
#define _test_Matrix_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_fwd.hpp>

/** Unit-tester for snl_fei::Matrix.
*/
class test_Matrix : public tester {
 public:
  test_Matrix(MPI_Comm comm);
  virtual ~test_Matrix();

  const char* getName()
    {
      static const char name[] = "fei::Matrix";
      return((const char*)name);
    }

  int runtests();

  fei::SharedPtr<fei::Matrix> create_matrix(fei::SharedPtr<fei::Factory> factory);

  void matrix_test1(fei::SharedPtr<fei::Matrix> mat);

  int test1();

  int test2();

  int test3();

  int test4();

  int serialtest1();

  int serialtest2();

  int serialtest3();
};


#endif // _test_Matrix_hpp_
