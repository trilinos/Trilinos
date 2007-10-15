/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_Vector_hpp_
#define _test_Vector_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_SharedPtr.hpp>
#include <fei_Vector.hpp>
#include <fei_Factory.hpp>

/** Unit-tester for fei::Vector.
*/
class test_Vector : public tester {
 public:
  test_Vector(MPI_Comm comm);
  virtual ~test_Vector();

  const char* getName()
    {
      static const char name[] = "fei::Vector";
      return((const char*)name);
    }

  int runtests();

  fei::SharedPtr<fei::Vector> create_vector(fei::SharedPtr<fei::Factory> factory);

  void vector_test1(fei::SharedPtr<fei::Vector> fei_vec);
};


#endif // _test_Vector_hpp_
