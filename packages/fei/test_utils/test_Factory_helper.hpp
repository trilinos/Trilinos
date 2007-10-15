#ifndef _test_Factory_helper_h_
#define _test_Factory_helper_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_Factory.hpp>

class test_Factory_helper {
 public:
  test_Factory_helper(){}
  virtual ~test_Factory_helper(){}

  int dyncastMatrix(fei::Matrix* matrix, const char* libname);

  int dyncastVector(fei::Vector* vector, const char* libname);
};

#endif

