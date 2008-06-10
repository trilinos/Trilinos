/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_FEI_h_
#define _test_FEI_h_

#include <fei_macros.hpp>

#include <fei_defs.h>
#include <fei_mpi.h>

#include <test_utils/tester.hpp>

#include <fei_fwd.hpp>

#include <fei_SharedPtr.hpp>

class test_FEI : public tester {
 public:
  test_FEI(MPI_Comm comm);
  virtual ~test_FEI();

  const char* getName()
    {
      static const char name[] = "FEI";
      return((const char*)name);
    }

  void setFileName(const char* filename)
  { fileName_ = filename; }

  int runtests();

  int test1();

  int test2();

  int test3();

  int test4();

 private:
  fei::SharedPtr<fei::ParameterSet> get_test_parameters();

  std::string fully_qualified_name(const std::string& fileName);

  std::string fileName_;
};


#endif // _test_FEI_h_
