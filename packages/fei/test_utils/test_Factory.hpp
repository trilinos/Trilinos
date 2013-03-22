/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _test_Factory_hpp_
#define _test_Factory_hpp_

#include <fei_macros.hpp>

#include <test_utils/tester.hpp>
#include <fei_SharedPtr.hpp>

#include <fei_fwd.hpp>

/** Tester for fei::Factory.
    The runtests() method constructs and tests a couple of factory 
    implementations that are contained in the fei source distribution.

    This class can also be used to test an arbitrary fei::Factory
    implementation as follows:

       //construct your specialized fei::Factory:
       fei::Factory* factory = new my_special_factory(...);

       //construct the test_Factory class:
       test_Factory factory_tester(comm);

       //run the test method:
       try {
         factory_tester.factory_test1(factory);
       }
       catch(std::runtime_error& exc) {
         std::cout << "factory test failed, exception: " << exc.what()
                   <<std::endl;
       }

    The factory_test1 method will print a small amount of information to
    cout, describing the tests that it is performing.
*/
class test_Factory : public tester {
 public:
  test_Factory(MPI_Comm comm);
  virtual ~test_Factory();

  const char* getName()
    {
      static const char name[] = "fei::Factory";
      return((const char*)name);
    }

  int runtests();

  void factory_test1(fei::SharedPtr<fei::Factory> factory);
};


#endif // _test_Factory_hpp_
