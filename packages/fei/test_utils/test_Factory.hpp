/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

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
