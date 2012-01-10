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

#include <fei_macros.hpp>
#include <fei_mpi.h>

#include <test_utils/test_Factory.hpp>

#include <test_utils/LibraryFactory.hpp>

#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>

#include <fei_Factory.hpp>
#include <snl_fei_Factory.hpp>

#include <test_utils/test_Factory_helper.hpp>

#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
#endif

#include <fei_Factory_Trilinos.hpp>

#undef fei_file
#define fei_file "test_Factory.cpp"
#include <fei_ErrMacros.hpp>

test_Factory::test_Factory(MPI_Comm comm)
  : tester(comm)
{
}

test_Factory::~test_Factory()
{
}

int test_Factory::runtests()
{
//////////// Factory_Trilinos test ///////////////////////////////
  {
  if (localProc_==0) FEI_COUT << "constructing Factory_Trilinos...";

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm_)); 

  if (localProc_==0) FEI_COUT << "ok" << FEI_ENDL;

  factory_test1(factory);

  if (localProc_==0) FEI_COUT << "testing fei::Factory::clone..." << FEI_ENDL;

  fei::SharedPtr<fei::Factory> clone(factory->clone());

  factory_test1(clone);

  FEI_COUT << FEI_ENDL;
  }

#ifdef HAVE_FEI_AZTECOO
//////////////// snl_fei::Factory(Aztec) test //////////////////////
  {
  if (localProc_==0) FEI_COUT << "constructing snl_fei::Factory(Aztec)...";

  fei::SharedPtr<LinearSystemCore> az_lsc(new fei_trilinos::Aztec_LinSysCore(comm_));

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, az_lsc));

  if (localProc_==0) FEI_COUT << "ok" << FEI_ENDL;

  factory_test1(factory);

  if (localProc_==0) FEI_COUT << "testing fei::Factory::clone..." << FEI_ENDL;

  fei::SharedPtr<fei::Factory> clone(factory->clone());

  factory_test1(clone);

  FEI_COUT << FEI_ENDL;
  }

#endif
  return(0);
}

void test_Factory::factory_test1(fei::SharedPtr<fei::Factory> factory)
{
  if (localProc_==0) FEI_COUT << "  testing factory->createVectorSpace...";

  fei::SharedPtr<fei::VectorSpace> vecspace =
    factory->createVectorSpace(comm_, "dummy_Name");

  if (vecspace.get() == 0) {
    FEI_COUT << "no"<<FEI_ENDL;
    throw std::runtime_error("factory failed to create a fei::VectorSpace");
  }

  //do an extremely simple test to make sure the vector-space
  //is 'alive'.
  int fieldID = 0;
  int fieldSize = 3;
  vecspace->defineFields(1, &fieldID, &fieldSize);

  if (vecspace->getNumFields() != 1) {
    FEI_COUT << "no"<<FEI_ENDL;
    throw std::runtime_error("vecspace->defineFields/getNumFields failed.");
  }

  if (localProc_==0) FEI_COUT << "ok"<<FEI_ENDL;

  if (localProc_==0) FEI_COUT << "  testing factory->createFEI...";

  fei::SharedPtr<FEI> fei = factory->createFEI(comm_);

  //again, do a simple test to make sure the FEI instance is alive...

  int err = fei->initFields(1, &fieldSize, &fieldID);
  if (err != 0) {
    FEI_COUT << "failed"<<FEI_ENDL;
    throw std::runtime_error("fei->initFields() failed.");
  }

  int testFieldSize = -1;
  err = fei->getFieldSize(fieldID, testFieldSize);
  if (err != 0 || testFieldSize != fieldSize) {
    FEI_COUT << "failed"<<FEI_ENDL;
    throw std::runtime_error("fei->getFieldSize() failed.");
  }

  if (localProc_==0) FEI_COUT << "ok"<<FEI_ENDL;
}
