/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
