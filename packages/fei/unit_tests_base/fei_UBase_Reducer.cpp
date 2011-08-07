
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_mpi.h>

#if defined(HAVE_FEI_EPETRA) && defined(HAVE_FEI_AZTECOO)
#include <test_utils/LibraryFactory.hpp>
#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>
#endif

#include <fei_MatrixGraph_Impl2.hpp>
#include <fei_Reducer.hpp>
#include <fei_Factory.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_impl_utils.hpp>
#include <fei_defs.h>
#include <snl_fei_Factory.hpp>

#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
#endif

#include <fei_Factory_Trilinos.hpp>

#ifdef HAVE_FEI_FETI
#include <FETI_DP_FiniteElementData.h>
#endif

#undef fei_file
#define fei_file "fei_unit_Reducer.cpp"
#include <fei_ErrMacros.hpp>

#include <vector>
#include <cmath>

namespace {

TEUCHOS_UNIT_TEST(Reducer, reducer_unit1)
{
  //define equation-space to be 0 .. 4
  std::vector<int> eqns(5);
  for(unsigned i=0; i<5; ++i) {
    eqns[i] = i;
  }

  fei::SharedPtr<fei::FillableMat> D(new fei::FillableMat);

  //create slave 1 = 0.5*(2+3)

  D->putCoef(1, 2, 0.5);
  D->putCoef(1, 3, 0.5);

  fei::SharedPtr<fei::CSVec> g;
  fei::Reducer reducer(D, g, MPI_COMM_WORLD);

  reducer.setLocalUnreducedEqns(eqns);

  TEUCHOS_TEST_EQUALITY(reducer.isSlaveEqn(1), true, out, success);

  bool exception_caught = false;
  //calling translateToReducedEqn with a slave eqn should cause reducer
  //to throw an exception.
  try {
    reducer.translateToReducedEqn(1);
  }
  catch(...) {
    exception_caught = true;
  }

  TEUCHOS_TEST_EQUALITY(exception_caught, true, out, success);

  int reducedEqn = reducer.translateToReducedEqn(3);
  TEUCHOS_TEST_EQUALITY(reducedEqn, 2, out, success);
}

void fill_matrices(fei::FillableMat& Kii, fei::FillableMat& Kid,
                   fei::FillableMat& Kdi, fei::FillableMat& Kdd)
{
  //make a tri-diagonal matrix:
  //    cols   0  1  2  3  4
  //        ---------------
  // row 0: |  2 -1
  // row 1: | -1  2 -1
  // row 2: |    -1  2 -1
  // row 3: |       -1  2 -1
  // row 4: |          -1  2
  //
  //Partition it like this:
  //| Kii Kid |
  //| Kdi Kdd |
  //
  //such that Kii is 4x4.
  // i.e.,
  //Kii will contain [0:3,0:3]
  //Kid will contain [  4,0:3]
  //Kdi will contain [0:3,  4]
  //Kdd will contain [  4,  4]

  Kii.putCoef(0,0, 2.0);  Kii.putCoef(0,1, -1.0);
  Kii.putCoef(1,0, -1.0); Kii.putCoef(1,1, 2.0);  Kii.putCoef(1,2, -1.0);
  Kii.putCoef(2,1, -1.0); Kii.putCoef(2,2, 2.0);  Kii.putCoef(2,3, -1.0);
                          Kii.putCoef(3,2, -1.0); Kii.putCoef(3,3, 2.0);

  Kid.putCoef(3,4, -1.0);

  Kdi.putCoef(4,3, -1.0);

  Kdd.putCoef(4,4, 2.0);
}

void addFillableMatToReducer(fei::FillableMat& mat, fei::Reducer& reducer,
                             fei::Matrix& feimat)
{
  fei::FillableMat::iterator
    iter = mat.begin(), iter_end = mat.end();
  for(; iter!=iter_end; ++iter) {
    int row = iter->first;
    fei::FillableVec* rowvec = iter->second;
    fei::FillableVec::iterator
      viter = rowvec->begin(), viter_end = rowvec->end();
    for(; viter!=viter_end; ++viter) {
      int col = viter->first;
      double coef = viter->second;
      double* coefPtr = &coef;
      reducer.addMatrixValues(1, &row, 1, &col, &coefPtr, true, feimat, FEI_DENSE_ROW);
    }
  }
}

#if defined(HAVE_FEI_EPETRA) && defined(HAVE_FEI_AZTECOO)

TEUCHOS_UNIT_TEST(Reducer, test_Reducer_test1)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs > 1) return;

  fei::SharedPtr<fei::FillableMat> D(new fei::FillableMat);

  //create slave 4 = 0.5*(2+3)

  D->putCoef(4, 2, 0.5);
  D->putCoef(4, 3, 0.5);

  fei::SharedPtr<fei::CSVec> g;
  fei::Reducer reducer(D, g, comm);

  //define equation-space to be 0 .. 4
  std::vector<int> eqns(5);
  for(unsigned i=0; i<5; ++i) {
    eqns[i] = i;
  }

  reducer.setLocalUnreducedEqns(eqns);

  fei::FillableMat Kii, Kid, Kdi, Kdd;

  fill_matrices(Kii, Kid, Kdi, Kdd);

  //Now form the reduced matrix Kr = Kii + Kid*D + D^T*Kdi + D^T*Kdd*D
  fei::CSRMat csrD(*D);
  fei::CSRMat Kr, tmpMat1, tmpMat2;
  fei::CSRMat csrKdi(Kdi), csrKid(Kid), csrKdd(Kdd);

  Kr = Kii;

  fei::multiply_CSRMat_CSRMat(csrKid, csrD, tmpMat1);
  fei::multiply_trans_CSRMat_CSRMat(csrD, csrKdi, tmpMat2);

  Kr += tmpMat1;
  Kr += tmpMat2;

  fei::multiply_CSRMat_CSRMat(csrKdd, csrD, tmpMat1);
  fei::multiply_trans_CSRMat_CSRMat(csrD, tmpMat1, tmpMat2);

  Kr += tmpMat2;

  fei::impl_utils::translate_to_reduced_eqns(reducer, Kr);

  fei::SharedPtr<fei::Factory> factory;
  try {
    factory = fei::create_fei_Factory(comm, "Trilinos");
  }
  catch(...) {
    FEI_COUT << "\ncouldn't create Trilinos factory."<<FEI_ENDL;
    return;
  }

  fei::SharedPtr<fei::VectorSpace> vspace =
    factory->createVectorSpace(comm, "reducer_test");
  fei::SharedPtr<fei::MatrixGraph> mgraph =
    factory->createMatrixGraph(vspace, vspace, "reducer_test");

  int fieldID = 0;
  int fieldSize = 1;
  vspace->defineFields(1, &fieldID, &fieldSize);
  int idType = 0;
  vspace->defineIDTypes(1, &idType);

  vspace->addDOFs(fieldID, idType, eqns.size(), &eqns[0]);

  mgraph->initComplete();

  fei::SharedPtr<fei::FillableMat> target(new fei::FillableMat);
  fei::Matrix_Impl<fei::FillableMat> feimat(target, mgraph, Kr.getNumRows());

  //now add the unreduced K** matrices to the reducer
  addFillableMatToReducer(Kii, reducer, feimat);
  addFillableMatToReducer(Kid, reducer, feimat);
  addFillableMatToReducer(Kdi, reducer, feimat);
  addFillableMatToReducer(Kdd, reducer, feimat);

  reducer.assembleReducedMatrix(feimat);

  fei::CSRMat csrtarget;
  csrtarget = *target;
  TEUCHOS_TEST_EQUALITY(Kr == csrtarget, true, out, success);
}

#endif

TEUCHOS_UNIT_TEST(Reducer, test_Reducer_test2)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return;

  fei::SharedPtr<fei::FillableMat> D(new fei::FillableMat);
  D->putCoef(0, 1, 0.5); D->putCoef(0, 2, 0.5);
  D->putCoef(3, 1, 0.5); D->putCoef(3, 2, 0.5);

  fei::SharedPtr<fei::CSVec> g;
  fei::Reducer reducer(D, g, comm);

  int num = 4;
  std::vector<int> eqns1(num), eqns2(num);
  for(int i=0; i<num; ++i) {
    eqns1[i] = i;
    eqns2[i] = i;
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateToReducedEqn(eqns1[i]);
    }
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateFromReducedEqn(eqns1[i]);
    }
  }

  TEUCHOS_TEST_EQUALITY(eqns1 == eqns2, true, out, success);
}

TEUCHOS_UNIT_TEST(Reducer, test_Reducer_test3)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return;

  fei::SharedPtr<fei::FillableMat> D(new fei::FillableMat);
  D->putCoef(1, 0, 0.5); D->putCoef(1, 2, 0.5);
  D->putCoef(3, 4, 0.5); D->putCoef(3, 2, 0.5);
  D->putCoef(5, 6, 0.5); D->putCoef(5, 4, 0.5);

  fei::SharedPtr<fei::CSVec> g;
  fei::Reducer reducer(D, g, comm);

  int num = 7;
  std::vector<int> eqns1(num), eqns2(num);
  for(int i=0; i<num; ++i) {
    eqns1[i] = i;
    eqns2[i] = i;
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateToReducedEqn(eqns1[i]);
    }
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateFromReducedEqn(eqns1[i]);
    }
  }

  TEUCHOS_TEST_EQUALITY(eqns1 == eqns2, true, out, success);
}

TEUCHOS_UNIT_TEST(Reducer, test_Reducer_test4)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return;

  fei::SharedPtr<fei::FillableMat> D(new fei::FillableMat);
  D->putCoef(3, 0, 0.5); D->putCoef(3, 2, 0.5);
  D->putCoef(4, 1, 0.5); D->putCoef(4, 2, 0.5);
  D->putCoef(5, 6, 0.5); D->putCoef(5, 7, 0.5);

  fei::SharedPtr<fei::CSVec> g;
  fei::Reducer reducer(D, g, comm);

  int num = 9;
  std::vector<int> eqns1(num), eqns2(num);
  for(int i=0; i<num; ++i) {
    eqns1[i] = i;
    eqns2[i] = i;
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateToReducedEqn(eqns1[i]);
    }
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateFromReducedEqn(eqns1[i]);
    }
  }

  TEUCHOS_TEST_EQUALITY(eqns1 == eqns2, true, out, success);
}

TEUCHOS_UNIT_TEST(Reducer, test_Reducer_test5)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return;

  fei::SharedPtr<fei::FillableMat> D(new fei::FillableMat);
//  D.putCoef(28, 0, 1.0);
//  D.putCoef(30, 1, 1.0);
//  D.putCoef(34, 2, 1.0);
//  D.putCoef(37, 3, 1.0);
//  D.putCoef(57, 4, 1.0);
//  D.putCoef(60, 5, 1.0);
//  D.putCoef(63, 6, 1.0);
//  D.putCoef(66, 7, 1.0);
//  D.putCoef(68, 8, 1.0);
//  D.putCoef(71, 9, 1.0);
//  D.putCoef(72, 10, 1.0);
//  D.putCoef(74, 11, 1.0);
//  D.putCoef(76, 12, 1.0);
//  D.putCoef(79, 13, 1.0);
//  D.putCoef(81, 14, 1.0);
//  D.putCoef(83, 15, 1.0);
  D->putCoef(1, 0, 1.0);
  D->putCoef(3, 2, 1.0);
  D->putCoef(4, 5, 1.0);
  D->putCoef(6, 7, 1.0);

  fei::SharedPtr<fei::CSVec> g;
  fei::Reducer reducer(D, g, comm);

  int num = 8;
  std::vector<int> eqns1(num), eqns2(num);
  for(int i=0; i<num; ++i) {
    eqns1[i] = i;
    eqns2[i] = i;
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateToReducedEqn(eqns1[i]);
    }
  }

  for(int i=0; i<num; ++i) {
    if (!reducer.isSlaveEqn(eqns2[i])) {
      eqns1[i] = reducer.translateFromReducedEqn(eqns1[i]);
    }
  }

  TEUCHOS_TEST_EQUALITY(eqns1 == eqns2, true, out, success);
}

}//namespace <anonymous>

