
#include <fei_iostream.hpp>
#include <fei_mpi.h>
#include <test_utils/LibraryFactory.hpp>
#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>
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

#include <fei_unit_Reducer.hpp>

#include <vector>
#include <cmath>

int test_reducer_unit1()
{
  FEI_COUT << "testing fei::Reducer: basic test 1...";

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

  if (!reducer.isSlaveEqn(1)) {
    ERReturn(-1);
  }

  bool exception_caught = false;
  //calling translateToReducedEqn with a slave eqn should cause reducer
  //to throw an exception.
  try {
    reducer.translateToReducedEqn(1);
  }
  catch(std::runtime_error& exc) {
    exception_caught = true;
  }

  if (!exception_caught) {
    FEI_COUT << " ERROR, expected exception not thrown."<<FEI_ENDL;
    ERReturn(-2);
  }

  int reducedEqn = reducer.translateToReducedEqn(3);
  if (reducedEqn != 2) {
    ERReturn(-3);
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
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

void addFillableMatToReducer(fei::FillableMat& mat, fei::Reducer& reducer)
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
      reducer.addMatrixValues(1, &row, 1, &col, &coefPtr, true);
    }
  }
}

int test_Reducer_test1(MPI_Comm comm)
{
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs > 1) return(0);

  FEI_COUT << "testing fei::Reducer matrix assembly...";

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

  //now add the unreduced K** matrices to the reducer
  addFillableMatToReducer(Kii, reducer);
  addFillableMatToReducer(Kid, reducer);
  addFillableMatToReducer(Kdi, reducer);
  addFillableMatToReducer(Kdd, reducer);

  fei::SharedPtr<fei::Factory> factory;
  try {
    factory = fei::create_fei_Factory(comm, "Trilinos");
  }
  catch(std::runtime_error& exc) {
    FEI_COUT << "couldn't create Trilinos factory."<<FEI_ENDL;
    return(0);
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

  vspace->addDOFs(fieldID, 1, idType, eqns.size(), &eqns[0]);

  mgraph->initComplete();

  fei::SharedPtr<fei::FillableMat> target(new fei::FillableMat);
  fei::Matrix_Impl<fei::FillableMat> feimat(target, mgraph, Kr.getNumRows());

  reducer.assembleReducedMatrix(feimat);

  fei::CSRMat csrtarget;
  csrtarget = *target;
  if (Kr != csrtarget) {
    FEI_COUT << "reducer produced incorrect matrix."<<FEI_ENDL;
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_Reducer_test2(MPI_Comm comm)
{
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return(0);

  FEI_COUT << "testing fei::Reducer::translate[To/From]ReducedEqn case 1...";

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

  if (eqns1 != eqns2) {
    FEI_COUT << "failed."<<FEI_ENDL;
    ERReturn(-1);
  }

  FEI_COUT << "ok" <<FEI_ENDL;
  return(0);
}

int test_Reducer_test3(MPI_Comm comm)
{
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return(0);

  FEI_COUT << "testing fei::Reducer::translate[To/From]ReducedEqn case 2...";

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

  if (eqns1 != eqns2) {
    FEI_COUT << "failed."<<FEI_ENDL;
    ERReturn(-1);
  }

  FEI_COUT << "ok" <<FEI_ENDL;
  return(0);
}

int test_Reducer_test4(MPI_Comm comm)
{
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return(0);

  FEI_COUT << "testing fei::Reducer::translate[To/From]ReducedEqn case 3...";

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

  if (eqns1 != eqns2) {
    FEI_COUT << "failed."<<FEI_ENDL;
    ERReturn(-1);
  }

  FEI_COUT << "ok" <<FEI_ENDL;
  return(0);
}

int test_Reducer_test5(MPI_Comm comm)
{
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs>1) return(0);

  FEI_COUT << "testing fei::Reducer::translate[To/From]ReducedEqn case 4...";

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

  if (eqns1 != eqns2) {
    FEI_COUT << "failed."<<FEI_ENDL;
    ERReturn(-1);
  }

  FEI_COUT << "ok" <<FEI_ENDL;
  return(0);
}

bool test_Reducer::run(MPI_Comm comm)
{
  test_reducer_unit1();

  if (test_Reducer_test1(comm) != 0) {
    throw std::runtime_error("test_Reducer_test1 failed.");
  }

  if (test_Reducer_test2(comm) != 0) {
    throw std::runtime_error("test_Reducer_test2 failed.");
  }

  if (test_Reducer_test3(comm) != 0) {
    throw std::runtime_error("test_Reducer_test3 failed.");
  }

  if (test_Reducer_test4(comm) != 0) {
    throw std::runtime_error("test_Reducer_test4 failed.");
  }

  if (test_Reducer_test5(comm) != 0) {
    throw std::runtime_error("test_Reducer_test5 failed.");
  }

  return true;
}

