
#include <fei_iostream.hpp>
#include <fei_DirichletBCRecord.hpp>
#include <fei_DirichletBCManager.hpp>

#include <fei_MatrixGraph_Impl2.hpp>
#include <fei_Matrix_Impl.hpp>

#include <fei_unit_DirBC.hpp>

#include <vector>
#include <cmath>

void test_DirBCRecord()
{
  FEI_COUT << "testing fei::DirichletBCRecord...";

  fei::DirichletBCRecord dbc1;
  fei::DirichletBCRecord dbc2;

  dbc1.IDType = 0;
  dbc1.ID = 19;
  dbc1.fieldID = 2;
  dbc1.whichComponent = 1;

  dbc2.IDType = 0;
  dbc2.ID = 20;
  dbc2.fieldID = 2;
  dbc2.whichComponent = 1;

  if (!(dbc1 != dbc2)) {
    throw std::runtime_error("DirBCRecord::operator!= test 1 failed.");
  }

  fei::less_DirichletBCRecord lessdbc;

  if (!lessdbc(dbc1,dbc2)) {
    throw std::runtime_error("DirBCRecord less test 1 failed.");
  }

  dbc2.ID = 19;
  dbc2.whichComponent = 2;

  if (!(dbc1 != dbc2)) {
    throw std::runtime_error("DirBCRecord::operator!= test 2 failed.");
  }

  if (!lessdbc(dbc1,dbc2)) {
    throw std::runtime_error("DirBCRecord less test 2 failed.");
  }

  dbc2.whichComponent = 1;

  if (dbc1 != dbc2) {
    throw std::runtime_error("DirBCRecord::operator!= test 3 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_DirBCManager_addBCRecords()
{
  FEI_COUT << "testing fei::DirichletBCManager::addBCRecords...";

  fei::DirichletBCManager bcmgr;

  int idtype = 0;
  int fieldID = 0;
  int offsetIntoField = 1;
  std::vector<int> ids(5);
  std::vector<double> vals(5);

  ids[0] = 2; vals[0] = 2.0;
  ids[1] = 1; vals[1] = 1.0;
  ids[2] = 3; vals[2] = 3.0;
  ids[3] = 4; vals[3] = 4.0;
  ids[4] = 0; vals[4] = 0.0;

  bcmgr.addBCRecords(5, idtype, fieldID, offsetIntoField,
                     &ids[0], &vals[0]);
 
  if (bcmgr.getNumBCRecords() != 5) {
    throw std::runtime_error("test_DirBCManager test 1 failed.");
  }

  std::vector<int> offsetsIntoField(5, 1);

  bcmgr.addBCRecords(5, idtype, fieldID, &ids[0],
                     &offsetsIntoField[0], &vals[0]);
 
  if (bcmgr.getNumBCRecords() != 5) {
    throw std::runtime_error("test_DirBCManager test 2 failed.");
  }

  offsetsIntoField[1] = 0;
  offsetsIntoField[3] = 0;
  offsetsIntoField[4] = 0;

  bcmgr.addBCRecords(5, idtype, fieldID, &ids[0],
                     &offsetsIntoField[0], &vals[0]);
 
  if (bcmgr.getNumBCRecords() != 8) {
    throw std::runtime_error("test_DirBCManager test 3 failed.");
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

void test_DirBCManager_finalizeBCEqns(MPI_Comm comm)
{
  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs > 1) {
    FEI_COUT << "skipping test of fei::DirichletBCManager::finalizeBCEqn, which only"
     << " runs on 1 proc." << FEI_ENDL;
    return;
  }

  FEI_COUT << "testing fei::DirichletBCManager::finalizeBCEqn...";

  fei::DirichletBCManager bcmgr;

  int idtype = 0;
  int fieldID = 0;
  int fieldSize = 2;
  int offsetIntoField = 1;
  std::vector<int> ids(5);
  std::vector<double> vals(5);

  ids[0] = 2; vals[0] = 2.0;
  ids[1] = 1; vals[1] = 1.0;
  ids[2] = 3; vals[2] = 3.0;
  ids[3] = 4; vals[3] = 4.0;
  ids[4] = 0; vals[4] = 0.0;

  bcmgr.addBCRecords(5, idtype, fieldID, offsetIntoField,
                     &ids[0], &vals[0]);
 
  fei::SharedPtr<fei::VectorSpace> vspace(new fei::VectorSpace(comm));

  vspace->defineFields(1, &fieldID, &fieldSize);
  vspace->defineIDTypes(1, &idtype);

  fei::SharedPtr<fei::MatrixGraph> mgraph(new fei::MatrixGraph_Impl2(vspace, vspace));

  int patternID = 0;
  int numIDs = 1;

  mgraph->definePattern(patternID, numIDs, idtype);

  int blockID = 0;

  mgraph->initConnectivityBlock(blockID, ids.size(), patternID);

  for(size_t i = 0; i<ids.size(); ++i) {
    mgraph->initConnectivity(blockID, i, &ids[i]);
  }

  mgraph->initComplete();

  fei::SharedPtr<fei::FillableMat> inner(new fei::FillableMat);
  fei::SharedPtr<fei::Matrix_Impl<fei::FillableMat> > feimat(new fei::Matrix_Impl<fei::FillableMat>(inner, mgraph, ids.size()));

  int err = bcmgr.finalizeBCEqns(*feimat, false);
  if (err != 0) {
    throw std::runtime_error("test_DirBCManager test 4 failed.");
  }

  if (feimat->getGlobalNumRows() != feimat->getLocalNumRows() ||
      feimat->getGlobalNumRows() != (int)ids.size()) {
    std::ostringstream osstr;
    osstr << "test_DirBCManager test 5 failed, num-rows == " << feimat->getGlobalNumRows()
         << ", expected " << ids.size();
    throw std::runtime_error(osstr.str());
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

bool test_DirBC::run(MPI_Comm comm)
{
  test_DirBCRecord();

  test_DirBCManager_addBCRecords();

  test_DirBCManager_finalizeBCEqns(comm);

  return true;
}

