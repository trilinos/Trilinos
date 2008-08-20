
#include <fei_iostream.hpp>
#include <fei_DirichletBCRecord.hpp>
#include <fei_DirichletBCManager.hpp>

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

void test_DirBCManager()
{
  FEI_COUT << "testing fei::DirichletBCManager...";

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

bool test_DirBC::run(MPI_Comm /*comm*/)
{
  test_DirBCRecord();

  test_DirBCManager();

  return true;
}

