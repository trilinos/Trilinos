/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/fei_test_utils.hpp>

#include <test_utils/test_misc.hpp>

#include <test_utils/test_Factory_helper.hpp>

#include <fei_FieldMask.hpp>
#include <snl_fei_RecordCollection.hpp>

#include <fei_VectorSpace.hpp>
#include <fei_Vector_Impl.hpp>

#undef fei_file
#define fei_file "test_misc.cpp"
#include <fei_ErrMacros.hpp>

test_misc::test_misc(MPI_Comm comm)
  : tester(comm)
{
}

test_misc::~test_misc()
{
}

void test_misc_FieldMask()
{
  FEI_COUT << "testing fei::FieldMask...";

  //A general test of fei::FieldMask.

  unsigned numFields = 5;
  std::vector<int> fieldIDs(numFields);
  std::vector<int> fieldSizes(numFields);
  int checkNumIndices = 0;
  for(unsigned i=0; i<numFields; ++i) {
    fieldIDs[i] = i;
    fieldSizes[i] = i;
    checkNumIndices += i;
  }

  fei::FieldMask fieldMask;

  for(int i=fieldIDs.size()-1; i>= 0; --i) {
    fieldMask.addField(fieldIDs[i], fieldSizes[i]);
  }

  std::vector<int>& maskFields = fieldMask.getFieldIDs();
  std::vector<int>& maskFieldSizes = fieldMask.getFieldSizes();

  if (maskFields != fieldIDs) {
    throw std::runtime_error("FieldMask test failed.");
  }

  if (maskFieldSizes != fieldSizes) {
    throw std::runtime_error("FieldMask size test failed.");
  }

  int checkOffset = 0;
  for(unsigned j=0; j<fieldIDs.size(); ++j) {
    int offset = -1;
    fieldMask.getFieldEqnOffset(fieldIDs[j], offset);
    if (offset != checkOffset) {
      throw std::runtime_error("FieldMask offset test failed.");
    }
    checkOffset += j;
  }

  int numIndices = fieldMask.getNumIndices();
  if (numIndices != checkNumIndices) {
    throw std::runtime_error("FieldMask numIndices test failed.");
  }

  bool exc_caught = false;
  try {
    fieldMask.addField(-1, 0);
  }
  catch (...) {
    exc_caught = true;
  }

  if (!exc_caught) {
    throw std::runtime_error("FieldMask failed to throw on negative fieldID.");
  }

  fieldMask.addField(2, 2);

  if (fieldMask.getNumFields() != numFields) {
    throw std::runtime_error("FieldMask getNumFields test failed.");
  }

  int fieldid1 = 0;
  int fieldid2 = 1;
  int fieldid3 = 2;
  int fieldsize = 1;

  fei::FieldMask fm1(1, &fieldid1, &fieldsize);
  fei::FieldMask fm2(1, &fieldid2, &fieldsize);
  fei::FieldMask fm3(1, &fieldid3, &fieldsize);

  fei::FieldMask fm12(fm1);
  fm12.addField(fieldid2, fieldsize);

  fei::FieldMask fm123(fm1);
  fm123.addField(fieldid2, fieldsize);
  fm123.addField(fieldid3, fieldsize);

  if (fm1.getMaskID() == fm2.getMaskID()) {
    throw std::runtime_error("FieldMask getMaskID test failed.");
  }

  if (fm2.getMaskID() == fm12.getMaskID()) {
    throw std::runtime_error("FieldMask getMaskID2 test failed.");
  }

  if (fm12.getMaskID() !=
      fei::FieldMask::calculateMaskID(fm1, fieldid2)){
    throw std::runtime_error("FieldMask getMaskID3 test failed.");
  }

  if (fm12.getMaskID() == fm3.getMaskID()) {
    throw std::runtime_error("FieldMask getMaskID4 test failed.");
  }

  if (fm123.getMaskID() != 
      fei::FieldMask::calculateMaskID(fm12, fieldid3)){
    throw std::runtime_error("FieldMask getMaskID5 test failed.");
  }

  if (fm3.getMaskID() == fm123.getMaskID()) {
    throw std::runtime_error("FieldMask getMaskID6 test failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_misc_RecordCollection()
{
  FEI_COUT << "testing snl_fei::RecordCollection...";

  int fieldID0 = 0;
  int fieldID1 = 1;
  int fieldID2 = 2;
  int fieldSize = 1;
  int ID0 = 0;
  int ID1 = 1;

  std::vector<fei::FieldMask*> fieldMasks;

  snl_fei::RecordCollection recColl(0);

  int* records = new int[1];

  recColl.initRecords(fieldID0, fieldSize, 1, &ID0,
		      fieldMasks, records);

  recColl.initRecords(fieldID1, fieldSize, 1, &ID0,
		      fieldMasks, records);

  recColl.initRecords(fieldID2, fieldSize, 1, &ID0,
		      fieldMasks, records);

  recColl.initRecords(fieldID0, fieldSize, 1, &ID1,
		      fieldMasks, records);

  recColl.initRecords(fieldID1, fieldSize, 1, &ID1,
		      fieldMasks, records);

  recColl.initRecords(fieldID2, fieldSize, 1, &ID1,
		      fieldMasks, records);

  if (fieldMasks.size() != 5) {
    throw std::runtime_error("RecordCollection fieldMasks.length test failed.");
  }

  std::vector<fei::Record<int> >& rvec = recColl.getRecords();

  std::vector<fei::Record<int> >::iterator
    r_iter = rvec.begin(),
    r_end = rvec.end();

  int numIndices = 0;
  for(; r_iter != r_end; ++r_iter) {
    numIndices += (*r_iter).getFieldMask()->getNumIndices();
  }

  if (numIndices != 6) {
    throw std::runtime_error("RecordCollection numIndices test failed.");
  }

  delete [] records;
  for(unsigned i=0; i<fieldMasks.size(); ++i) delete fieldMasks[i];

  FEI_COUT << "ok"<<FEI_ENDL;
}

int test_misc::runtests()
{
  if (numProcs_ < 2) {
    test_misc_FieldMask();
    test_misc_RecordCollection();

    CHK_ERR( serialtest1() );
    CHK_ERR( serialtest2() );
    CHK_ERR( serialtest3() );
  }

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_misc::serialtest1()
{
  FEI_COUT << "testing fei_test_utils::within_percentage_margin...";
  double value1 = 65000.0;
  double value2 = 6500.0;
  bool result = fei_test_utils::within_percentage_margin(value1, value2, 10);
  if (result == true) {
    ERReturn(-1);
  }

  value1 = 6500.0;
  value2 = 6500.1;
  result = fei_test_utils::within_percentage_margin(value1,value2,1);
  if (result != true) {
    ERReturn(-1);
  }

  value1 = -10.0;
  value2 = 0.0;
  result = fei_test_utils::within_percentage_margin(value1,value2,30);
  if (result == true) {
    ERReturn(-1);
  }

  value1 = -1.e-18;
  value2 = 1.e-15;
  result = fei_test_utils::within_percentage_margin(value1,value2,10);
  if (result != true) {
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_misc::serialtest2()
{
  FEI_COUT << "testing fei::lowerBound...";
  std::vector<int> list(5);

  list[0] = 1;
  list[1] = 4;
  list[2] = 6;
  list[3] = 7;
  list[4] = 11;

  int item = 0;
  int lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 0) {
    throw std::runtime_error("failed test 1");
  }

  item = 1;
  lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 0) {
    throw std::runtime_error("failed test 2");
  }

  item = 2;
  lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 1) {
    throw std::runtime_error("failed test 3");
  }

  item = 7;
  lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 3) {
    throw std::runtime_error("failed test 4");
  }

  item = 9;
  lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 4) {
    throw std::runtime_error("failed test 5");
  }

  item = 11;
  lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 4) {
    throw std::runtime_error("failed test6");
  }

  item = 12;
  lowerbound = fei::lowerBound<int>(item, &list[0], list.size());

  if (lowerbound != 5) {
    throw std::runtime_error("failed test 7");
  }

  lowerbound = fei::lowerBound<int>(item, (int*)0, (int)0);

  if (lowerbound != 0) {
    throw std::runtime_error("failed test 8");
  }

  std::vector<int> list2;
  list2.push_back(2);

  item = 2;
  lowerbound = fei::lowerBound<int>(item, &list2[0], list2.size());

  if (lowerbound != 0) {
    throw std::runtime_error("failed test 9");
  }

  item = 5;
  lowerbound = fei::lowerBound<int>(item, &list2[0], list2.size());

  if (lowerbound != 1) {
    throw std::runtime_error("failed test 10");
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}

int test_misc::serialtest3()
{
  FEI_COUT << "testing snl_fei::MapContig<fei::ctg_set<int>*>...";

  snl_fei::MapContig<fei::ctg_set<int>*> mc(1,5);
  fei_Pool_alloc<fei::ctg_set<int> > pool_alloc;

  static fei::ctg_set<int> dummy;

  for(int i=1; i<6; ++i) {
    fei::ctg_set<int>* newset = pool_alloc.allocate(1);
    pool_alloc.construct(newset,dummy);
  

    for(int j=0; j<3; ++j) {
      newset->insert2(j);
    }

    std::pair<int,fei::ctg_set<int>*> newpair(i,newset);
    mc.insert(newpair);
  }
  
  snl_fei::MapContig<fei::ctg_set<int>*> m_copy(mc);

  if (m_copy.size() != mc.size()) {
    throw std::runtime_error("failed test 1.");
  }

  snl_fei::MapContig<fei::ctg_set<int>*>::iterator
    mc_iter = mc.begin(),
    mc_end = mc.end();

  snl_fei::MapContig<fei::ctg_set<int>*>::iterator
    c_iter = m_copy.begin(),
    c_end = m_copy.end();

  for(; mc_iter != mc_end; ++mc_iter) {
    std::pair<int,fei::ctg_set<int>*> mc_pair = *mc_iter;
    std::pair<int,fei::ctg_set<int>*> c_pair = *c_iter;

    if (mc_pair.first != c_pair.first) {
      throw std::runtime_error("failed test 2.");
    }

    if (*(mc_pair.second) != *(c_pair.second)) {
      throw std::runtime_error("failed test 3.");
    }
    ++c_iter;
  }

  mc_iter = mc.begin();
  for(; mc_iter != mc_end; ++mc_iter) {
    pool_alloc.destroy((*mc_iter).second);
  }

  FEI_COUT << "ok" << FEI_ENDL;

  return(0);
}

int test_misc::test1()
{

  return(0);
}

int test_misc::test2()
{

 return(0);
}

int test_misc::test3()
{

  return(0);
}

int test_misc::test4()
{
  return(0);
}
