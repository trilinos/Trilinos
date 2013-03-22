/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/test_Utils.hpp>

#include <fei_ArrayUtils.hpp>
#include <fei_utils.hpp>
#include <fei_CommUtils.hpp>
#include <snl_fei_Utils.hpp>
#include <fei_Param.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_SharedPtr.hpp>
#include <cmath>

#undef fei_file
#define fei_file "test_Utils.cpp"
#include <fei_ErrMacros.hpp>

test_Utils::test_Utils(MPI_Comm comm)
  : tester(comm)
{
}

test_Utils::~test_Utils()
{
}

void test_Utils_binarySearch()
{
  std::vector<int> intarray;
  intarray.push_back(1);
  intarray.push_back(2);
  intarray.push_back(5);
  intarray.push_back(6);
  intarray.push_back(9);

  int offset = 0;
  int insertPoint = -1;

  FEI_COUT << "testing correctness of fei::binarySearch(int,int*,int,int)...";

  offset = fei::binarySearch(0, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != -1 || insertPoint != 0) {
    throw std::runtime_error("fei::binarySearch test failed 1.");
  }

  offset = fei::binarySearch(2, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != 1) {
    throw std::runtime_error("fei::binarySearch test failed 2.");
  }

  offset = fei::binarySearch(3, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != -1 || insertPoint != 2) {
    throw std::runtime_error("fei::binarySearch test failed 3.");
  }

  offset = fei::binarySearch(4, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != -1 || insertPoint != 2) {
    throw std::runtime_error("fei::binarySearch test failed 4.");
  }

  offset = fei::binarySearch(9, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != 4) {
    throw std::runtime_error("fei::binarySearch test failed 5.");
  }

  offset = fei::binarySearch(8, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != -1 || insertPoint != 4) {
    throw std::runtime_error("fei::binarySearch test failed 6.");
  }

  offset = fei::binarySearch(10, &intarray[0], intarray.size(),
				 insertPoint);
  if (offset != -1 || insertPoint != 5) {
    throw std::runtime_error("fei::binarySearch test failed 7.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

int test_Utils::runtests()
{
  if (numProcs_ < 2) {
    test_Utils_binarySearch();

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

int test_Utils::serialtest1()
{
  FEI_COUT << "testing snl_fei::leading_substring_length...";

  static char string1[] = "test ";
  string1[4] = '\0';
  if (snl_fei::leading_substring_length(string1) != 4) {
    ERReturn(-1);
  }

  static char string2[] = "second test";
  if (snl_fei::leading_substring_length(string2) != 6) {
    ERReturn(-1);
  }

  static char string3[] = "third test";
  string3[5] = '\t';
  if (snl_fei::leading_substring_length(string3) != 5) {
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}

int test_Utils::serialtest2()
{
  FEI_COUT << "testing snl_fei::getDoubleParamValue...";

  static char string1[] = "DOUBLE1 1.0";
  static char string2[] = "DOUBLE2 1.0e+0";
  static char string3[] = "DOUBLE3 1.0E+0";
  static char string4[] = "DOUBLE4 1";

  std::vector<char*> params;
  params.push_back(string1);
  params.push_back(string2);
  params.push_back(string3);
  params.push_back(string4);

  double d1,d2,d3,d4;

  CHK_ERR( snl_fei::getDoubleParamValue("DOUBLE1",
					params.size(), &params[0],d1));
  CHK_ERR( snl_fei::getDoubleParamValue("DOUBLE2",
					params.size(), &params[0],d2));
  CHK_ERR( snl_fei::getDoubleParamValue("DOUBLE3",
					params.size(), &params[0],d3));
  CHK_ERR( snl_fei::getDoubleParamValue("DOUBLE4",
					params.size(), &params[0],d4));

  if (std::abs(d1 - 1.0) > 1.e-49 || std::abs(d2 - 1.0) > 1.e-49 ||
      std::abs(d3 - 1.0) > 1.e-49 || std::abs(d4 - 1.0) > 1.e-49) {
    ERReturn(-1);
  }

  FEI_COUT <<"ok"<<FEI_ENDL;

  return(0);
}

int test_Utils::serialtest3()
{
  FEI_COUT << "testing fei::Param and fei::ParameterSet...";

  fei::Param param1("string-param", "garbage value");
  fei::Param param2("double-param", 2.5);
  fei::Param param3("int-param", 1);

  if (param1.getType() != fei::Param::STRING) {
    ERReturn(-1);
  }

  if (param2.getType() != fei::Param::DOUBLE) {
    ERReturn(-1);
  }

  if (param3.getType() != fei::Param::INT) {
    ERReturn(-1);
  }

  fei::ParameterSet paramset;
  paramset.add(fei::Param("string-param", "garbage value"));
  paramset.add(param2);
  paramset.add(param3);

  if (paramset.size() != 3) {
    ERReturn(-1);
  }

  fei::ParameterSet::const_iterator
    iter = paramset.begin(),
    iter_end = paramset.end();

  int i=0;
  for(; iter != iter_end; ++iter) {
    if (i==3) {
      ERReturn(-1);
    }
    ++i;
  }
 
  if (paramset.get("int-param") == NULL) {
    ERReturn(-1);
  }

  int dummy;
  int err = paramset.getIntParamValue("int-param", dummy);
  if (err != 0) {
    ERReturn(-1);
  }

  if (dummy != 1) {
    ERReturn(-1);
  }

  std::string dummychars;
  err = paramset.getStringParamValue("string-param", dummychars);
  if (err != 0) {
    ERReturn(-1);
  }

  if ("garbage value" != dummychars) {
    ERReturn(-1);
  }

  //if (!snl_fei::leadingSubstring("garbage-value", "garbage")) {
  //  ERReturn(-1);
  //}

  //if (snl_fei::leadingSubstring("garb-value", "garbage")) {
  //  ERReturn(-1);
  //}

  std::vector<std::string> stdstrings;
  std::string tempstr;

  tempstr = "string-param garbage value";
  stdstrings.push_back(tempstr);

  tempstr = "int-param 58";
  stdstrings.push_back(tempstr);

  tempstr = "real-param 45.e-2";
  stdstrings.push_back(tempstr);

  fei::ParameterSet pset;
  fei::utils::parse_strings(stdstrings, " ", pset);

  err = pset.getStringParamValue("string-param", dummychars);
  if ("garbage value" != dummychars) {
    ERReturn(-1);
  }

  err = pset.getIntParamValue("int-param", dummy);
  if (dummy != 58) {
    ERReturn(-1);
  }

  double ddummy;
  err = pset.getDoubleParamValue("real-param", ddummy);
  if (std::abs(ddummy - 45.e-2) > 1.e-49) {
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}

void test_Utils_function_that_throws()
{
  throw std::runtime_error("testing...");
}

int test_Utils::test1()
{
  FEI_COUT << "testing std::runtime_error...";

  bool exc_thrown_and_caught = false;

  try {
    test_Utils_function_that_throws();
  }
  catch(std::runtime_error& exc) {
    std::string str(exc.what());
    if (str == "testing...") {
      exc_thrown_and_caught = true;
    }
  }

  if (!exc_thrown_and_caught) {
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

bool test_Utils_dummy_destroyed = true;

class test_Utils_dummy {
public:
  test_Utils_dummy() {test_Utils_dummy_destroyed = false;}
  ~test_Utils_dummy()
  {
    test_Utils_dummy_destroyed = true;
  }
};

int test_Utils_test_SharedPtr()
{
  //In this function, make sure the global bool is set to true, then create
  //the fei::SharedPtr and make sure that the global bool has been set to false.
  //If so, return 0, otherwise return -1.
  //When we return, the SharedPtr goes out of scope which should destroy the
  //test-dummy and cause the global bool to get set back to true. The code 
  //that's calling this function will verify that.

  test_Utils_dummy_destroyed = true;
  fei::SharedPtr<test_Utils_dummy> ptr(new test_Utils_dummy);
  if (test_Utils_dummy_destroyed == true) return(-1);
  else return(0);
}

int test_Utils::test2()
{
  FEI_COUT << "testing fei::SharedPtr...";
  int err = test_Utils_test_SharedPtr();
  if (err != 0) {
    ERReturn(-1);
  }

  if (test_Utils_dummy_destroyed != true) {
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;
 return(0);
}

int test_Utils::test3()
{
  FEI_COUT << "testing snl_fei::copy2DToColumnContig...";

  int numrows1 = 3;
  int numcols1 = 4;
  int numrows2 = 4;
  int numcols2 = 3;

  int i, j;
  int len1 = numrows1*numcols1;
  int len2 = numrows2*numcols2;

  double** table2d_1 = new double*[numrows1];
  for(i=0; i<numrows1; ++i) {
    table2d_1[i] = new double[numcols1];
    for(j=0; j<numcols1; ++j) {
      table2d_1[i][j] = j*numrows1+i;
    }
  }

  double** table2d_2 = new double*[numcols2];
  for(j=0; j<numcols2; ++j) {
    table2d_2[j] = new double[numrows2];
    for(i=0; i<numrows2; ++i) {
      table2d_2[j][i] = j*numrows2+i;
    }
  }

  double* cc1 = new double[len1];
  double* cc2 = new double[len2];

  snl_fei::copy2DToColumnContig(numrows1, numcols1, table2d_1,
				FEI_DENSE_ROW, cc1);

  snl_fei::copy2DToColumnContig(numrows2, numcols2, table2d_2,
				FEI_DENSE_COL, cc2);

  for(i=0; i<len1; ++i) {
    if (std::abs(cc1[i] - cc2[i]) > 1.e-49) {
      throw std::runtime_error("column-contig arrays not equal.");
    }
  }

  for(j=0; j<numrows1; ++j) delete [] table2d_1[j];
  delete [] table2d_1;
  delete [] cc1;
  delete [] cc2;

  FEI_COUT << "ok"<<FEI_ENDL;

  FEI_COUT << "testing snl_fei::copy2DBlockDiagToColumnContig...";

  numrows1 = 12;
  int numBlocks = 3;
  int* blockSizes = new int[numBlocks];
  for(i=0; i<numBlocks; ++i) {
    blockSizes[i] = 4;
  }

  table2d_1 = new double*[numrows1];
  for(i=0; i<numrows1; ++i) {
    table2d_1[i] = new double[4];
    for(j=0; j<4; ++j) {
      table2d_1[i][j] = 1.0*i*4+j;
    }
  }

  len1 = numrows1*4;
  cc1 = new double[len1];

  snl_fei::copy2DBlockDiagToColumnContig(numBlocks, blockSizes, table2d_1,
					 FEI_BLOCK_DIAGONAL_ROW, cc1);

  for(i=0; i<len1; ++i) {
    if (std::abs(1.0*i - cc1[i]) > 1.e-49) {
      throw std::runtime_error("copy2DBlockDiagToColumnContig row test failed.");
    }
  }

  for(j=0; j<numrows1; ++j) delete [] table2d_1[j];
  delete [] table2d_1;
  for(j=0; j<numcols2; ++j) delete [] table2d_2[j];
  delete [] table2d_2;

  delete [] cc1;
  delete [] blockSizes;

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_Utils::test4()
{
  return(0);
}
