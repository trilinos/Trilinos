/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_utils.hpp>
#include <test_utils/fei_test_utils.hpp>
#include <test_utils/test_benchmarks.hpp>
#include <snl_fei_Utils.hpp>
#include <fei_ctg_set.hpp>
#include <snl_fei_RaggedTable.hpp>
#include <snl_fei_RaggedTable_specialize.hpp>
#include <test_utils/HexBeam.hpp>

#undef fei_file
#define fei_file "test_benchmarks.cpp"

#include <fei_ErrMacros.hpp>

test_benchmarks::test_benchmarks(MPI_Comm comm)
  : tester(comm)
{
}

test_benchmarks::~test_benchmarks()
{
}

template<typename MAP_TYPE, typename SET_TYPE>
double time_raggedtable_insert(int len)
{
  double start_time = fei::utils::cpu_time();

  HexBeam hexcube(10, len, 1, HexBeam::OneD, 1, 0);

  snl_fei::RaggedTable<MAP_TYPE,SET_TYPE> table(0, hexcube.numLocalNodes());

  int numIndices = hexcube.numNodesPerElem();

  int* indices = new int[numIndices];

  int first = hexcube.firstLocalElem();

  for(int n=0; n<hexcube.numLocalElems(); ++n) {
    int elem = first+n;

    hexcube.getElemConnectivity(elem, indices);

    table.addIndices(numIndices, indices, numIndices, indices);
  }

  delete [] indices;

  double elapsed_time = fei::utils::cpu_time() - start_time;
  return(elapsed_time);
}

template<typename MAP_TYPE, typename SET_TYPE>
double benchmark_raggedtable()
{
  int len = 80;

  //first find a len such that time_taken is at least 1 second
  double time_taken = time_raggedtable_insert<MAP_TYPE,SET_TYPE>(len);
  while(time_taken < 1.0) {
    len *= 2;
    time_taken = time_raggedtable_insert<MAP_TYPE,SET_TYPE>(len);
  }

  //now repeat until time_taken passes 5 seconds
  time_taken = 0.0;
  int i=0;
  while(time_taken<5.0) {
    time_taken += time_raggedtable_insert<MAP_TYPE,SET_TYPE>(len);
    ++i;
  }

  return((double)(i*len)/time_taken);
}

void print_benchmark_banner()
{

  FEI_COUT.width(38);
  FEI_COUT << " Benchmark name         ";
  FEI_COUT.width(10);
  FEI_COUT << "Value  ";
  FEI_COUT.width(12);
  FEI_COUT << "gold-copy";
  FEI_COUT.width(10);
  FEI_COUT <<" Result "<<FEI_ENDL;

  FEI_COUT.width(38);
  FEI_COUT << " -----------------------";
  FEI_COUT.width(10);
  FEI_COUT << "-----  ";
  FEI_COUT.width(12);
  FEI_COUT << "---------";
  FEI_COUT.width(10);
  FEI_COUT <<" ------ "<<FEI_ENDL;
}

void print_benchmark_line(const char* name,
			  double value,
			  double goldvalue,
			  const char* passfail)
{
  FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
  FEI_COUT.precision(1);

  FEI_COUT.width(38);
  FEI_COUT << name;
  FEI_COUT.width(10);
  FEI_COUT << value;
  FEI_COUT.width(12);
  if (goldvalue < 0.0) FEI_COUT << "n/a";
  else FEI_COUT << goldvalue;
  FEI_COUT.width(10);
  FEI_COUT << passfail << FEI_ENDL;
}

std::string add_macro_values(const char* name)
{
  FEI_OSTRINGSTREAM osstr;
  osstr << name;

#if defined(FEI_PLATFORM) && defined(FEI_OPT_LEVEL)
  osstr << "_" << FEI_PLATFORM << "_" << FEI_OPT_LEVEL;
#else
  osstr << "_unknown_unknown";
#endif

  return(osstr.str());
}

int test_benchmarks::runtests()
{
  if (numProcs_ > 1) return(0);

  //CHK_ERR( test3() );

  FEI_COUT << FEI_ENDL
	   << "  ***** Benchmarks pass if within 10% of 'gold-copy' *****"
	   <<FEI_ENDL<<FEI_ENDL;

#if defined(FEI_PLATFORM) && defined(FEI_OPT_LEVEL)
  FEI_COUT << "  FEI_PLATFORM: "<<FEI_PLATFORM
	   <<", FEI_OPT_LEVEL: "<<FEI_OPT_LEVEL
	   <<FEI_ENDL<<FEI_ENDL;
  FEI_COUT << "  'gold-copy' benchmark values will be searched for in ./fei_utest_timings.txt"<<FEI_ENDL;
  FEI_COUT <<FEI_ENDL;
#else
  FEI_COUT << "  preprocessor macros FEI_PLATFORM and FEI_OPT_LEVEL aren't defined, so"<<FEI_ENDL;
  FEI_COUT << "  ./fei_utest_timings.txt will not be searched for 'gold-copy' benchmark values"<<FEI_ENDL<<FEI_ENDL;
#endif

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test4() );
  CHK_ERR( test5() );
  CHK_ERR( test6() );
  CHK_ERR( test7() );
  CHK_ERR( test8() );

  return(0);
}

int test_benchmarks::test1()
{
  FEI_COUT << "Following group of benchmarks inserts integers into ragged tables"
	   << " (simulating"<<FEI_ENDL
	   << "matrix-graph construction) using various data structures."<<FEI_ENDL
	   << "A higher number is better, indicating more insertions"
	   << " in fixed amount of time." << FEI_ENDL<<FEI_ENDL;

  print_benchmark_banner();

  int returnValue = 0;
  double value, goldvalue;
  std::string passfail;
  std::string testname;

  value = benchmark_raggedtable<std::map<int,std::set<int>*>,std::set<int> >();
  goldvalue = -1.0;
  passfail = " ";
  print_benchmark_line("std::map<std::set>", value, goldvalue, passfail.c_str());

  testname = add_macro_values("std::map<fei::ctg_set>");
  value = benchmark_raggedtable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> >();
  try {
    goldvalue = fei_test_utils::get_file_benchmark("./fei_utest_timings.txt",
					    testname.c_str());
    passfail = fei_test_utils::check_test_result(value, goldvalue, 10);
    if (passfail != "passed") returnValue = -1;
  }
  catch(...) {
    goldvalue = -1.0;
    passfail = " ";
  }

  print_benchmark_line("std::map<fei::ctg_set>", value, goldvalue, passfail.c_str());



  testname = add_macro_values("snl_fei::MapContig<fei::ctg_set>");
  value = benchmark_raggedtable<snl_fei::MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >();
  try {
    goldvalue = fei_test_utils::get_file_benchmark("./fei_utest_timings.txt",
					    testname.c_str());
    passfail = fei_test_utils::check_test_result(value, goldvalue, 10);
    if (passfail != "passed") returnValue = -1;
  }
  catch(...) {
    goldvalue = -1.0;
    passfail = " ";
  }

  print_benchmark_line("snl_fei::MapContig<fei::ctg_set>", value, goldvalue, passfail.c_str());



#ifdef FEI_HASH_MAP
  value = benchmark_raggedtable<FEI_HASH_MAP<int,FEI_HASH_SET<int>*>,FEI_HASH_SET<int> >();
  goldvalue = -1.0;
  passfail = " ";
  print_benchmark_line("hash_map<hash_set>", value, goldvalue, passfail.c_str());


#endif


  FEI_COUT << FEI_ENDL;
  if (returnValue != 0) {
    FEI_COUT << "at least 1 benchmark failed."<< FEI_ENDL << FEI_ENDL;
  }
  return(returnValue);
}

template<typename SET_TYPE>
double time_set_insert(int len)
{
  double start_time = fei::utils::cpu_time();

  SET_TYPE* set_objs = new SET_TYPE[len];

  int inner = 24;
  int outer = 8;

  int first = 9000;
  int inner_2 = inner/2;

  for(int n=0; n<len; ++n) {
    int col_n = first+n;
    SET_TYPE& set_ref = set_objs[n];

    for(int i=0; i<outer; ++i) {
      int col_i = col_n+i*outer;

      for(int j=0; j<inner_2; ++j) {
	set_ref.insert(col_i+j);
	set_ref.insert(col_i+j+inner_2);
      }
    }
  }

  delete [] set_objs;

  double elapsed_time = fei::utils::cpu_time() - start_time;
  return(elapsed_time);
}

template<typename SET_TYPE>
double time_set_insert2(int len)
{
  double start_time = fei::utils::cpu_time();

  SET_TYPE* set_objs = new SET_TYPE[len];

  int inner = 24;
  int outer = 8;

  int first = 9000;

  for(int n=0; n<len; ++n) {
    int col_n = first+n;
    SET_TYPE& set_ref = set_objs[n];

    for(int i=0; i<outer; ++i) {
      int col_i = col_n+i*outer;

      for(int j=0; j<inner/2; ++j) {
	set_ref.insert2(col_i+j);
	set_ref.insert2(col_i+j+inner/2);
      }
    }
  }

  delete [] set_objs;

  double elapsed_time = fei::utils::cpu_time() - start_time;
  return(elapsed_time);
}

template<typename SET_TYPE>
double benchmark_set()
{
  int len = 1200;

  //first find a len such that time_taken is at least 1 second
  double time_taken = time_set_insert<SET_TYPE>(len);
  while(time_taken < 1.0) {
    len *= 2;
    time_taken = time_set_insert<SET_TYPE>(len);
  }

  //now repeat until time_taken passes 5 seconds
  time_taken = 0.0;
  int i=0;
  while(time_taken<5.0) {
    time_taken += time_set_insert<SET_TYPE>(len);
    ++i;
  }

  return((double)(i*len)/time_taken);
}

template<typename SET_TYPE>
double benchmark_set2()
{
  int len = 1200;

  //first find a len such that time_taken is at least 1 second
  double time_taken = time_set_insert2<SET_TYPE>(len);
  while(time_taken < 1.0) {
    len *= 2;
    time_taken = time_set_insert2<SET_TYPE>(len);
  }

  //now repeat until time_taken passes 5 seconds
  time_taken = 0.0;
  int i=0;
  while(time_taken<5.0) {
    time_taken += time_set_insert2<SET_TYPE>(len);
    ++i;
  }

  return((double)(i*len)/time_taken);
}


int test_benchmarks::test2()
{
  int returnValue = 0;

  FEI_COUT<<FEI_ENDL;
  FEI_COUT << "Following group of benchmarks inserts integers into sorted lists"
	   << " (actually"<<FEI_ENDL
	   << "sets), which is a sub-task of the ragged-table benchmarks..."<<FEI_ENDL
	   << "A higher number is better."<<FEI_ENDL<<FEI_ENDL;

  print_benchmark_banner();

  double value, goldvalue;
  std::string passfail;
  std::string testname;


  value = benchmark_set<fei::ctg_set<int> >();
  testname = add_macro_values("fei::ctg_set");
  try {
    goldvalue = fei_test_utils::get_file_benchmark("./fei_utest_timings.txt",
					    testname.c_str());
    passfail = fei_test_utils::check_test_result(value, goldvalue, 10);
    if (passfail != "passed") returnValue = -1;
  }
  catch(...) {
    goldvalue = -1.0;
    passfail = " ";
  }

  print_benchmark_line("fei::ctg_set::insert", value, goldvalue, passfail.c_str());


#ifndef FEI_NO_STL_SET

  value = benchmark_set<std::set<int> >();
  goldvalue = -1.0;
  passfail = " ";
  print_benchmark_line("std::set::insert", value, goldvalue, passfail.c_str());

#endif

#ifdef FEI_HASH_SET

  value = benchmark_set<FEI_HASH_SET<int> >();
  goldvalue = -1.0;
  passfail = " ";
  print_benchmark_line("hash_set::insert", value, goldvalue, passfail.c_str());

#endif
  FEI_COUT << FEI_ENDL;
  FEI_COUT << "More list/set insertions..." << FEI_ENDL << FEI_ENDL;

  print_benchmark_banner();


  value = benchmark_set2<fei::ctg_set<int> >();
  testname = add_macro_values("fei::ctg_set2");
  try {
    goldvalue = fei_test_utils::get_file_benchmark("./fei_utest_timings.txt",
					    testname.c_str());
    passfail = fei_test_utils::check_test_result(value, goldvalue, 10);
    if (passfail != "passed") returnValue = -1;
  }
  catch(...) {
    goldvalue = -1.0;
    passfail = " ";
  }

  print_benchmark_line("fei::ctg_set::insert2", value, goldvalue, passfail.c_str());

  FEI_COUT << FEI_ENDL;
  if (returnValue != 0) {
    FEI_COUT << "at least 1 benchmark failed."<< FEI_ENDL << FEI_ENDL;
  }

  return(returnValue);
}

int test_benchmarks::test3()
{
  int len = 100000;
  int   n = 100000;

  std::vector<int> stdvector(len);

  std::vector<int> stdv_dest;

  int* stdvptr = &(stdvector[0]);

  for(int i=0; i<len; ++i) {
    stdvptr[i] = i*2;
  }

  FEI_COUT << FEI_ENDL << "time to perform " << n
           << " binary-searches and inserts on an std::vector" << FEI_ENDL
	   << " of length " << len << ": " << FEI_ENDL;

  double start_time = fei::utils::cpu_time();

  stdvector.reserve(n*2);

  std::vector<int>::iterator
    v_iter,
    v_beg = stdvector.begin(),
    v_end = stdvector.end();

  for(int k=0; k<n; ++k) {
    v_iter = std::lower_bound(v_beg, v_end, k*2);
    stdvector.insert(v_iter, k*2-1);
    v_beg = stdvector.begin();
    v_end = stdvector.end();
  }

  double elapsed_time = fei::utils::cpu_time() - start_time;

  FEI_COUT << elapsed_time << FEI_ENDL;

  return(0);
}

int test_benchmarks::test4()
{
  return(0);
}

int test_benchmarks::test5()
{
  return(0);
}

int test_benchmarks::test6()
{
  return(0);
}

int test_benchmarks::test7()
{
  return(0);
}

int test_benchmarks::test8()
{
  return(0);
}

