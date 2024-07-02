// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test tests the Anasazi BasicSort sort manager
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"

#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// finish this test

using namespace Teuchos;
using namespace Anasazi;

typedef double                    MT;
typedef ScalarTraits<MT>         MTTraits;

class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};

template <class MT>
bool checkValsLM(int howmany, const std::vector<MT> &vals) {
  // largest magnitude to smallest magnitude: |vals[i]| >= |vals[i+1]|
  for (int i=0; i<howmany-1; ++i) {
    if ( MTTraits::magnitude(vals[i]) < MTTraits::magnitude(vals[i+1]) ) return false;
  }
  return true;
}

template <class MT>
bool checkValsSM(int howmany, const std::vector<MT> &vals) {
  // smallest magnitude to largest magnitude: |vals[i]| <= |vals[i+1]|
  for (int i=0; i<howmany-1; ++i) {
    if ( MTTraits::magnitude(vals[i]) > MTTraits::magnitude(vals[i+1]) ) return false;
  }
  return true;
}

template <class MT>
bool checkValsLM(int howmany, const std::vector<MT> &r_vals, const std::vector<MT> &i_vals) {
  // largest magnitude to smallest magnitude: |vals[i]| >= |vals[i+1]|
  for (int i=0; i<howmany-1; ++i) {
    if ( r_vals[i]*r_vals[i]+i_vals[i]*i_vals[i] < r_vals[i+1]*r_vals[i+1]+i_vals[i+1]*i_vals[i+1] ) return false;
  }
  return true;
}

template <class MT>
bool checkValsSM(int howmany, const std::vector<MT> &r_vals, const std::vector<MT> &i_vals) {
  // smallest magnitude to largest magnitude: |vals[i]| <= |vals[i+1]|
  for (int i=0; i<howmany-1; ++i) {
    if ( r_vals[i]*r_vals[i]+i_vals[i]*i_vals[i] > r_vals[i+1]*r_vals[i+1]+i_vals[i+1]*i_vals[i+1] ) return false;
  }
  return true;
}

template <class MT>
bool checkValsLA(int howmany, const std::vector<MT> &vals) {
  // largest to smallest: vals[i] >= vals[i+1]
  for (int i=0; i<howmany-1; ++i) {
    if ( vals[i] < vals[i+1] ) return false;
  }
  return true;
}

template <class MT>
bool checkValsSA(int howmany, const std::vector<MT> &vals) {
  // smallest to largest: vals[i] <= vals[i+1]
  for (int i=0; i<howmany-1; ++i) {
    if ( vals[i] > vals[i+1] ) return false;
  }
  return true;
}

bool checkPermValid(int howmany, const std::vector<int> &perm) {
  std::vector<bool> hit(howmany,false);
  for (int i=0; i<howmany; i++) {
    if (perm[i] < 0 || perm[i] >= howmany) return false;
    hit[perm[i]] = true;
  }
  for (int i=0; i<howmany; i++) {
    if (hit[i] == false) return false;
  }
  // check that the rest are "unchanged"
  // in the tests, perm was always initialized so that perm[i] == i
  for (int i=howmany; i < (int) perm.size(); ++i) {
    if (perm[i] != i) return false;
  }
  return true;
}

template <class MT>
bool checkPermMatch(int howmany, const std::vector<int> &perm, const std::vector<MT> &orig, const std::vector<MT> &sorted) {
  for (int i=0; i< (int)perm.size(); ++i) {
    if ( orig[perm[i]] != sorted[i] ) return false;
  }
  return true;
}

int main(int argc, char *argv[])
{

  using std::endl;
  using std::string;
  using std::pair;
  using std::vector;
  using std::copy;
  using std::sort;
  using std::ostream_iterator;
  using Teuchos::null;
  using Teuchos::rcp;

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  bool debug = false;
  bool verbose = false;
  bool testFailed = false;
  int numVals = 10;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging info.");
  cmdp.setOption("numVals",&numVals,"Number of values for testing sorting.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (numVals < 0) {
    std::cerr << "--numVals must be >= 0. Setting to 10." << endl;
    numVals = 10;
  }

  //
  // Create an output manager
  RCP<OutputManager<MT> > printer 
    = rcp( new BasicOutputManager<MT>() );
  int verbosity = Errors;
  if (verbose || debug) {
    verbosity += Warnings;
  }
  if (debug) {
    verbosity += Debug;
  }
  printer->setVerbosity( verbosity );

  printer->stream(Warnings) << Anasazi_Version() << endl << endl;

  // 
  // seed random number generator
  int walltime = (int)Time::wallTime();
  printer->stream(Warnings) << "Seeding PRNG with " << walltime << endl << endl;
  MTTraits::seedrandom(walltime);

  // 
  // create the array of values to be sorted
  vector<MT> unsorted_r(numVals), unsorted_i(numVals);
  vector<int> pureperm(numVals);
  for (int i=0; i<numVals; ++i)
  {
    unsorted_r[i] = MTTraits::random();
    unsorted_i[i] = MTTraits::random();
    pureperm[i] = i;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  // Perform tests
  //
  try {

    // 
    // Create a sort manager with invalid argument: string
    {
      printer->print(Warnings,"*** Initialization with invalid sort string.\n");
      bool caught_expected_exception;
      try {
        BasicSort<MT> sorter("??");
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort accepted invalid sort string without throwing exception.");
    }

    // 
    // Create a sort manager with invalid argument: ParameterList
    {
      printer->print(Warnings,"*** Initialization with invalid sort string.\n");
      bool caught_expected_exception;
      try {
        Teuchos::ParameterList pl;
        pl.set("Sort Strategy","??");
        BasicSort<MT> sorter(pl);
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort accepted invalid sort string without throwing exception.");
    }

    // 
    // Try "LI" sort with real values
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with \"LI\".\n");
      BasicSort<MT> sorter("LI");
      vector<MT> vals(5,0.0);
      try {
        sorter.sort(vals);
        caught_expected_exception = false;
      }
      catch (const SortManagerError &sme) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted sort string \"LI\" without throwing exception.");
    }

    // 
    // Try "SI" sort with real values
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with \"SI\".\n");
      BasicSort<MT> sorter("SI");
      vector<MT> vals(5,0.0);
      try {
        sorter.sort(vals);
        caught_expected_exception = false;
      }
      catch (const SortManagerError &sme) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted sort string \"SI\" without throwing exception.");
    }

    // 
    // Try real values sort with too small value vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with too small value vector.\n");
      BasicSort<MT> sorter("SR");
      vector<MT> vals(5);
      try {
        sorter.sort(vals,Teuchos::null,6);
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted too small value vector without throwing exception.");
    }

    // 
    // Try real values sort with too small perm vector (specified)
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with too small perm vector.\n");
      BasicSort<MT> sorter("SR");
      vector<MT> vals(5);
      vector<int> perm(4);
      try {
        sorter.sort(vals,Teuchos::rcp(&perm,false),5);
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted too small perm vector without throwing exception.");
    }

    // 
    // Try real values sort with too small perm vector (default)
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with too small perm vector.\n");
      BasicSort<MT> sorter("SR");
      vector<MT> vals(5);
      vector<int> perm(4);
      try {
        sorter.sort(vals,Teuchos::rcp(&perm,false));
        caught_expected_exception = false;
      }
      catch (std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted too small perm vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small rvalue vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small rvalue vector.\n");
      BasicSort<MT> sorter("SR");
      vector<MT> rvals(5), ivals(6);
      try {
        sorter.sort(rvals,ivals,Teuchos::null,6);
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small rvalue vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small ivalue vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small ivalue vector.\n");
      BasicSort<MT> sorter("SR");
      vector<MT> rvals(6), ivals(5);
      try {
        sorter.sort(rvals,ivals,Teuchos::null,6);
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small ivalue vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small perm vector (specified)
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small perm vector.\n");
      BasicSort<MT> sorter("SR");
      vector<MT> rvals(6), ivals(6);
      vector<int> perm(5);
      try {
        sorter.sort(rvals,ivals,Teuchos::rcp(&perm,false),6);
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small perm vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small perm vector (default)
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small perm vector.\n");
      BasicSort<MT> sorter("SR");
      // min( rvals.length(), ivals.length() ) = min(7,6) = 6
      // this is the default sort length
      vector<MT> rvals(7), ivals(6);
      vector<int> perm(5);
      try {
        sorter.sort(rvals,ivals,Teuchos::rcp(&perm,false));
        caught_expected_exception = false;
      }
      catch (const std::invalid_argument &ia) {
        caught_expected_exception = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small perm vector without throwing exception.");
    }

    // sorter to use for rest of tests
    // the sort strategy will be modified; must be set to something in the meantime.
    BasicSort<MT> sorter("SR");

    // 
    // "SR" sorting test with real values
    {
      string which("SR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSA(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSA(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
    }

    // 
    // "LR" sorting test with real values
    {
      string which("LR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLA(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLA(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
    }

    // 
    // "SM" sorting test with real values
    {
      string which("SM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSM(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSM(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
    }

    // 
    // "LM" sorting test with real values
    {
      string which("LM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLM(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted(unsorted_r);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted.begin(), sorted.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLM(i,sorted) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
      }
    }

    // 
    // "SR" sorting test with complex values
    {
      string which("SR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSA(i,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSA(i,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
    }

    // 
    // "LR" sorting test with complex values
    {
      string which("LR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLA(i,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLA(i,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
    }

    // 
    // "SM" sorting test with complex values
    {
      string which("SM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSM(i,sorted_r,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSM(i,sorted_r,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
    }

    // 
    // "LM" sorting test with complex values
    {
      string which("LM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLM(i,sorted_r,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLM(i,sorted_r,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
    }

    // 
    // "SI" sorting test with complex values
    {
      string which("SI");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSA(i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsSA(i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
    }

    // 
    // "LI" sorting test with complex values
    {
      string which("LI");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      // try for each length without permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,null,i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLA(i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
      // try for each length with permutation
      for (int i=0; i<=numVals; ++i) {
        vector<MT> sorted_r(unsorted_r), sorted_i(unsorted_i);
        vector<int> perm(pureperm);
        printer->print(Debug,">> Before sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        sorter.sort(sorted_r,sorted_i,rcp(&perm,false),i);
        printer->print(Debug,">>  After sort: "); copy(sorted_r.begin(), sorted_r.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>      (imag): "); copy(sorted_i.begin(), sorted_i.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        printer->print(Debug,">>  Permutation: "); copy(perm.begin(), perm.end(), ostream_iterator<MT>(printer->stream(Debug), " ")); printer->print(Debug,"\n");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermValid(i,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_r,sorted_r) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkPermMatch(i,perm,unsorted_i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
        TEUCHOS_TEST_FOR_EXCEPTION( checkValsLA(i,sorted_i) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
      }
    }

  } // end of try
  catch (const get_out &go) {
    printer->stream(Warnings) << go.what() << endl;
    testFailed = true;
  }

  printer->print(Warnings,"\n");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (testFailed) {
    printer->print(Warnings,"End Result: TEST FAILED\n");
    return -1;
  }
  //
  // Default return value
  //
  printer->print(Warnings,"End Result: TEST PASSED\n");
  return 0;
}
