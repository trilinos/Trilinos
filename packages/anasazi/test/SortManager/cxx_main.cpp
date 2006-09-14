//@HEADER
// ************************************************************************
// 
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
//
//  This test tests the Anasazi BasicSort sort manager
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_LAPACK.hpp"

#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "MyMultiVec.hpp"
#include "MyOperator.hpp"

using namespace Teuchos;
using namespace Anasazi;
using namespace std;

typedef double                    ST;
typedef double                    MT;
typedef MultiVec<ST>              MV;
typedef Operator<ST>              OP;
typedef ScalarTraits<ST>         SCT;


class get_out : public logic_error {
  public: get_out(const string &whatarg) : logic_error(whatarg) {}
};

template <class MT>
class MySWO {
  public:
  MySWO(string which) : which_(which) {}
  bool operator()(MT p1, MT p2) {
    if ( !which_.compare("SM") ) {
      return SCT::magnitude(p1) < SCT::magnitude(p2);
    }
    if ( !which_.compare("LM") ) {
      return SCT::magnitude(p1) > SCT::magnitude(p2);
    }
    if ( !which_.compare("SR") ) {
      return p1 < p2;
    }
    if ( !which_.compare("LR") ) {
      return p1 > p2;
    }
    throw std::logic_error("logic error");
    return false;
  }
  private: 
  string which_;
};

template <class MT>
class MyPairSWO {
  public:
  MyPairSWO(string which) : which_(which) {}
  bool operator()(pair<MT,MT> p1, pair<MT,MT> p2) {
    LAPACK<int,MT> lapack;
    if ( !which_.compare("SM") ) {
      return lapack.LAPY2(p1.first,p1.second) < lapack.LAPY2(p2.first,p2.second);
    }
    if ( !which_.compare("LM") ) {
      return lapack.LAPY2(p1.first,p1.second) > lapack.LAPY2(p2.first,p2.second);
    }
    if ( !which_.compare("SR") ) {
      return p1.first < p2.first;
    }
    if ( !which_.compare("LR") ) {
      return p1.first > p2.first;
    }
    if ( !which_.compare("SI") ) {
      return p1.second < p2.second;
    }
    if ( !which_.compare("LI") ) {
      return p1.second > p2.second;
    }
    throw std::logic_error("logic error");
    return false;
  }
  private: 
  string which_;
};

namespace std {
  ostream &operator<<(ostream &out, const pair<MT,MT> &p) {
    out << "(" << p.first << ", " << p.second << ")";
    return out;
  }
}

bool checkPermValid(int howmany, const vector<int> &perm) {
  vector<bool> hit(howmany,false);
  for (int i=0; i<howmany; i++) {
    if (perm[i] < 0 || perm[i] >= howmany) return false;
    hit[perm[i]] = true;
  }
  for (int i=0; i<howmany; i++) {
    if (hit[i] == false) return false;
  }
  return true;
}

template <class MT>
bool checkPermMatch(int howmany, const vector<int> &perm, const vector<MT> &vref, const vector<MT> &vind) {
  for (int i=0; i<howmany; i++) {
    if ( vref[perm[i]] != vind[i] ) return false;
  }
  return true;
}

template <class MT>
bool checkValsMatch(const vector<MT> &vals1, const vector<MT> &vals2) {
  for (unsigned int i=0; i<vals1.size(); i++) {
    if ( vals1[i] != vals2[i] ) return false;
  }
  return true;
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  int MyPID;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#else 
  MyPID = 0;
#endif
  bool debug = false;
  bool verbose = false;
  bool testFailed = false;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging info.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  //
  // Create an output manager
  RefCountPtr<OutputManager<ST> > printer 
    = rcp( new BasicOutputManager<ST>() );
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
  SCT::seedrandom(walltime);

  ////////////////////////////////////////////////////////////////////////////////
  //
  // Perform tests
  //
  try {

    // 
    // Create a sort manager with invalid argument
    {
      printer->print(Warnings,"*** Initialization with invalid sort string.\n");
      bool caught_expected_exception;
      try {
        BasicSort<ST,MV,OP> sorter("??");
        caught_expected_exception = false;
      }
      catch (invalid_argument ia) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort accepted invalid sort string without throwing exception.");
    }

    // 
    // Try "LI" sort with real values
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with \"LI\".\n");
      BasicSort<ST,MV,OP> sorter("LI");
      vector<MT> vals(5);
      for (int i=0; i<5; i++) vals[i] = SCT::random();
      try {
        sorter.sort(NULL,vals.size(),vals);
        caught_expected_exception = false;
      }
      catch (SortManagerError sme) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted sort string \"LI\" without throwing exception.");
    }

    // 
    // Try "SI" sort with real values
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with \"SI\".\n");
      BasicSort<ST,MV,OP> sorter("SI");
      vector<MT> vals(5);
      try {
        sorter.sort(NULL,vals.size(),vals);
        caught_expected_exception = false;
      }
      catch (SortManagerError sme) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted sort string \"SI\" without throwing exception.");
    }

    // 
    // Try real values sort with too small value vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with too small value vector.\n");
      BasicSort<ST,MV,OP> sorter("SR");
      vector<MT> vals(5);
      try {
        sorter.sort(NULL,6,vals);
        caught_expected_exception = false;
      }
      catch (invalid_argument ia) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted too small value vector without throwing exception.");
    }

    // 
    // Try real values sort with too small perm vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort real values with too small perm vector.\n");
      BasicSort<ST,MV,OP> sorter("SR");
      vector<MT> vals(5);
      vector<int> perm(4);
      try {
        sorter.sort(NULL,5,vals,&perm);
        caught_expected_exception = false;
      }
      catch (invalid_argument ia) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(real) accepted too small perm vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small rvalue vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small rvalue vector.\n");
      BasicSort<ST,MV,OP> sorter("SR");
      vector<MT> rvals(5), ivals(6);
      try {
        sorter.sort(NULL,6,rvals,ivals);
        caught_expected_exception = false;
      }
      catch (invalid_argument ia) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small rvalue vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small ivalue vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small ivalue vector.\n");
      BasicSort<ST,MV,OP> sorter("SR");
      vector<MT> rvals(6), ivals(5);
      try {
        sorter.sort(NULL,6,rvals,ivals);
        caught_expected_exception = false;
      }
      catch (invalid_argument ia) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small ivalue vector without throwing exception.");
    }

    // 
    // Try complex values sort with too small perm vector
    {
      bool caught_expected_exception;
      printer->print(Warnings,"*** Invalid: Sort complex values with too small perm vector.\n");
      BasicSort<ST,MV,OP> sorter("SR");
      vector<MT> rvals(6), ivals(6);
      vector<int> perm(5);
      try {
        sorter.sort(NULL,6,rvals,ivals,&perm);
        caught_expected_exception = false;
      }
      catch (invalid_argument ia) {
        caught_expected_exception = true;
      }
      TEST_FOR_EXCEPTION(caught_expected_exception == false,get_out,"BasicSort::sort(complex) accepted too small perm vector without throwing exception.");
    }

    // sorter to use for rest of tests
    BasicSort<ST,MV,OP> sorter("SR");

    // 
    // "SR" sorting test with real values
    {
      string which("SR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      MySWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<MT> valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) valsorig[i] = SCT::random();
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsbsm(valsorig);
      sorter.sort(NULL,5,valsbsm,&perm);
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
    }

    // 
    // "SM" sorting test with real values
    {
      string which("SM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      MySWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<MT> valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) valsorig[i] = SCT::random();
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsbsm(valsorig);
      sorter.sort(NULL,5,valsbsm,&perm);
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
    }

    // 
    // "LR" sorting test with real values
    {
      string which("LR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      MySWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<MT> valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) valsorig[i] = SCT::random();
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsbsm(valsorig);
      sorter.sort(NULL,5,valsbsm,&perm);
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
    }

    // 
    // "LM" sorting test with real values
    {
      string which("LM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort real values by \"" << which << "\"" << endl;
      MySWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<MT> valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) valsorig[i] = SCT::random();
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<MT> valsbsm(valsorig);
      sorter.sort(NULL,5,valsbsm,&perm);
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<MT>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(real) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(real) returned incorrect sort.");
    }

    // 
    // "SR" sorting test with complex values
    {
      string which("SR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      MyPairSWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<pair<MT,MT> > valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) { valsorig[i].first = SCT::random(); valsorig[i].second = SCT::random(); }
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsbsm(valsorig);
      {
        vector<MT> bsm_indr(valsbsm.size()), bsm_indi(valsbsm.size());
        for (unsigned int i=0; i<valsbsm.size(); i++) { bsm_indr[i] = valsbsm[i].first; bsm_indi[i] = valsbsm[i].second; }
        sorter.sort(NULL,5,bsm_indr,bsm_indi,&perm);
        for (unsigned int i=0; i<valsbsm.size(); i++) { valsbsm[i].first = bsm_indr[i]; valsbsm[i].second = bsm_indi[i]; }
      }
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
    }

    // 
    // "SM" sorting test with complex values
    {
      string which("SM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      MyPairSWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<pair<MT,MT> > valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) { valsorig[i].first = SCT::random(); valsorig[i].second = SCT::random(); }
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsbsm(valsorig);
      {
        vector<MT> bsm_indr(valsbsm.size()), bsm_indi(valsbsm.size());
        for (unsigned int i=0; i<valsbsm.size(); i++) { bsm_indr[i] = valsbsm[i].first; bsm_indi[i] = valsbsm[i].second; }
        sorter.sort(NULL,5,bsm_indr,bsm_indi,&perm);
        for (unsigned int i=0; i<valsbsm.size(); i++) { valsbsm[i].first = bsm_indr[i]; valsbsm[i].second = bsm_indi[i]; }
      }
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
    }

    // 
    // "SI" sorting test with complex values
    {
      string which("SI");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      MyPairSWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<pair<MT,MT> > valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) { valsorig[i].first = SCT::random(); valsorig[i].second = SCT::random(); }
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsbsm(valsorig);
      {
        vector<MT> bsm_indr(valsbsm.size()), bsm_indi(valsbsm.size());
        for (unsigned int i=0; i<valsbsm.size(); i++) { bsm_indr[i] = valsbsm[i].first; bsm_indi[i] = valsbsm[i].second; }
        sorter.sort(NULL,5,bsm_indr,bsm_indi,&perm);
        for (unsigned int i=0; i<valsbsm.size(); i++) { valsbsm[i].first = bsm_indr[i]; valsbsm[i].second = bsm_indi[i]; }
      }
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
    }

    // 
    // "LR" sorting test with complex values
    {
      string which("LR");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      MyPairSWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<pair<MT,MT> > valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) { valsorig[i].first = SCT::random(); valsorig[i].second = SCT::random(); }
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsbsm(valsorig);
      {
        vector<MT> bsm_indr(valsbsm.size()), bsm_indi(valsbsm.size());
        for (unsigned int i=0; i<valsbsm.size(); i++) { bsm_indr[i] = valsbsm[i].first; bsm_indi[i] = valsbsm[i].second; }
        sorter.sort(NULL,5,bsm_indr,bsm_indi,&perm);
        for (unsigned int i=0; i<valsbsm.size(); i++) { valsbsm[i].first = bsm_indr[i]; valsbsm[i].second = bsm_indi[i]; }
      }
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
    }

    // 
    // "LM" sorting test with complex values
    {
      string which("LM");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      MyPairSWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<pair<MT,MT> > valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) { valsorig[i].first = SCT::random(); valsorig[i].second = SCT::random(); }
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsbsm(valsorig);
      {
        vector<MT> bsm_indr(valsbsm.size()), bsm_indi(valsbsm.size());
        for (unsigned int i=0; i<valsbsm.size(); i++) { bsm_indr[i] = valsbsm[i].first; bsm_indi[i] = valsbsm[i].second; }
        sorter.sort(NULL,5,bsm_indr,bsm_indi,&perm);
        for (unsigned int i=0; i<valsbsm.size(); i++) { valsbsm[i].first = bsm_indr[i]; valsbsm[i].second = bsm_indi[i]; }
      }
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
    }

    // 
    // "LI" sorting test with complex values
    {
      string which("LI");
      sorter.setSortType(which);
      printer->stream(Warnings) << "*** Sort complex values by \"" << which << "\"" << endl;
      MyPairSWO<MT> myswo(which);   // StrictWeakOrdering for use with std::sort
      vector<pair<MT,MT> > valsorig(10);      // values to be sorted
      vector<int> perm(8);      // permutation vector
      for (unsigned int i=0; i<valsorig.size(); i++) { valsorig[i].first = SCT::random(); valsorig[i].second = SCT::random(); }
      printer->print(Debug,">> Original vals                       : ");
      copy(valsorig.begin(), valsorig.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsstl(valsorig);
      sort(valsstl.begin(),valsstl.begin()+5,myswo);
      printer->print(Debug,">> After std::sort of first five       : ");
      copy(valsstl.begin(), valsstl.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      vector<pair<MT,MT> > valsbsm(valsorig);
      {
        vector<MT> bsm_indr(valsbsm.size()), bsm_indi(valsbsm.size());
        for (unsigned int i=0; i<valsbsm.size(); i++) { bsm_indr[i] = valsbsm[i].first; bsm_indi[i] = valsbsm[i].second; }
        sorter.sort(NULL,5,bsm_indr,bsm_indi,&perm);
        for (unsigned int i=0; i<valsbsm.size(); i++) { valsbsm[i].first = bsm_indr[i]; valsbsm[i].second = bsm_indi[i]; }
      }
      printer->print(Debug,">> After BasicSort::sort of first five : ");
      copy(valsbsm.begin(), valsbsm.end(), ostream_iterator<pair<MT,MT> >(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      printer->print(Debug,">> BasicSort::sort returns perm        : ");
      copy(perm.begin(), perm.end(), ostream_iterator<int>(printer->stream(Debug), " "));
      printer->print(Debug,"\n");
      TEST_FOR_EXCEPTION( checkPermValid(5,perm) == false, get_out, "BasicSort::sort(complex) returned invalid permutation vector.");
      TEST_FOR_EXCEPTION( checkPermMatch(5,perm,valsorig,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect permutation vector.");
      TEST_FOR_EXCEPTION( checkValsMatch(valsstl,valsbsm) == false, get_out, "BasicSort::sort(complex) returned incorrect sort.");
    }

  } // end of try
  catch (get_out go) {
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
