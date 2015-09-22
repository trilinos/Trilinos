// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// TODO: doxygen comments

// TODO: just testing that the test compiles and looks reasonable
//   We need to test the validity of the values returned in check_traits.

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_IdentifierTraits.hpp>
#include <ostream>
#include <string>
#include <algorithm>
#include <utility>
#include <Teuchos_GlobalMPISession.hpp>   

using namespace std;

// TODO - add answers to check_traits to verify the results
template <typename T> 
void check_traits(
  Teuchos::RCP<const Teuchos::Comm<int> > &tcomm , T &val, T &higherVal)
{
  const Teuchos::Comm<int> *comm = tcomm.get();
  typedef Zoltan2::IdentifierTraits<T> id;
  T valueList[2] = {val, higherVal};

  double k(0);
  if (id::hasUniqueKey())
    k = id::key(val);

  std::cout << "ID type: " << id::name() << std::endl;
  std::cout << "Value: " << id::stringify(val) << std::endl;
  std::cout << "\"Higher\" value: " << id::stringify(higherVal) << std::endl;
  std::cout << "Int hash code (non-unique): " << id::hashCode(val) << std::endl;
  std::cout << "Supports unique key: " << id::hasUniqueKey() << std::endl;
  if (id::hasUniqueKey())
    std::cout << "Key value: " << k << std::endl;
  std::cout << "Is Teuchos Global Ordinal: " << id::isGlobalOrdinal() << std::endl;
  std::cout << "Is valid id type: " << id::is_valid_id_type() << std::endl;

  //////////

  bool supported = true;
  T lo, hi;
  std::cout << "minMax: ";
  try{
    id::minMax(valueList, 2, lo, hi);
  }
  catch(...){
    supported = false;
  }

  if (supported)
    std::cout << id::stringify(lo) << ", " << id::stringify(hi) << std::endl;
  else
    std::cout << "not supported" << std::endl;

  //////////

  supported = true;
  T gmin, gmax;
  std::cout << "global minMax: ";
  try{
    id::globalMinMax(*comm, false, lo, hi, gmin, gmax);
  }
  catch(...){
    supported = false;
  }

  if (supported)
    std::cout << id::stringify(lo) << ", " << id::stringify(hi) << std::endl;
  else
    std::cout << "not supported" << std::endl;

  //////////

  supported = true;
  T diff;
  std::cout << "Difference : ";
  try{
    diff = id::difference(val, higherVal);
  }
  catch(...){
    supported = false;
  }

  if (supported)
    std::cout << id::stringify(diff) << std::endl;
  else
    std::cout << "not supported" << std::endl;

  //////////

  supported = true;
  bool consec=false; 
  std::cout << "Are consecutive: ";
  try{
    consec = id::areConsecutive(valueList, 2);
  }
  catch(...){
    supported = false;
  }

  if (supported)
    std::cout << consec << std::endl;
  else
    std::cout << "not supported" << std::endl;
  
  std::cout << std::endl;
}
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > nullComm;
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getDefaultSerialComm(nullComm);

  int rank = session.getRank();

  if (rank == 0){
    char c='a';
    char c_other ='v';
    unsigned char uc='A';
    unsigned char uc_other ='Z';
    short int si=1024;
    short int si_other =11024;
    unsigned short int usi=1024;
    unsigned short int usi_other =11024;
    int i=1024;
    int i_other=11024;
    unsigned int ui=1024;
    unsigned int ui_other=11024;
    long int li=1024;
    long int li_other=11024;
    long unsigned int lui =1024;
    long unsigned int lui_other =11024;
#ifdef HAVE_ZOLTAN2_LONG_LONG
    long long int lli = 3000000000;
    long long int lli_other = 3102400000000;
    unsigned long long int ulli = 3000000000;
    unsigned long long int ulli_other = 3102400000000;
#endif
    std::pair<int, int> pairVals(1024, 1024);
    std::pair<int, int> pairVals_other(11024, 11024);

    std::pair<long, long> pairValsLong(1024, 1024);
    std::pair<long, long> pairValsLong_other(11024, 11024);

    string strVal("right front wheel");
    string strVal_other("left front wheel");

    check_traits(comm, c, c_other);
    check_traits(comm, uc, uc_other);
    check_traits(comm, si, si_other);
    check_traits(comm, usi, usi_other);
    check_traits(comm, i, i_other);
    check_traits(comm, ui, ui_other);
    check_traits(comm, li, li_other);
    check_traits(comm, lui, lui_other);
#ifdef HAVE_ZOLTAN2_LONG_LONG
    check_traits(comm, lli, lli_other);
    check_traits(comm, ulli, ulli_other);
#endif
    check_traits(comm, pairVals, pairVals_other);
    check_traits(comm, pairValsLong, pairValsLong_other);
    check_traits(comm, strVal, strVal_other);

    std::cout << "PASS" << std::endl;
  }

  return 0;
}

