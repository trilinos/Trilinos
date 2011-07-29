// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

// TODO: doxygen comments

// TODO: just testing that the test compiles and looks reasonable
//   We need to test the validity of the values returned in check_traits.

#include <ostream>
#include <string>
#include <algorithm>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Teuchos_GlobalMPISession.hpp>   

using namespace std;

template <typename T> 
void check_traits(T &val, T &compval)
{
  typedef Z2::IdentifierTraits<T> id;

  double k(id::key(val));

  std::cout << "ID type: " << id::name() << std::endl;
  std::cout << "Hash key (unique): " << k << std::endl;
  std::cout << "Int hash code (non-unique): " << id::hashCode(val) << std::endl;
  std::cout << "Is Teuchos hash key type: " << id::isHashKeyType() << std::endl;
  std::cout << "Is Teuchos Global Ordinal: " << id::isGlobalOrdinalType() << std::endl;
  std::cout << "Is Teuchos Packet type: " << id::isPacketType() << std::endl;
  std::cout << "Equal to self: " << id::equal(val, val) << std::endl;
  std::cout << "Equal to other: " << id::equal(val, compval) << std::endl;
  std::cout << "Key to original ID produces original ID: ";
    std::cout << id::equal(val, id::keyToGid(k)) << std::endl;
  std::cout << std::endl;
}
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int rank = session.getRank();

  if (rank == 0){
    char c='a';
    char c_other ='v';
    short int si=1024;
    short int si_other =11024;
    int i=1024;
    int i_other=11024;
    unsigned int ui=1024;
    unsigned int ui_other=11024;
    long int li=1024;
    long int li_other=11024;
    long unsigned int lui =1024;
    long unsigned int lui_other =11024;
    long long int lli = 3000000000;
    long long int lli_other = 3102400000000;
    std::pair<int, int> pairVals(1024, 1024);
    std::pair<int, int> pairVals_other(1024, 11024);

    check_traits(c, c_other);
    check_traits(si, si_other);
    check_traits(i, i_other);
    check_traits(ui, ui_other);
    check_traits(li, li_other);
    check_traits(lui, lui_other);
    check_traits(lli, lli_other);
    check_traits(pairVals, pairVals_other);
  }

  std::cout << "PASS" << std::endl;
  return 0;
}

