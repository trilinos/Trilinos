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

#include <vector>
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

  std::string s(id::key(val));

  std::cout << "ID type: " << id::name() << std::endl;
  std::cout << "String hash key (unique): " << s << std::endl;
  std::cout << "Int hash code (non-unique): " << id::hashCode(val) << std::endl;
  std::cout << "Is Teuchos hash key type: " << id::isHashKeyType() << std::endl;
  std::cout << "Is Teuchos Global Ordinal: " << id::isGlobalOrdinalType() << std::endl;
  std::cout << "Is Teuchos Packet type: " << id::isPacketType() << std::endl;
  std::cout << "Equal to self: " << id::equal(val, val) << std::endl;
  std::cout << "Equal to other: " << id::equal(val, compval) << std::endl;
  std::cout << "Key to original ID produces original ID: ";
    std::cout << id::equal(val, id::keyToGid(s)) << std::endl;
  std::cout << std::endl;
}
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int rank = session.getRank();

  if (rank == 0){
    char c='9';
    char c_other ='v';
    short int si=9;
    short int si_other =19;
    int i=9;
    int i_other=19;
    unsigned int ui=9;
    unsigned int ui_other=19;
    long int li=9;
    long int li_other=19;
    long unsigned int lui =9;
    long unsigned int lui_other =19;
    long long int lli = 3000000000;
    long long int lli_other = 3900000000;
    std::pair<int, int> pairVals(9, 9);
    std::pair<int, int> pairVals_other(9, 19);
    std::vector<int> vecVals(10, 9);
    std::vector<int> vecVals_other(10, 19);

    check_traits(c, c_other);
    check_traits(si, si_other);
    check_traits(i, i_other);
    check_traits(ui, ui_other);
    check_traits(li, li_other);
    check_traits(lui, lui_other);
    check_traits(lli, lli_other);
    check_traits(pairVals, pairVals_other);
    check_traits(vecVals, vecVals_other);
  }

  std::cout << "PASS" << std::endl;
}

