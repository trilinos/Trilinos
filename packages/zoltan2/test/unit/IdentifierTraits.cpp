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

#include <vector>
#include <ostream>
#include <string>
#include <algorithm>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Teuchos_GlobalMPISession.hpp>   

template <typename T> 
void check_traits(T &val)
{
  std::string key = Z2::IdentifierTraits<T>::key(val);
  std::string name = Z2::IdentifierTraits<T>::name();

  std::cout << name << ", " << key << std::endl;
}
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int rank = session.getRank();

  if (rank == 0){
    char c='9';
    short int si=9;
    int i=9;
    unsigned int ui=9;
    long int li=9;
    long unsigned int lui =9;
    long long int lli = 3000000000;
    std::pair<int, int> pairVals(9, 9);
#ifdef SERIALIZATION_SUPPORTS_VECTORS
    std::vector<int> vecVals(10, 9);
#endif



    check_traits(c);
    check_traits(si);
    check_traits(i);
    check_traits(ui);
    check_traits(li);
    check_traits(lui);
    check_traits(lli);    // TODO only do this when HAVE_LONG_LONG
    check_traits(pairVals);
#ifdef SERIALIZATION_SUPPORTS_VECTORS
    check_traits(vecVals);
#endif

    // TODO: When we can test that a source file fails to compile, we can
    //       do this test.
    //double d=9.0;  
    //check_traits(d);

  }

  std::cout << "PASS" << std::endl;
}
