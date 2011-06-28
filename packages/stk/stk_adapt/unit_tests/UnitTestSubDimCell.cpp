#include <stk_percept/Util.hpp>
#include <stk_percept/Hashtable.hpp>

#include <stk_adapt/SubDimCell.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#define USE_SPARSEHASH 0
#if USE_SPARSEHASH
#include <google/sparse_hash_map>
#include <google/dense_hash_map>
#endif

#include <stdexcept>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

namespace Teuchos {

template <>
inline
int hashCode(const stk::adapt::SubDimCell<int>& x)
{
  return (int)x.getHash();
}

} // namespace Teuchos

namespace stk {
namespace adapt {
namespace unit_tests {

using namespace std;

#define OUT(expr) cout << QUOTE(expr) " = " << expr << endl

#if USE_SPARSEHASH

STKUNIT_UNIT_TEST(SubDimCell, test_google_dense_hash_map)
{
  google::dense_hash_map<SubDimCell<int>, int, my_hash<int, 4>, my_equal_to<int, 4> > map_;

  SubDimCell<int> empty_key;
  empty_key.insert( std::numeric_limits<int>::max() );
  map_.set_empty_key(empty_key);

  SubDimCell<int> deleted_key;
  deleted_key.insert( std::numeric_limits<int>::max()-1 );
  map_.set_deleted_key(deleted_key);

  //percept::Hashtable<SubDimCell<int>, int, less<int> >::iterator map_it;
  SubDimCell<int> set1;
  SubDimCell<int> set2;
  SubDimCell<int> set3;

  set1.insert(1);  set1.insert(2);  set1.insert(3);
  set2.insert(3);  set2.insert(1);  set2.insert(2);
  set3.insert(4);  set3.insert(1);  set3.insert(2);
  typedef pair<SubDimCell<int>::iterator, bool> result_type;

  std::ostringstream strs1;
  strs1 << set1;

  map_[set1] = 1;
  map_[set2] = 2;
  map_[set3] = 3;

  cout << "map_[set1]= " << map_[set1] << endl;
  cout << "map_[set2]= " << map_[set2] << endl;
  cout << "map_[set3]= " << map_[set3] << endl;

  STKUNIT_EXPECT_EQ(map_[set1], map_[set2]);
  STKUNIT_EXPECT_TRUE(set1 == set2);

  cout << "set1 = " << endl;
  for (SubDimCell<int>::iterator si1 = set1.begin(); si1 != set1.end(); si1++)
  {
    cout <<  *si1 << " ";
  }
  cout << endl;
  cout << "set2 = " << endl;
  for (SubDimCell<int>::iterator si2 = set2.begin(); si2 != set2.end(); si2++)
  {
    cout <<  *si2 << " ";
  }
  cout << endl;
  cout << "set3 = " << endl;
  for (SubDimCell<int>::iterator si3 = set3.begin(); si3 != set3.end(); si3++)
  {
    cout <<  *si3 << " ";
  }
  cout << endl;

  //#define QUOTE(expr) #expr
  struct less<SubDimCell<int> > less_set;
  OUT(less_set(set1,set2));
  OUT(less_set(set2,set1));
  OUT(less_set(set1,set3));
  OUT(less_set(set3,set1));

  copy(set1.begin(), set1.end(), ostream_iterator<int>(cout, " "));
  cout << endl;
  //cout << set1 << endl;
  cout << set1 << endl;
  vector<int> aaa1(3);
  aaa1[0] = 0;
  aaa1[1] = 1;
  aaa1[2] = 2;
  cout << aaa1 << endl;
  //aaa(cout, aaa1);
  //cout << map_ << endl;

  map_.clear();
}
#endif

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(SubDimCell, test_percept_hashtable)
{
  percept::Hashtable<SubDimCell<int>, int> map_;
  //percept::Hashtable<SubDimCell<int>, int, less<int> >::iterator map_it;
  SubDimCell<int> set1;
  SubDimCell<int> set2;
  SubDimCell<int> set3;

  set1.insert(1);  set1.insert(2);  set1.insert(3);
  set2.insert(3);  set2.insert(1);  set2.insert(2);
  set3.insert(4);  set3.insert(1);  set3.insert(2);
  typedef pair<SubDimCell<int>::iterator, bool> result_type;

  cout << "set1 = " << set1 << " hashCode(set1) = " << Teuchos::hashCode(set1) << endl;
  cout << "set2 = " << set2 << " hashCode(set2) = " << Teuchos::hashCode(set2) << endl;
  cout << "set3 = " << set3 << " hashCode(set3) = " << Teuchos::hashCode(set3) << endl;

  std::ostringstream strs1;
  strs1 << set1;

#if 0
  map_[set1] = 1;
  map_[set2] = 2;
  map_[set3] = 3;
#else
  map_.put(set1, 1, true);
  map_.put(set2, 2, true);
  map_.put(set3, 3, true);
  cout << "map_= " << map_ << endl;
#endif
  cout << "map_[set1]= " << map_[set1] << endl;
  cout << "map_[set2]= " << map_[set2] << endl;
  cout << "map_[set3]= " << map_[set3] << endl;

  cout << "map_.contains(set1) = " << map_.containsKey(set1, true)  << endl;
  cout << "map_.contains(set2) = " << map_.containsKey(set2, true)  << endl;
  cout << "map_.contains(set3) = " << map_.containsKey(set3, true)  << endl;

  STKUNIT_EXPECT_EQ(map_[set1], map_[set2]);
  STKUNIT_EXPECT_TRUE(set1 == set2);

  cout << "set1 = " << endl;
  for (SubDimCell<int>::iterator si1 = set1.begin(); si1 != set1.end(); si1++)
  {
    cout <<  *si1 << " ";
  }
  cout << endl;
  cout << "set2 = " << endl;
  for (SubDimCell<int>::iterator si2 = set2.begin(); si2 != set2.end(); si2++)
  {
    cout <<  *si2 << " ";
  }
  cout << endl;
  cout << "set3 = " << endl;
  for (SubDimCell<int>::iterator si3 = set3.begin(); si3 != set3.end(); si3++)
  {
    cout <<  *si3 << " ";
  }
  cout << endl;

  //#define QUOTE(expr) #expr
  struct less<SubDimCell<int> > less_set;
  OUT(less_set(set1,set2));
  OUT(less_set(set2,set1));
  OUT(less_set(set1,set3));
  OUT(less_set(set3,set1));

  copy(set1.begin(), set1.end(), ostream_iterator<int>(cout, " "));
  cout << endl;
  //cout << set1 << endl;
  cout << set1 << endl;
  vector<int> aaa1(3);
  aaa1[0] = 0;
  aaa1[1] = 1;
  aaa1[2] = 2;
  cout << aaa1 << endl;
  //aaa(cout, aaa1);
  cout << map_ << endl;
}

//=============================================================================
//=============================================================================
//=============================================================================


STKUNIT_UNIT_TEST(SubDimCell, test1)
{
  map<SubDimCell<int>, int> map_;
  map<SubDimCell<int>, int, less<int> >::iterator map_it;
  SubDimCell<int> set1;
  SubDimCell<int> set2;
  SubDimCell<int> set3;

  set1.insert(1);  set1.insert(2);  set1.insert(3);
  set2.insert(3);  set2.insert(1);  set2.insert(2);
  set3.insert(4);  set3.insert(1);  set3.insert(2);
  typedef pair<SubDimCell<int>::iterator, bool> result_type;


  std::cout << "tmp srk set1= " << set1[0] << " " << set1[1] << " " << set1[2] << std::endl;

  std::ostringstream strs1;
  strs1 << set1;

  map_[set1] = 1;
  map_[set2] = 2;
  map_[set3] = 3;

  cout << "map_[set1]= " << map_[set1] << endl;
  cout << "map_[set2]= " << map_[set2] << endl;
  cout << "map_[set3]= " << map_[set3] << endl;

  STKUNIT_EXPECT_EQ(map_[set1], map_[set2]);
  STKUNIT_EXPECT_TRUE(set1 == set2);

  cout << "set1 = " << endl;
  for (SubDimCell<int>::iterator si1 = set1.begin(); si1 != set1.end(); si1++)
  {
    cout <<  *si1 << " ";
  }
  cout << endl;
  cout << "set2 = " << endl;
  for (SubDimCell<int>::iterator si2 = set2.begin(); si2 != set2.end(); si2++)
  {
    cout <<  *si2 << " ";
  }
  cout << endl;
  cout << "set3 = " << endl;
  for (SubDimCell<int>::iterator si3 = set3.begin(); si3 != set3.end(); si3++)
  {
    cout <<  *si3 << " ";
  }
  cout << endl;

  //#define QUOTE(expr) #expr
  struct less<SubDimCell<int> > less_set;
  OUT(less_set(set1,set2));
  OUT(less_set(set2,set1));
  OUT(less_set(set1,set3));
  OUT(less_set(set3,set1));

  copy(set1.begin(), set1.end(), ostream_iterator<int>(cout, " "));
  cout << endl;
  //cout << set1 << endl;
  cout << set1 << endl;
  vector<int> aaa1(3);
  aaa1[0] = 0;
  aaa1[1] = 1;
  aaa1[2] = 2;
  cout << aaa1 << endl;
  //aaa(cout, aaa1);
  //cout << map_ << endl;
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(SubDimCell, test2)
{
#if 0
  map<SubDimCell<int>, int> map_;
  map<SubDimCell<int>, int, less<int> >::iterator map_it;
  SubDimCell<int> set1;
  SubDimCell<int> set2;
  SubDimCell<int> set3;
  set1.insert(1);  set1.insert(2);  set1.insert(3);
  set2.insert(3);  set2.insert(1);  set2.insert(2);
  set3.insert(4);  set3.insert(1);  set3.insert(2);
  typedef pair<SubDimCell<int>::iterator, bool> result_type;

  map_[set1] = 1;
  map_[set2] = 2;
  map_[set3] = 3;

  cout << "map_[set1]= " << map_[set1] << endl;
  cout << "map_[set2]= " << map_[set2] << endl;
  cout << "map_[set3]= " << map_[set3] << endl;

  STKUNIT_EXPECT_EQ(map_[set1], map_[set2]);
  STKUNIT_EXPECT_TRUE(set1 == set2);

  cout << "set1 = " << endl;
  for (SubDimCell<int>::iterator si1 = set1.begin(); si1 != set1.end(); si1++)
  {
    cout <<  *si1 << " ";
  }
  cout << endl;
  cout << "set2 = " << endl;
  for (SubDimCell<int>::iterator si2 = set2.begin(); si2 != set2.end(); si2++)
  {
    cout <<  *si2 << " ";
  }
  cout << endl;
  cout << "set3 = " << endl;
  for (SubDimCell<int>::iterator si3 = set3.begin(); si3 != set3.end(); si3++)
  {
    cout <<  *si3 << " ";
  }
  cout << endl;

  //#define QUOTE(expr) #expr

  struct less<SubDimCell<int> > less_set;
  OUT(less_set(set1,set2));
  OUT(less_set(set2,set1));
  OUT(less_set(set1,set3));
  OUT(less_set(set3,set1));

  copy(set1.begin(), set1.end(), ostream_iterator<int>(cout, " "));
  cout << endl;
  //cout << set1 << endl;
  cout << set1 << endl;
  vector<int> aaa1(3);
  aaa1[0] = 0;
  aaa1[1] = 1;
  aaa1[2] = 2;
  cout << aaa1 << endl;
  //aaa(cout, aaa1);
  //cout << map_ << endl;
#endif
}

}
}
}
