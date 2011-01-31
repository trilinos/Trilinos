#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <stk_percept/Util.hpp>
#include <stk_adapt/SubDimCell.hpp>

#include <stk_percept/Hashtable.hpp>

namespace Teuchos
{
  template <> 
  inline
  int hashCode(const stk::adapt::SubDimCell<int>& x)
   {
     return (int)x.getHash();
   }

}

namespace stk
{
  namespace adapt
  {
    namespace unit_tests
    {
      using namespace std;

      TEST(SubDimCell, test0)
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

        std::ostringstream strs1;
        strs1 << set1;

        map_[set1] = 1;
        map_[set2] = 2;
        map_[set3] = 3;

        cout << "map_[set1]= " << map_[set1] << endl;
        cout << "map_[set2]= " << map_[set2] << endl;
        cout << "map_[set3]= " << map_[set3] << endl;

        EXPECT_EQ(map_[set1], map_[set2]);
        EXPECT_TRUE(set1 == set2);

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
#define OUT(expr) cout << QUOTE(expr) " = " << expr << endl
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

      
      TEST(SubDimCell, test1)
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

        std::ostringstream strs1;
        strs1 << set1;

        map_[set1] = 1;
        map_[set2] = 2;
        map_[set3] = 3;

        cout << "map_[set1]= " << map_[set1] << endl;
        cout << "map_[set2]= " << map_[set2] << endl;
        cout << "map_[set3]= " << map_[set3] << endl;

        EXPECT_EQ(map_[set1], map_[set2]);
        EXPECT_TRUE(set1 == set2);

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
#define OUT(expr) cout << QUOTE(expr) " = " << expr << endl
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

      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      

      TEST(SubDimCell, test2)
      {
        if (1) return;
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

        EXPECT_EQ(map_[set1], map_[set2]);
        EXPECT_TRUE(set1 == set2);

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

    }
  }
}
