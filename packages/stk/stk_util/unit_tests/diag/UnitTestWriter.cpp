/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <boost/regex.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>

#include <utility>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <bitset>
#include <memory>

#include <stk_util/util/IndentStreambuf.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterManip.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

using namespace stk::diag;

enum LogMask {
  LOG_ALWAYS		= stk::LOG_ALWAYS,
  LOG_TRACE		= stk::LOG_TRACE,
  LOG_TRACE_STATS	= stk::LOG_TRACE_STATS,
  LOG_TRACE_SUB_CALLS	= stk::LOG_TRACE_SUB_CALLS,
  LOG_MEMBERS		= stk::LOG_MEMBERS,

  LOG_TEST1             = 0x0000010,
  LOG_TEST2             = 0x0000020
};


std::ostringstream &
oss()
{
  static std::ostringstream     s_oss;

  return s_oss;
}


std::ostream &
dwout()
{
  static stk::indent_streambuf s_dwoutStreambuf(oss().rdbuf());
  static std::ostream s_dwout(&s_dwoutStreambuf);
  
  return s_dwout;
}


// Diagnostic writer
stk::diag::Writer &
dw()
{
  static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);

  return s_diagWriter;
}

STKUNIT_UNIT_TEST(UnitTestWriter, UnitTest)
{
  dw() << "This is a test" << dendl << dflush;

  STKUNIT_ASSERT_EQUAL(std::string("This is a test\n"), oss().str());
  
  oss().str("");
  dw() << "Level 0" << push << dendl
       << "This is a test" << dendl
       << pop << dendl;

  STKUNIT_ASSERT_EQUAL(std::string("Level 0 {\n  This is a test\n}\n"), oss().str());

  oss().str("");
  {
    stk::diag::WriterThrowSafe throw_safe__(dw());
    
    dw() << "Level 0" << push << dendl
         << "Level 1" << push << dendl
         << "Level 2" << push << dendl
         << "Level 3" << push << dendl
         << "This is a test" << dendl
         << pop << dendl;  
  }
  dw() << dendl;
  STKUNIT_ASSERT_EQUAL(std::string("Level 0 {\n  Level 1 {\n    Level 2 {\n      Level 3 {\n        This is a test\n      }\n    }\n  }\n}\n"), oss().str());

  oss().str("");
  {
    int                 x1 = 7;
    unsigned            x2 = 7;
    short               x3 = 7;
    unsigned short      x4 = 7;
    long                x5 = 7;
    unsigned long       x6 = 7;
    long long           x7 = 7;
    unsigned long long  x8 = 7;
    char                x9 = '\x7';
    signed char         x10 = '\x7';
    unsigned char       x11 = '\x7';
    
    dw() << x1 << dendl;
    dw() << x2 << dendl;
    dw() << x3 << dendl;
    dw() << x4 << dendl;
    dw() << x5 << dendl;
    dw() << x6 << dendl;
    dw() << x7 << dendl;
    dw() << x8 << dendl;
    dw() << x9 << dendl;
    dw() << x10 << dendl;
    dw() << x11 << dendl;
    dw() << 7.0 << dendl;
    dw() << 7.0f << dendl;
    dw() << 7.0l << dendl;
    dw() << "This is a test" << dendl;
    dw() << std::string("This is a test") << dendl;    
  }
  STKUNIT_ASSERT_EQUAL(std::string("7\n7\n7\n7\n7\n7\n7\n7\n7\n7\n7\n7\n7\n7\nThis is a test\nThis is a test\n"), oss().str());
  
  oss().str("");
  dw() << std::hex << 16 << dendl;
  dw() << std::oct << 16 << dendl;
  dw() << std::dec << 16 << dendl;
  dw() << std::fixed << 3.14159265 << dendl;
  dw() << std::scientific << 3.14159265 << dendl;
  dw() << stk::diag::hex << 16 << dendl;
  dw() << stk::diag::oct << 16 << dendl;
  dw() << stk::diag::dec << 16 << dendl;
  dw() << stk::diag::fixed << 3.14159265 << dendl;
  dw() << stk::diag::scientific << 3.14159265 << dendl;
  dw() << stk::diag::setw(5) << 3.14159265 << dendl;
  dw() << stk::diag::setprecision(5) << 3.14159265 << dendl;
  dw() << stk::diag::setiosflags(std::ios::fixed) << 3.14159265 << dendl;
  dw() << stk::diag::resetiosflags(std::ios::fixed) << 3.14159265 << dendl;
  dw() << stk::diag::setfill('#') << stk::diag::setw(10) << "x" << dendl;
  STKUNIT_ASSERT_EQUAL(std::string("10\n20\n16\n3.141593\n3.141593e+00\n10\n20\n16\n3.141593\n3.141593e+00\n3.141593e+00\n3.14159e+00\n3.1416\n3.14159e+00\n#########x\n"), oss().str());

  oss().str("");
  {
    std::vector<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    
    std::vector<int> vl;
    for (int i = 0; i < 20; ++i)
      vl.push_back(i);
    
    std::vector<int *> vp;
    vp.push_back(new int(1));
    vp.push_back(new int(2));
    vp.push_back(new int(3));
    vp.push_back(0);
    
    std::list<int> l;
    l.push_back(1);
    l.push_back(2);
    l.push_back(3);
    
    std::list<int *> lp;
    lp.push_back(new int(1));
    lp.push_back(new int(2));
    lp.push_back(new int(3));
    
    std::map<int, int> m;
    m[1] = 2;
    m[2] = 3;
    m[3] = 4;
    
    std::map<int, int *> mp;
    mp[1] = new int(2);
    mp[2] = new int(3);
    mp[3] = new int(4);

    std::multimap<int, int> mm;
    mm.insert(std::multimap<int, int>::value_type(1, 2));
    mm.insert(std::multimap<int, int>::value_type(1, 3));
    mm.insert(std::multimap<int, int>::value_type(2, 4));
    
    std::multimap<int, int *> mmp;
    mmp.insert(std::multimap<int, int *>::value_type(1, new int(2)));
    mmp.insert(std::multimap<int, int *>::value_type(1, new int(3)));
    mmp.insert(std::multimap<int, int *>::value_type(2, new int(4)));

    std::set<int> s;
    s.insert(2);
    s.insert(3);
    s.insert(4);
    
    std::set<int *> sp;
    sp.insert(new int(2));
    sp.insert(new int(3));
    sp.insert(new int(4));

    std::multiset<int> ms;
    ms.insert(2);
    ms.insert(2);
    ms.insert(4);
    
    std::multiset<int *> msp;
    msp.insert(new int(2));
    msp.insert(new int(3));
    msp.insert(new int(4));

    std::bitset<8> b;
    b[1] = 1;
    b[3] = 1;
    
    dw() << typeid(int) << dendl;
    dw() << std::pair<int, int>(5, 7) << dendl;
    dw() << v << dendl;
    dw() << vl << dendl;
    dw() << vp << dendl;
    dw() << l << dendl;
    dw() << lp << dendl;
    dw() << m << dendl;
    dw() << mp << dendl;
    dw() << mm << dendl;
    dw() << mmp << dendl;
    dw() << s << dendl;
    dw() << sp << dendl;
    dw() << ms << dendl;
    dw() << msp << dendl;
    dw() << b << dendl;

    for ( std::multimap<int , int *>::iterator curmmp = mmp.begin() ; curmmp != mmp.end() ; curmmp++ )
      delete curmmp->second;
    for ( std::map<int,int *>::iterator curmp = mp.begin() ; curmp != mp.end() ; curmp++ )
      delete curmp->second;
    for ( std::list<int *>::iterator curp = lp.begin() ; curp != lp.end() ; curp++ )
      delete *curp;
    for ( std::set<int *>::iterator cursp = sp.begin() ; cursp != sp.end() ; cursp++ )
      delete *cursp;
    for ( std::multiset<int *>::iterator curmsp = msp.begin() ; curmsp != msp.end() ; curmsp++ )
      delete *curmsp;
    
  }
//  STKUNIT_ASSERT_EQUAL(std::string("int\n(5:7)\nstd::vector<int, std::allocator<int> >, size 3 {\n  1 2 3 \n}\nstd::vector<int*, std::allocator<int*> >, size 3 {\n  [0] (pointer 0x53b040), 1\n  [1] (pointer 0x53b770), 2\n  [2] (pointer 0x53b750), 3\n}\nstd::list<int, std::allocator<int> >, size 3 {\n  [0] 1\n  [1] 2\n  [2] 3\n}\nstd::list<int*, std::allocator<int*> >, size 3 {\n  [0] (pointer 0x53b820), 1\n  [1] (pointer 0x53b8a0), 2\n  [2] (pointer 0x53b8e0), 3\n}\nstd::map<int, int, std::less<int>, std::allocator<std::pair<int const, int> > >, size 3 {\n  [1] 2\n  [2] 3\n  [3] 4\n}\nstd::map<int, int*, std::less<int>, std::allocator<std::pair<int const, int*> > >, size 3 {\n  [1] 0x53b9f0\n  [2] 0x53ba50\n  [3] 0x53bab0\n}\nstd::multimap<int, int, std::less<int>, std::allocator<std::pair<int const, int> > >, size 3 {\n  [1] 2\n  [1] 3\n  [2] 4\n}\nstd::multimap<int, int*, std::less<int>, std::allocator<std::pair<int const, int*> > >, size 3 {\n  [1] 0x53bb60\n  [1] 0x53bbc0\n  [2] 0x53bc20\n}\nstd::set<int, std::less<int>, std::allocator<int> >, size 3 {\n  2\n  3\n  4\n}\nstd::set<int*, std::less<int*>, std::allocator<int*> >, size 3 {\n  0x53bd10\n  0x53bd60\n  0x53bdb0\n}\nstd::multiset<int, std::less<int>, std::allocator<int> >, size 3 {\n  2\n  2\n  4\n}\nstd::multiset<int*, std::less<int*>, std::allocator<int*> >, size 3 {\n  0x53be90\n  0x53bee0\n  0x53bf30\n}\n00001010\n"), oss().str());  

  oss().str("");
  {
    std::auto_ptr<int> a0(new int(1));
    std::auto_ptr<int> a1;
    dw() << a0 << dendl;
    dw() << a1 << dendl;
    a1 = a0;
    dw() << a0 << dendl;
    dw() << a1 << dendl;    
  }
//  STKUNIT_ASSERT_EQUAL(std::string("std::auto_ptr<int>, 0x53b8c0, 1\n std::auto_ptr<int>, <not created or not owner>\n std::auto_ptr<int>, <not created or not owner>\n std::auto_ptr<int>, 0x53b8c0, 1\n"), oss().str());
}
