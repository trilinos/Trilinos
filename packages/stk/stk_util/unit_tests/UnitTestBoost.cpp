#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <utility>
#include <iostream>
#include <cstring>
#include <string>

//boost tr1 headers...
//On the sun, couldn't get '#include <memory>' to work, so we're using the boost
//form instead...
#include <boost/tr1/memory.hpp>
#include <boost/tr1/array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/shared_array.hpp>

#include <stk_util/util/ci_string.hpp>

#include <boost/program_options.hpp>

namespace {
  char * my_strdup(const char *s) {
    return std::strcpy(new char[std::strlen(s) + 1], s);
  }
}

class UnitTestBoost : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestBoost );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();
  char *m_argv[2];
public:
  UnitTestBoost() : CppUnit::TestCase() {
    m_argv[0] = my_strdup("UnitTestBoost");
    m_argv[1] = my_strdup("--compression=1");
  }
  ~UnitTestBoost() {
    delete [] m_argv[0];
    delete [] m_argv[1];
  }
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION( UnitTestBoost);



namespace boost {

template <>
struct hash<ci_string>
{
  std::size_t operator()(const ci_string &s) const {
    std::size_t seed = 0;
    
    for(ci_string::const_iterator first = s.begin(); first != s.end(); ++first) {
      boost::hash<char> hasher;
      seed ^= hasher(std::tolower(*first)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    
    return seed;
  }
};
  
} // namespace boost


void UnitTestBoost::testUnit()
{
  {  
    double* d = new double;
    boost::shared_ptr<double> dptr(d);

    CPPUNIT_ASSERT_EQUAL( dptr.get(), d);

    double* d2 = new double[1];
    boost::shared_array<double> dptr2(d2);

    CPPUNIT_ASSERT_EQUAL( dptr2.get(), d2);
  }
  
  boost::array<double,5> my_array;

  my_array[0] = 5.0;

  CPPUNIT_ASSERT_EQUAL( my_array[0], 5.0 );
  CPPUNIT_ASSERT_EQUAL( my_array.size(), (boost::array<double,5>::size_type)5 );

  boost::unordered_set<int> int_set;

  int_set.insert(5);

  CPPUNIT_ASSERT_EQUAL( int_set.size(), (boost::unordered_set<int>::size_type)1 );

  boost::unordered_set<ci_string> ci_string_set;

  ci_string_set.insert("Test");
  std::pair<boost::unordered_set<ci_string>::iterator, bool> res = ci_string_set.insert("test");

  CPPUNIT_ASSERT_EQUAL( ci_string_set.size(), (boost::unordered_set<ci_string>::size_type)1 );
  CPPUNIT_ASSERT_EQUAL( res.second, false );

  ci_string s("This is a test");

  CPPUNIT_ASSERT( s == "this is a test" );
  
  std::cout << s << std::endl;
  
  namespace po = boost::program_options;

// Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("compression", po::value<int>(), "set compression level")
    ;

  int argc = sizeof(m_argv)/sizeof(m_argv[0]);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, m_argv, desc), vm);
  po::notify(vm);    

  CPPUNIT_ASSERT_EQUAL(vm["compression"].as<int>(), 1);
}

