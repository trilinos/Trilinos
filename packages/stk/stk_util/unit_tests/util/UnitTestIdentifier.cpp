#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <utility>
#include <map>
#include <iostream>

#include <stk_util/util/Identifier.hpp>

class UnitTestIdentifier : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestIdentifier );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION( UnitTestIdentifier);

void UnitTestIdentifier::testUnit()
{
  {
    stk::IdentifierA id1("identifier1");

    CPPUNIT_ASSERT(id1 == "identifier1");
    CPPUNIT_ASSERT(id1 == "IDENTIFIER1");
    CPPUNIT_ASSERT(id1 == std::string("identifier1"));
    CPPUNIT_ASSERT(std::string("identifier1") == id1);  
    CPPUNIT_ASSERT(id1 < "idf");
    CPPUNIT_ASSERT(id1 < "IDF");
    CPPUNIT_ASSERT(id1 < std::string("idf"));
    CPPUNIT_ASSERT(id1 <= "idf");
    CPPUNIT_ASSERT(id1 <= "IDF");
    CPPUNIT_ASSERT(id1 <= std::string("idf"));
    CPPUNIT_ASSERT(id1 > "idd");
    CPPUNIT_ASSERT(id1 > "IDD");
    CPPUNIT_ASSERT(id1 > std::string("idd"));
    CPPUNIT_ASSERT(id1 >= "idd");
    CPPUNIT_ASSERT(id1 >= "IDD");
    CPPUNIT_ASSERT(id1 >= std::string("idd"));

    CPPUNIT_ASSERT(id1 <= "identifier1");
    CPPUNIT_ASSERT(id1 <= "IDENTIFIER1");
    CPPUNIT_ASSERT(id1 <= std::string("identifier1"));
    CPPUNIT_ASSERT(id1 >= "identifier1");
    CPPUNIT_ASSERT(id1 >= "IDENTIFIER1");
    CPPUNIT_ASSERT(id1 >= std::string("identifier1"));
  
    stk::IdentifierA id2(id1);

    CPPUNIT_ASSERT(id1 == id2);

    std::cout << id1 << std::endl;

    stk::IdentifierA id3 = id1 + "test1";
    CPPUNIT_ASSERT(id3 == "identifier1test1");

    id3 += "test2";
    CPPUNIT_ASSERT(id3 == "identifier1test1test2");

    id3 = "identifier3";
    CPPUNIT_ASSERT(id3 == "identifier3");

    typedef std::map<stk::IdentifierA, int> IdIntMap;

    IdIntMap id_int_map;

    id_int_map[stk::IdentifierA("identifier1")] = 1;
    id_int_map[stk::IdentifierA("IDENTIFIER1")] = 2;

    CPPUNIT_ASSERT(id_int_map[stk::IdentifierA("identifier1")] == 2);    
  }


  {
    stk::IdentifierB id1("identifier1");


    CPPUNIT_ASSERT(id1 == "identifier1");
    CPPUNIT_ASSERT(id1 == "IDENTIFIER1");
    CPPUNIT_ASSERT(id1 == std::string("identifier1"));
    CPPUNIT_ASSERT(std::string("identifier1") == id1);  
    CPPUNIT_ASSERT(id1 < "idf");
    CPPUNIT_ASSERT(id1 < "IDF");
    CPPUNIT_ASSERT(id1 < std::string("idf"));
    CPPUNIT_ASSERT(id1 <= "idf");
    CPPUNIT_ASSERT(id1 <= "IDF");
    CPPUNIT_ASSERT(id1 <= std::string("idf"));
    CPPUNIT_ASSERT(id1 > "idd");
    CPPUNIT_ASSERT(id1 > "IDD");
    CPPUNIT_ASSERT(id1 > std::string("idd"));
    CPPUNIT_ASSERT(id1 >= "idd");
    CPPUNIT_ASSERT(id1 >= "IDD");
    CPPUNIT_ASSERT(id1 >= std::string("idd"));

    CPPUNIT_ASSERT(id1 <= "identifier1");
    CPPUNIT_ASSERT(id1 <= "IDENTIFIER1");
    CPPUNIT_ASSERT(id1 <= std::string("identifier1"));
    CPPUNIT_ASSERT(id1 >= "identifier1");
    CPPUNIT_ASSERT(id1 >= "IDENTIFIER1");
    CPPUNIT_ASSERT(id1 >= std::string("identifier1"));
  
    stk::IdentifierB id2(id1);

    CPPUNIT_ASSERT(id1 == id2);

    std::cout << id1 << std::endl;

    stk::IdentifierB id3 = id1 + "test1";
    CPPUNIT_ASSERT(id3 == "identifier1test1");

    id3 += "test2";
    CPPUNIT_ASSERT(id3 == "identifier1test1test2");

    id3 = "identifier3";
    CPPUNIT_ASSERT(id3 == "identifier3");

    typedef std::map<stk::IdentifierB, int> IdIntMap;

    IdIntMap id_int_map;

    id_int_map[stk::IdentifierB("identifier1")] = 1;
    id_int_map[stk::IdentifierB("IDENTIFIER1")] = 2;

    CPPUNIT_ASSERT(id_int_map[stk::IdentifierB("identifier1")] == 2);    
  }
}

