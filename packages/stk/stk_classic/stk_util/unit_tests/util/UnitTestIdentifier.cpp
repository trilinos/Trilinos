/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <utility>
#include <map>
#include <iostream>

#include <stk_util/util/Identifier.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST( UnitTestIdentifier, UnitTest)
{
  {
    stk_classic::IdentifierA id1("identifier1");

    STKUNIT_ASSERT(id1 == "identifier1");
    STKUNIT_ASSERT(id1 == "IDENTIFIER1");
    STKUNIT_ASSERT(id1 == std::string("identifier1"));
    STKUNIT_ASSERT(std::string("identifier1") == id1);
    STKUNIT_ASSERT(id1 < "idf");
    STKUNIT_ASSERT(id1 < "IDF");
    STKUNIT_ASSERT(id1 < std::string("idf"));
    STKUNIT_ASSERT(id1 <= "idf");
    STKUNIT_ASSERT(id1 <= "IDF");
    STKUNIT_ASSERT(id1 <= std::string("idf"));
    STKUNIT_ASSERT(id1 > "idd");
    STKUNIT_ASSERT(id1 > "IDD");
    STKUNIT_ASSERT(id1 > std::string("idd"));
    STKUNIT_ASSERT(id1 >= "idd");
    STKUNIT_ASSERT(id1 >= "IDD");
    STKUNIT_ASSERT(id1 >= std::string("idd"));

    STKUNIT_ASSERT(id1 <= "identifier1");
    STKUNIT_ASSERT(id1 <= "IDENTIFIER1");
    STKUNIT_ASSERT(id1 <= std::string("identifier1"));
    STKUNIT_ASSERT(id1 >= "identifier1");
    STKUNIT_ASSERT(id1 >= "IDENTIFIER1");
    STKUNIT_ASSERT(id1 >= std::string("identifier1"));

    stk_classic::IdentifierA id2(id1);

    STKUNIT_ASSERT(id1 == id2);

    std::cout << id1 << std::endl;

    stk_classic::IdentifierA id3 = id1 + "test1";
    STKUNIT_ASSERT(id3 == "identifier1test1");

    id3 += "test2";
    STKUNIT_ASSERT(id3 == "identifier1test1test2");

    id3 = "identifier3";
    STKUNIT_ASSERT(id3 == "identifier3");

    typedef std::map<stk_classic::IdentifierA, int> IdIntMap;

    IdIntMap id_int_map;

    id_int_map[stk_classic::IdentifierA("identifier1")] = 1;
    id_int_map[stk_classic::IdentifierA("IDENTIFIER1")] = 2;

    STKUNIT_ASSERT(id_int_map[stk_classic::IdentifierA("identifier1")] == 2);
  }


  {
    stk_classic::IdentifierB id1("identifier1");


    STKUNIT_ASSERT(id1 == "identifier1");
    STKUNIT_ASSERT(id1 == "IDENTIFIER1");
    STKUNIT_ASSERT(id1 == std::string("identifier1"));
    STKUNIT_ASSERT(std::string("identifier1") == id1);
    STKUNIT_ASSERT(id1 < "idf");
    STKUNIT_ASSERT(id1 < "IDF");
    STKUNIT_ASSERT(id1 < std::string("idf"));
    STKUNIT_ASSERT(id1 <= "idf");
    STKUNIT_ASSERT(id1 <= "IDF");
    STKUNIT_ASSERT(id1 <= std::string("idf"));
    STKUNIT_ASSERT(id1 > "idd");
    STKUNIT_ASSERT(id1 > "IDD");
    STKUNIT_ASSERT(id1 > std::string("idd"));
    STKUNIT_ASSERT(id1 >= "idd");
    STKUNIT_ASSERT(id1 >= "IDD");
    STKUNIT_ASSERT(id1 >= std::string("idd"));

    STKUNIT_ASSERT(id1 <= "identifier1");
    STKUNIT_ASSERT(id1 <= "IDENTIFIER1");
    STKUNIT_ASSERT(id1 <= std::string("identifier1"));
    STKUNIT_ASSERT(id1 >= "identifier1");
    STKUNIT_ASSERT(id1 >= "IDENTIFIER1");
    STKUNIT_ASSERT(id1 >= std::string("identifier1"));

    stk_classic::IdentifierB id2(id1);

    STKUNIT_ASSERT(id1 == id2);

    std::cout << id1 << std::endl;

    stk_classic::IdentifierB id3 = id1 + "test1";
    STKUNIT_ASSERT(id3 == "identifier1test1");

    id3 += "test2";
    STKUNIT_ASSERT(id3 == "identifier1test1test2");

    id3 = "identifier3";
    STKUNIT_ASSERT(id3 == "identifier3");

    typedef std::map<stk_classic::IdentifierB, int> IdIntMap;

    IdIntMap id_int_map;

    id_int_map[stk_classic::IdentifierB("identifier1")] = 1;
    id_int_map[stk_classic::IdentifierB("IDENTIFIER1")] = 2;

    STKUNIT_ASSERT(id_int_map[stk_classic::IdentifierB("identifier1")] == 2);
  }
}

