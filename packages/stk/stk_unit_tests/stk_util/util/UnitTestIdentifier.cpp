/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>                     // for ostream, cout, endl
#include <map>                          // for map, map<>::value_compare
#include <gtest/gtest.h>
#include <stk_util/util/Identifier.hpp>  // for IdentifierA, IdentifierB, etc
#include <string>                       // for string, basic_string



TEST( UnitTestIdentifier, UnitTest)
{
  {
    stk::IdentifierA id1("identifier1");

    ASSERT_TRUE(id1 == "identifier1");
    ASSERT_TRUE(id1 == "IDENTIFIER1");
    ASSERT_TRUE(id1 == std::string("identifier1"));
    ASSERT_TRUE(std::string("identifier1") == id1);
    ASSERT_TRUE(id1 < "idf");
    ASSERT_TRUE(id1 < "IDF");
    ASSERT_TRUE(id1 < std::string("idf"));
    ASSERT_TRUE(id1 <= "idf");
    ASSERT_TRUE(id1 <= "IDF");
    ASSERT_TRUE(id1 <= std::string("idf"));
    ASSERT_TRUE(id1 > "idd");
    ASSERT_TRUE(id1 > "IDD");
    ASSERT_TRUE(id1 > std::string("idd"));
    ASSERT_TRUE(id1 >= "idd");
    ASSERT_TRUE(id1 >= "IDD");
    ASSERT_TRUE(id1 >= std::string("idd"));

    ASSERT_TRUE(id1 <= "identifier1");
    ASSERT_TRUE(id1 <= "IDENTIFIER1");
    ASSERT_TRUE(id1 <= std::string("identifier1"));
    ASSERT_TRUE(id1 >= "identifier1");
    ASSERT_TRUE(id1 >= "IDENTIFIER1");
    ASSERT_TRUE(id1 >= std::string("identifier1"));

    stk::IdentifierA id2(id1);

    ASSERT_TRUE(id1 == id2);

    std::cout << id1 << std::endl;

    stk::IdentifierA id3 = id1 + "test1";
    ASSERT_TRUE(id3 == "identifier1test1");

    id3 += "test2";
    ASSERT_TRUE(id3 == "identifier1test1test2");

    id3 = "identifier3";
    ASSERT_TRUE(id3 == "identifier3");

    typedef std::map<stk::IdentifierA, int> IdIntMap;

    IdIntMap id_int_map;

    id_int_map[stk::IdentifierA("identifier1")] = 1;
    id_int_map[stk::IdentifierA("IDENTIFIER1")] = 2;

    ASSERT_TRUE(id_int_map[stk::IdentifierA("identifier1")] == 2);
  }


  {
    stk::IdentifierB id1("identifier1");


    ASSERT_TRUE(id1 == "identifier1");
    ASSERT_TRUE(id1 == "IDENTIFIER1");
    ASSERT_TRUE(id1 == std::string("identifier1"));
    ASSERT_TRUE(std::string("identifier1") == id1);
    ASSERT_TRUE(id1 < "idf");
    ASSERT_TRUE(id1 < "IDF");
    ASSERT_TRUE(id1 < std::string("idf"));
    ASSERT_TRUE(id1 <= "idf");
    ASSERT_TRUE(id1 <= "IDF");
    ASSERT_TRUE(id1 <= std::string("idf"));
    ASSERT_TRUE(id1 > "idd");
    ASSERT_TRUE(id1 > "IDD");
    ASSERT_TRUE(id1 > std::string("idd"));
    ASSERT_TRUE(id1 >= "idd");
    ASSERT_TRUE(id1 >= "IDD");
    ASSERT_TRUE(id1 >= std::string("idd"));

    ASSERT_TRUE(id1 <= "identifier1");
    ASSERT_TRUE(id1 <= "IDENTIFIER1");
    ASSERT_TRUE(id1 <= std::string("identifier1"));
    ASSERT_TRUE(id1 >= "identifier1");
    ASSERT_TRUE(id1 >= "IDENTIFIER1");
    ASSERT_TRUE(id1 >= std::string("identifier1"));

    stk::IdentifierB id2(id1);

    ASSERT_TRUE(id1 == id2);

    std::cout << id1 << std::endl;

    stk::IdentifierB id3 = id1 + "test1";
    ASSERT_TRUE(id3 == "identifier1test1");

    id3 += "test2";
    ASSERT_TRUE(id3 == "identifier1test1test2");

    id3 = "identifier3";
    ASSERT_TRUE(id3 == "identifier3");

    typedef std::map<stk::IdentifierB, int> IdIntMap;

    IdIntMap id_int_map;

    id_int_map[stk::IdentifierB("identifier1")] = 1;
    id_int_map[stk::IdentifierB("IDENTIFIER1")] = 2;

    ASSERT_TRUE(id_int_map[stk::IdentifierB("identifier1")] == 2);
  }
}

