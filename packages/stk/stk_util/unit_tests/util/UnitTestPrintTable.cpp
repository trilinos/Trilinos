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

#include <stk_util/diag/PrintTable.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST( UnitTestPrintTable, UnitTest)
{
  std::ostringstream oss;

  {
    oss.str("");
    oss << std::endl;

    stk::PrintTable table;
    table.setAutoEndCol(false);

    table.setTitle("Test");

    table << "x" << stk::end_col
          << "y" << stk::end_col << stk::end_header
          << "2" << stk::end_row
          << "" << stk::end_col
          << "3" << stk::end_row
          << "" << stk::end_col
          << "" << stk::end_col
          << "4" << stk::end_row;
    oss << table;
  }
  STKUNIT_ASSERT_EQUAL((std::string("\n Test\nx y\n- - -\n2\n  3\n    4\n\n") == oss.str()), true);

  {
    oss.str("");    
    oss << std::endl;

    stk::PrintTable table;

    table.setTitle("A multiplication table, auto end column");

    table << "x" << "|";
    for (int i = 0; i < 10; ++i)
      table << i << "|";
    table << stk::end_row;
    for (int i = 0; i < 10; ++i) {
      table << i << "|";
      for (int j = 0; j < 10; ++j)
	table << i*j << "|";
      table << stk::end_row;
    }

    oss << table;
  }
  STKUNIT_ASSERT_EQUAL((std::string("\n      A multiplication table, auto end column\nx | 0 | 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |\n0 | 0 | 0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |\n1 | 0 | 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |\n2 | 0 | 2 |  4 |  6 |  8 | 10 | 12 | 14 | 16 | 18 |\n3 | 0 | 3 |  6 |  9 | 12 | 15 | 18 | 21 | 24 | 27 |\n4 | 0 | 4 |  8 | 12 | 16 | 20 | 24 | 28 | 32 | 36 |\n5 | 0 | 5 | 10 | 15 | 20 | 25 | 30 | 35 | 40 | 45 |\n6 | 0 | 6 | 12 | 18 | 24 | 30 | 36 | 42 | 48 | 54 |\n7 | 0 | 7 | 14 | 21 | 28 | 35 | 42 | 49 | 56 | 63 |\n8 | 0 | 8 | 16 | 24 | 32 | 40 | 48 | 56 | 64 | 72 |\n9 | 0 | 9 | 18 | 27 | 36 | 45 | 54 | 63 | 72 | 81 |\n\n") == oss.str()), true);

  {
    oss.str("");
    oss << std::endl;

    stk::PrintTable table;

    table.setTitle("A multiplication table, auto end column");

    table.setAutoEndCol(true);
    table.setCommaSeparatedValues(false);

    table << "x" << "|";
    for (int i = 0; i < 10; ++i)
      table << i << "|";
    table << stk::end_row;
    for (int i = 0; i < 10; ++i) {
      table << i << "|";
      for (int j = 0; j < 10; ++j)
	table << i*j << "|";
      table << stk::end_row;
    }

    oss << table;
  }
  STKUNIT_ASSERT_EQUAL((std::string("\n      A multiplication table, auto end column\nx | 0 | 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |\n0 | 0 | 0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |\n1 | 0 | 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |\n2 | 0 | 2 |  4 |  6 |  8 | 10 | 12 | 14 | 16 | 18 |\n3 | 0 | 3 |  6 |  9 | 12 | 15 | 18 | 21 | 24 | 27 |\n4 | 0 | 4 |  8 | 12 | 16 | 20 | 24 | 28 | 32 | 36 |\n5 | 0 | 5 | 10 | 15 | 20 | 25 | 30 | 35 | 40 | 45 |\n6 | 0 | 6 | 12 | 18 | 24 | 30 | 36 | 42 | 48 | 54 |\n7 | 0 | 7 | 14 | 21 | 28 | 35 | 42 | 49 | 56 | 63 |\n8 | 0 | 8 | 16 | 24 | 32 | 40 | 48 | 56 | 64 | 72 |\n9 | 0 | 9 | 18 | 27 | 36 | 45 | 54 | 63 | 72 | 81 |\n\n") == oss.str()), true);

  {
    oss.str("");
    oss << std::endl;

    stk::PrintTable table;

    table.setTitle("A multiplication table (in hex), no auto end column");

    table.setAutoEndCol(false);
    table.setCommentPrefix("# ");

    table << "x" << stk::end_col << "|" << stk::end_col;
    for (int i = 0; i < 10; ++i)
      table << i << stk::end_col << "|" << stk::end_col;
    table << stk::end_header;
    for (int i = 0; i < 10; ++i) {
      table << i << stk::end_col << "|" << stk::end_col;
      for (int j = 0; j < 10; ++j)
	table << std::hex << i*j << stk::end_col << "|" << stk::end_col;
      table << stk::end_row;
    }

    oss << table;
  }
  STKUNIT_ASSERT_EQUAL((std::string("\n# A multiplication table (in hex), no auto end column\n# x | 0 | 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |\n# - - - - - - -- - -- - -- - -- - -- - -- - -- - -- -\n  0 | 0 | 0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |\n  1 | 0 | 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |\n  2 | 0 | 2 |  4 |  6 |  8 |  a |  c |  e | 10 | 12 |\n  3 | 0 | 3 |  6 |  9 |  c |  f | 12 | 15 | 18 | 1b |\n  4 | 0 | 4 |  8 |  c | 10 | 14 | 18 | 1c | 20 | 24 |\n  5 | 0 | 5 |  a |  f | 14 | 19 | 1e | 23 | 28 | 2d |\n  6 | 0 | 6 |  c | 12 | 18 | 1e | 24 | 2a | 30 | 36 |\n  7 | 0 | 7 |  e | 15 | 1c | 23 | 2a | 31 | 38 | 3f |\n  8 | 0 | 8 | 10 | 18 | 20 | 28 | 30 | 38 | 40 | 48 |\n  9 | 0 | 9 | 12 | 1b | 24 | 2d | 36 | 3f | 48 | 51 |\n  \n") == oss.str()), true);

  {
    oss.str("");
    oss << std::endl;

    stk::PrintTable table;
    
    table.setTitle("A multiplication table, comma separated values");

    table.setAutoEndCol(true);
    table.setCommaSeparatedValues(true);

    table << "x";
    for (int i = 0; i < 10; ++i)
      table << i;
    table << stk::end_header;

    for (int i = 0; i < 10; ++i) {
      table << i;
      for (int j = 0; j < 10; ++j)
	table << i*j;
      table << stk::end_row;
    }
    
    oss << table;
  }
  STKUNIT_ASSERT_EQUAL((std::string("\nA multiplication table, comma separated values\nx,0,1,2,3,4,5,6,7,8,9\n0,0,0,0,0,0,0,0,0,0,0\n1,0,1,2,3,4,5,6,7,8,9\n2,0,2,4,6,8,10,12,14,16,18\n3,0,3,6,9,12,15,18,21,24,27\n4,0,4,8,12,16,20,24,28,32,36\n5,0,5,10,15,20,25,30,35,40,45\n6,0,6,12,18,24,30,36,42,48,54\n7,0,7,14,21,28,35,42,49,56,63\n8,0,8,16,24,32,40,48,56,64,72\n9,0,9,18,27,36,45,54,63,72,81\n\n") == oss.str()), true);
}
