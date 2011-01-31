/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <string>
#include <stk_util/util/IndentStreambuf.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

using stk::push;
using stk::pop;

namespace {

void plain(std::ostream &os)
{
  os << "This is an ostream test" << std::endl;
  os << "This is an ostream test" << stk::PUSH << std::endl;
  os << "This is an ostream test" << stk::POP << std::endl;
  os << "This is an ostream test" << std::endl;
}

void deep(std::ostream &log_stream, int depth)
{
  if (depth < 100) {
    log_stream << "now " << depth << " deep" << stk::PUSH << std::endl;
    deep(log_stream, depth + 1);
    log_stream << stk::POP << std::endl;
  }
}

} // namespace <empty>

STKUNIT_UNIT_TEST(UnitTestIndentStreambuf, UnitTest)
{
  {
    std::string result = 
      "--- indent_streambuf ---\n"
      "indented 0 {\n"
      "  indented 1\n"
      "  indented 2 {\n"
      "    indented 2\n"
      "no indentation\n"
      "    indented 1\n"
      "  }\n"
      "  This is an ostream test\n"
      "  This is an ostream test {\n"
      "    This is an ostream test\n"
      "  }\n"
      "  This is an ostream test\n"
      "  indented 1\n"
      "}\n"
      "indented 0\n";

    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf());
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << std::endl;
    log_stream << "--- indent_streambuf ---" << std::endl;

    log_stream << "indented 0" << stk::PUSH << std::endl;
    log_stream << "indented 1" << std::endl;
    log_stream << stk::PUSH << "indented 2" << std::endl;
    log_stream << "indented 2" << std::endl;
    log_stream << stk::LEFT << "no indentation" << std::endl;
    log_stream << "\017indented 1" << std::endl;

    plain(log_stream);

    log_stream << "indented 1\017" << std::endl;
    log_stream << "indented 0" << std::endl;

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }

  {
    std::string result =
      "--- No braces, no blank lines ---\n"
      "indented 0\n"
      "  indented 1\n"
      "  indented 2\n"
      "    indented 2\n"
      "    indented 1\n"
      "  This is an ostream test\n"
      "  This is an ostream test\n"
      "    This is an ostream test\n"
      "  This is an ostream test\n"
      "  indented 1\n"
      "indented 0\n";
    
    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf(), 2, stk::indent_streambuf::NO_BRACES | stk::indent_streambuf::NO_BLANK_LINES);
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << "--- No braces, no blank lines ---" << std::endl;

    log_stream << "indented 0" << stk::PUSH << std::endl;
    log_stream << "indented 1" << std::endl;
    log_stream << stk::PUSH << "indented 2" << std::endl;
    log_stream << "indented 2" << std::endl;
    log_stream << "\017indented 1" << std::endl;

    plain(log_stream);

    log_stream << "indented 1\017" << std::endl;
    log_stream << "indented 0" << std::endl;

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }

  {
    std::string result =
    "--- push, pop manipulators ---\n"
      "indented 0 {\n"
      "  indented 1 {\n"
      "  }\n"
      "}\n";
    
    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf());
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << "--- push, pop manipulators ---" << std::endl;

    log_stream << "indented 0" << push << std::endl;
    log_stream << "indented 1" << push << std::endl;
    log_stream << pop << std::endl;
    log_stream << pop << std::endl;

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }

  {
    std::string result =
      "--- double push, double pop, push pop, pop push ---\n"
      "push push {\n"
      " {\n"
      "    pop push\n"
      "  }\n"
      " {\n"
      "    pop pop\n"
      "  }\n"
      "}\n";
    
    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf());
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << "--- double push, double pop, push pop, pop push ---" << std::endl;

    log_stream << "push push" << push << push << std::endl;
    log_stream << "pop push" << pop << push << std::endl;
    log_stream << "pop pop" << pop << pop << std::endl;

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }
  
  {
    std::string result =
      "--- No braces, blank lines ---\n"
      "indented 0\n"
      "\n"
      "  indented 1\n"
      "\n"
      "  indented 2\n"
      "\n"
      "    indented 2\n"
      "\n"
      "    indented 1\n"
      "\n"
      "  This is an ostream test\n"
      "  This is an ostream test\n"
      "    This is an ostream test\n"
      "  This is an ostream test\n"
      "  indented 1\n"
      "\n"
      "indented 0\n"
      "\n";
    
    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf(), 2, stk::indent_streambuf::BLANK_LINES);
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << "--- No braces, blank lines ---" << std::endl;

    log_stream << "indented 0" << stk::PUSH << std::endl << std::endl;
    log_stream << "indented 1" << std::endl << std::endl;
    log_stream << stk::PUSH << "indented 2" << std::endl << std::endl;
    log_stream << "indented 2" << std::endl << std::endl;
    log_stream << "\017indented 1" << std::endl << std::endl;

    plain(log_stream);

    log_stream << "indented 1\017" << std::endl << std::endl;
    log_stream << "indented 0" << std::endl << std::endl;

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }

  {
    std::string result =
      "--- Braces, blank lines ---\n"
      "indented 0 {\n"
      "\n"
      "  indented 1\n"
      "\n"
      "  indented 2 {\n"
      "\n"
      "    indented 2\n"
      "\n"
      "    indented 1\n"
      "  }\n"
      "\n"
      "  This is an ostream test\n"
      "  This is an ostream test {\n"
      "    This is an ostream test\n"
      "  }\n"
      "  This is an ostream test\n"
      "  indented 1\n"
      "}\n"
      "\n"
      "indented 0\n"
      "\n";
    
    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf(), 2, stk::indent_streambuf::BRACES | stk::indent_streambuf::BLANK_LINES);
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << "--- Braces, blank lines ---" << std::endl;

    log_stream << "indented 0" << stk::PUSH << std::endl << std::endl;
    log_stream << "indented 1" << std::endl << std::endl;
    log_stream << stk::PUSH << "indented 2" << std::endl << std::endl;
    log_stream << "indented 2" << std::endl << std::endl;
    log_stream << "\017indented 1" << std::endl << std::endl;

    plain(log_stream);

    log_stream << "indented 1\017" << std::endl << std::endl;
    log_stream << "indented 0" << std::endl << std::endl;

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }

  {
    std::string result =
      "Depth test\n"
      "now 0 deep {\n"
      "  now 1 deep {\n"
      "    now 2 deep {\n"
      "      now 3 deep {\n"
      "        now 4 deep {\n"
      "          now 5 deep {\n"
      "            now 6 deep {\n"
      "              now 7 deep {\n"
      "                now 8 deep {\n"
      "                  now 9 deep {\n"
      "                    now 10 deep {\n"
      "                      now 11 deep {\n"
      "                        now 12 deep {\n"
      "                          now 13 deep {\n"
      "                            now 14 deep {\n"
      "                              now 15 deep {\n"
      "                                now 16 deep {\n"
      "                                  now 17 deep {\n"
      "                                    now 18 deep {\n"
      "                                      now 19 deep {\n"
      "                                        now 20 deep {\n"
      "                                          now 21 deep {\n"
      "                                            now 22 deep {\n"
      "                                              now 23 deep {\n"
      "                                                now 24 deep {\n"
      "                                                  now 25 deep {\n"
      "                                                    now 26 deep {\n"
      "                                                      now 27 deep {\n"
      "                                                        now 28 deep {\n"
      "                                                          now 29 deep {\n"
      "                                                            now 30 deep {\n"
      "                                                              now 31 deep {\n"
      "                                                                now 32 deep {\n"
      "                                                                  now 33 deep {\n"
      "                                                                    now 34 deep {\n"
      "                                                                      now 35 deep {\n"
      "                                                                        now 36 deep {\n"
      "                                                                          now 37 deep {\n"
      "                                                                            now 38 deep {\n"
      "                                                                              now 39 deep {\n"
      "                                                                                now 40 deep {\n"
      "                                                                                  now 41 deep {\n"
      "                                                                                    now 42 deep {\n"
      "                                                                                      now 43 deep {\n"
      "                                                                                        now 44 deep {\n"
      "                                                                                          now 45 deep {\n"
      "                                                                                            now 46 deep {\n"
      "                                                                                              now 47 deep {\n"
      "                                                                                                now 48 deep {\n"
      "                                                                                                  now 49 deep {\n"
      "                                                                                                    now 50 deep {\n"
      "                                                                                                    now 51 deep {\n"
      "                                                                                                    now 52 deep {\n"
      "                                                                                                    now 53 deep {\n"
      "                                                                                                    now 54 deep {\n"
      "                                                                                                    now 55 deep {\n"
      "                                                                                                    now 56 deep {\n"
      "                                                                                                    now 57 deep {\n"
      "                                                                                                    now 58 deep {\n"
      "                                                                                                    now 59 deep {\n"
      "                                                                                                    now 60 deep {\n"
      "                                                                                                    now 61 deep {\n"
      "                                                                                                    now 62 deep {\n"
      "                                                                                                    now 63 deep {\n"
      "                                                                                                    now 64 deep {\n"
      "                                                                                                    now 65 deep {\n"
      "                                                                                                    now 66 deep {\n"
      "                                                                                                    now 67 deep {\n"
      "                                                                                                    now 68 deep {\n"
      "                                                                                                    now 69 deep {\n"
      "                                                                                                    now 70 deep {\n"
      "                                                                                                    now 71 deep {\n"
      "                                                                                                    now 72 deep {\n"
      "                                                                                                    now 73 deep {\n"
      "                                                                                                    now 74 deep {\n"
      "                                                                                                    now 75 deep {\n"
      "                                                                                                    now 76 deep {\n"
      "                                                                                                    now 77 deep {\n"
      "                                                                                                    now 78 deep {\n"
      "                                                                                                    now 79 deep {\n"
      "                                                                                                    now 80 deep {\n"
      "                                                                                                    now 81 deep {\n"
      "                                                                                                    now 82 deep {\n"
      "                                                                                                    now 83 deep {\n"
      "                                                                                                    now 84 deep {\n"
      "                                                                                                    now 85 deep {\n"
      "                                                                                                    now 86 deep {\n"
      "                                                                                                    now 87 deep {\n"
      "                                                                                                    now 88 deep {\n"
      "                                                                                                    now 89 deep {\n"
      "                                                                                                    now 90 deep {\n"
      "                                                                                                    now 91 deep {\n"
      "                                                                                                    now 92 deep {\n"
      "                                                                                                    now 93 deep {\n"
      "                                                                                                    now 94 deep {\n"
      "                                                                                                    now 95 deep {\n"
      "                                                                                                    now 96 deep {\n"
      "                                                                                                    now 97 deep {\n"
      "                                                                                                    now 98 deep {\n"
      "                                                                                                    now 99 deep {\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                    }\n"
      "                                                                                                  }\n"
      "                                                                                                }\n"
      "                                                                                              }\n"
      "                                                                                            }\n"
      "                                                                                          }\n"
      "                                                                                        }\n"
      "                                                                                      }\n"
      "                                                                                    }\n"
      "                                                                                  }\n"
      "                                                                                }\n"
      "                                                                              }\n"
      "                                                                            }\n"
      "                                                                          }\n"
      "                                                                        }\n"
      "                                                                      }\n"
      "                                                                    }\n"
      "                                                                  }\n"
      "                                                                }\n"
      "                                                              }\n"
      "                                                            }\n"
      "                                                          }\n"
      "                                                        }\n"
      "                                                      }\n"
      "                                                    }\n"
      "                                                  }\n"
      "                                                }\n"
      "                                              }\n"
      "                                            }\n"
      "                                          }\n"
      "                                        }\n"
      "                                      }\n"
      "                                    }\n"
      "                                  }\n"
      "                                }\n"
      "                              }\n"
      "                            }\n"
      "                          }\n"
      "                        }\n"
      "                      }\n"
      "                    }\n"
      "                  }\n"
      "                }\n"
      "              }\n"
      "            }\n"
      "          }\n"
      "        }\n"
      "      }\n"
      "    }\n"
      "  }\n"
      "}\n";
    

    std::ostringstream    dest;
    stk::indent_streambuf dest_indent_streambuf(dest.rdbuf());
    std::ostream log_stream(&dest_indent_streambuf);

    log_stream << "Depth test" << std::endl;
    deep(log_stream, 0);

    STKUNIT_ASSERT_EQUAL(result, dest.str());
  }
}
