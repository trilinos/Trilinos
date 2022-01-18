// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "gtest/gtest.h"
#include "stk_util/util/string_utils.hpp"  // for make_vector_of_strings
#include <string>                          // for string
#include <vector>                          // for vector


TEST( UnitTestStringUtils, makeVectorOfStrings_1shortString)
{
  std::string str("12345");
  std::vector<std::string> expected = { "12345" };

  std::vector<std::string> vecStrs = stk::make_vector_of_strings(str, ' ', 9);

  EXPECT_EQ(expected, vecStrs);
}

TEST( UnitTestStringUtils, makeVectorOfStrings_2strings1separator)
{
  std::string str("12345 123");
  std::vector<std::string> expected = { "12345", "123" };

  std::vector<std::string> vecStrs = stk::make_vector_of_strings(str, ' ', 7);

  EXPECT_EQ(expected, vecStrs);
}

TEST( UnitTestStringUtils, makeVectorOfStrings_3strings2separators)
{
  std::string str("123456789 123456789 1234567890");
  std::vector<std::string> expected = {
    "123456789", "123456789", "123456789", "0" };

  std::vector<std::string> vecStrs = stk::make_vector_of_strings(str, ' ', 9);

  EXPECT_EQ(expected, vecStrs);
}

TEST( UnitTestStringUtils, makeVectorOfStrings_linebreaks)
{
  std::string str("123456789\n\n123\n123\n\n123456789 123456789");
  std::vector<std::string> expected = {
    "123456789", "", "123", "123", "", "123456789", "123456789" };

  std::vector<std::string> vecStrs = stk::make_vector_of_strings(str, ' ', 9);

  EXPECT_EQ(expected, vecStrs);
}

TEST(StringStartsWith, empty) {
  EXPECT_TRUE(stk::string_starts_with("", ""));
  EXPECT_TRUE(stk::string_starts_with("test", ""));
}

TEST(StringStartsWith, notFound) {
  EXPECT_FALSE(stk::string_starts_with("test", "x"));
  EXPECT_FALSE(stk::string_starts_with("test", "testx"));
  EXPECT_FALSE(stk::string_starts_with("test", "TEST"));
  EXPECT_FALSE(stk::string_starts_with("test", " test"));
  EXPECT_FALSE(stk::string_starts_with(" test", "test"));
}

TEST(StringStartsWith, found) {
  EXPECT_TRUE(stk::string_starts_with("test", "test"));
  EXPECT_TRUE(stk::string_starts_with("testx", "test"));
  EXPECT_TRUE(stk::string_starts_with("Test", "Test"));
  EXPECT_TRUE(stk::string_starts_with(" test", " test"));
  EXPECT_TRUE(stk::string_starts_with(" testx", " test"));
  EXPECT_TRUE(stk::string_starts_with(" ", " "));
  EXPECT_TRUE(stk::string_starts_with("   ", " "));
  EXPECT_TRUE(stk::string_starts_with("\ntest", "\n"));
}

TEST(LtrimString, empty)
{
  EXPECT_EQ(stk::ltrim_string(""), "");
}

TEST(LtrimString, allWhitespace)
{
  EXPECT_EQ(stk::ltrim_string(       " "), "");
  EXPECT_EQ(stk::ltrim_string(      "  "), "");
  EXPECT_EQ(stk::ltrim_string("        "), "");
  EXPECT_EQ(stk::ltrim_string(      "\t"), "");
  EXPECT_EQ(stk::ltrim_string(    "\t  "), "");
  EXPECT_EQ(stk::ltrim_string(    "  \t"), "");
  EXPECT_EQ(stk::ltrim_string(    " \t "), "");
  EXPECT_EQ(stk::ltrim_string(      "\n"), "");
  EXPECT_EQ(stk::ltrim_string(    "\n\t"), "");
  EXPECT_EQ(stk::ltrim_string(      "\r"), "");
  EXPECT_EQ(stk::ltrim_string(    "\r\n"), "");
  EXPECT_EQ(stk::ltrim_string(  "\r\n\t"), "");
  EXPECT_EQ(stk::ltrim_string( "\r\n\t "), "");
}

TEST(LtrimString, oneWordNoWhitespace)
{
  EXPECT_EQ(stk::ltrim_string("word"), "word");
}

TEST(LtrimString, oneWordTrailingWhitespace)
{
  EXPECT_EQ(stk::ltrim_string(       "word "), "word "       );
  EXPECT_EQ(stk::ltrim_string(      "word  "), "word  "      );
  EXPECT_EQ(stk::ltrim_string("word        "), "word        ");
  EXPECT_EQ(stk::ltrim_string(      "word\t"), "word\t"      );
  EXPECT_EQ(stk::ltrim_string(    "word\t  "), "word\t  "    );
  EXPECT_EQ(stk::ltrim_string(    "word  \t"), "word  \t"    );
  EXPECT_EQ(stk::ltrim_string(    "word \t "), "word \t "    );
  EXPECT_EQ(stk::ltrim_string(      "word\n"), "word\n"      );
  EXPECT_EQ(stk::ltrim_string(    "word\n  "), "word\n  "    );
  EXPECT_EQ(stk::ltrim_string(    "word  \n"), "word  \n"    );
  EXPECT_EQ(stk::ltrim_string(    "word \n "), "word \n "    );
  EXPECT_EQ(stk::ltrim_string(      "word\r"), "word\r"      );
  EXPECT_EQ(stk::ltrim_string(    "word\r  "), "word\r  "    );
  EXPECT_EQ(stk::ltrim_string(    "word  \r"), "word  \r"    );
  EXPECT_EQ(stk::ltrim_string(    "word \r "), "word \r "    );
}

TEST(LtrimString, oneWordLeadingWhitespace)
{
  EXPECT_EQ(stk::ltrim_string(       " word"), "word");
  EXPECT_EQ(stk::ltrim_string(      "  word"), "word");
  EXPECT_EQ(stk::ltrim_string("        word"), "word");
  EXPECT_EQ(stk::ltrim_string(      "\tword"), "word");
  EXPECT_EQ(stk::ltrim_string(    "\t  word"), "word");
  EXPECT_EQ(stk::ltrim_string(    "  \tword"), "word");
  EXPECT_EQ(stk::ltrim_string(    " \t word"), "word");
  EXPECT_EQ(stk::ltrim_string(    "\t\nword"), "word");
  EXPECT_EQ(stk::ltrim_string(      "\nword"), "word");
  EXPECT_EQ(stk::ltrim_string(    "\n  word"), "word");
  EXPECT_EQ(stk::ltrim_string(    "  \nword"), "word");
  EXPECT_EQ(stk::ltrim_string(    " \n word"), "word");
  EXPECT_EQ(stk::ltrim_string(      "\rword"), "word");
  EXPECT_EQ(stk::ltrim_string(    "\r  word"), "word");
  EXPECT_EQ(stk::ltrim_string(    "  \rword"), "word");
  EXPECT_EQ(stk::ltrim_string(    " \r word"), "word");
}

TEST(LtrimString, oneWordMixedWhitespace)
{
  EXPECT_EQ(stk::ltrim_string(      " word "), "word "   );
  EXPECT_EQ(stk::ltrim_string(    "  word  "), "word  "  );
  EXPECT_EQ(stk::ltrim_string(    "\tword\t"), "word\t"  );
  EXPECT_EQ(stk::ltrim_string("\t  word\t  "), "word\t  ");
  EXPECT_EQ(stk::ltrim_string("  \tword  \t"), "word  \t");
  EXPECT_EQ(stk::ltrim_string(" \t word \t "), "word \t ");
  EXPECT_EQ(stk::ltrim_string(    "\nword\n"), "word\n"  );
  EXPECT_EQ(stk::ltrim_string("\n  word\n  "), "word\n  ");
  EXPECT_EQ(stk::ltrim_string("  \nword  \n"), "word  \n");
  EXPECT_EQ(stk::ltrim_string(" \n word \n "), "word \n ");
  EXPECT_EQ(stk::ltrim_string(    "\rword\r"), "word\r"  );
  EXPECT_EQ(stk::ltrim_string("\r  word\r  "), "word\r  ");
  EXPECT_EQ(stk::ltrim_string("  \rword  \r"), "word  \r");
  EXPECT_EQ(stk::ltrim_string(" \r word \r "), "word \r ");
}

TEST(LtrimString, twoWordsMixedWhitespace)
{
  EXPECT_EQ(stk::ltrim_string(         "word word"), "word word"     );
  EXPECT_EQ(stk::ltrim_string(       " word word "), "word word "    );
  EXPECT_EQ(stk::ltrim_string(    "\tword\tword\t"), "word\tword\t"  );
  EXPECT_EQ(stk::ltrim_string( " \t word word \t "), "word word \t " );
  EXPECT_EQ(stk::ltrim_string(" \t word\tword \t "), "word\tword \t ");
  EXPECT_EQ(stk::ltrim_string( " \n word word \n "), "word word \n " );
  EXPECT_EQ(stk::ltrim_string(" \n word\nword \n "), "word\nword \n ");
  EXPECT_EQ(stk::ltrim_string( " \r word word \r "), "word word \r " );
  EXPECT_EQ(stk::ltrim_string(" \r word\nword \r "), "word\nword \r ");
}

TEST(RtrimString, empty)
{
  EXPECT_EQ(stk::rtrim_string(""), "");
}

TEST(RtrimString, allWhitespace)
{
  EXPECT_EQ(stk::rtrim_string(       " "), "");
  EXPECT_EQ(stk::rtrim_string(      "  "), "");
  EXPECT_EQ(stk::rtrim_string("        "), "");
  EXPECT_EQ(stk::rtrim_string(      "\t"), "");
  EXPECT_EQ(stk::rtrim_string(    "\t  "), "");
  EXPECT_EQ(stk::rtrim_string(    "  \t"), "");
  EXPECT_EQ(stk::rtrim_string(    " \t "), "");
  EXPECT_EQ(stk::rtrim_string(      "\n"), "");
  EXPECT_EQ(stk::rtrim_string(    "\n\t"), "");
  EXPECT_EQ(stk::rtrim_string(      "\r"), "");
  EXPECT_EQ(stk::rtrim_string(    "\r\n"), "");
  EXPECT_EQ(stk::rtrim_string(  "\r\n\t"), "");
  EXPECT_EQ(stk::rtrim_string( "\r\n\t "), "");
}

TEST(RtrimString, oneWordNoWhitespace)
{
  EXPECT_EQ(stk::rtrim_string("word"), "word");
}

TEST(RtrimString, oneWordLeadingWhitespace)
{
  EXPECT_EQ(stk::rtrim_string(       " word"),        " word");
  EXPECT_EQ(stk::rtrim_string(      "  word"),       "  word");
  EXPECT_EQ(stk::rtrim_string("        word"), "        word");
  EXPECT_EQ(stk::rtrim_string(      "\tword"),       "\tword");
  EXPECT_EQ(stk::rtrim_string(    "\t  word"),     "\t  word");
  EXPECT_EQ(stk::rtrim_string(    "  \tword"),     "  \tword");
  EXPECT_EQ(stk::rtrim_string(    " \t word"),     " \t word");
  EXPECT_EQ(stk::rtrim_string(   " \t\nword"),    " \t\nword");
  EXPECT_EQ(stk::rtrim_string(      "\nword"),       "\nword");
  EXPECT_EQ(stk::rtrim_string(    "\n  word"),     "\n  word");
  EXPECT_EQ(stk::rtrim_string(    "  \nword"),     "  \nword");
  EXPECT_EQ(stk::rtrim_string(    " \n word"),     " \n word");
  EXPECT_EQ(stk::rtrim_string(      "\rword"),       "\rword");
  EXPECT_EQ(stk::rtrim_string(    "\r  word"),     "\r  word");
  EXPECT_EQ(stk::rtrim_string(    "  \rword"),     "  \rword");
  EXPECT_EQ(stk::rtrim_string(    " \r word"),     " \r word");
}

TEST(RtrimString, oneWordTrailingWhitespace)
{
  EXPECT_EQ(stk::rtrim_string(       "word "), "word");
  EXPECT_EQ(stk::rtrim_string(      "word  "), "word");
  EXPECT_EQ(stk::rtrim_string("word        "), "word");
  EXPECT_EQ(stk::rtrim_string(      "word\t"), "word");
  EXPECT_EQ(stk::rtrim_string(    "word\t  "), "word");
  EXPECT_EQ(stk::rtrim_string(    "word  \t"), "word");
  EXPECT_EQ(stk::rtrim_string(    "word \t "), "word");
  EXPECT_EQ(stk::rtrim_string(      "word\n"), "word");
  EXPECT_EQ(stk::rtrim_string(    "word\n  "), "word");
  EXPECT_EQ(stk::rtrim_string(    "word  \n"), "word");
  EXPECT_EQ(stk::rtrim_string(    "word \n "), "word");
  EXPECT_EQ(stk::rtrim_string(      "word\r"), "word");
  EXPECT_EQ(stk::rtrim_string(    "word\r  "), "word");
  EXPECT_EQ(stk::rtrim_string(    "word  \r"), "word");
  EXPECT_EQ(stk::rtrim_string(    "word \r "), "word");
}

TEST(RtrimString, oneWordMixedWhitespace)
{
  EXPECT_EQ(stk::rtrim_string(      " word "),    " word");
  EXPECT_EQ(stk::rtrim_string(    "  word  "),   "  word");
  EXPECT_EQ(stk::rtrim_string(    "\tword\t"),   "\tword");
  EXPECT_EQ(stk::rtrim_string("\t  word\t  "), "\t  word");
  EXPECT_EQ(stk::rtrim_string("  \tword  \t"), "  \tword");
  EXPECT_EQ(stk::rtrim_string(" \t word \t "), " \t word");
  EXPECT_EQ(stk::rtrim_string(    "\nword\n"),   "\nword");
  EXPECT_EQ(stk::rtrim_string("\n  word\n  "), "\n  word");
  EXPECT_EQ(stk::rtrim_string("  \nword  \n"), "  \nword");
  EXPECT_EQ(stk::rtrim_string(" \n word \n "), " \n word");
  EXPECT_EQ(stk::rtrim_string(    "\rword\r"),   "\rword");
  EXPECT_EQ(stk::rtrim_string("\r  word\r  "), "\r  word");
  EXPECT_EQ(stk::rtrim_string("  \rword  \r"), "  \rword");
  EXPECT_EQ(stk::rtrim_string(" \r word \r "), " \r word");
}

TEST(RtrimString, twoWordsMixedWhitespace)
{
  EXPECT_EQ(stk::rtrim_string(         "word word"),      "word word");
  EXPECT_EQ(stk::rtrim_string(       " word word "),     " word word");
  EXPECT_EQ(stk::rtrim_string(    "\tword\tword\t"),   "\tword\tword");
  EXPECT_EQ(stk::rtrim_string( " \t word word \t "),  " \t word word");
  EXPECT_EQ(stk::rtrim_string(" \t word\tword \t "), " \t word\tword");
  EXPECT_EQ(stk::rtrim_string( " \n word word \n "),  " \n word word");
  EXPECT_EQ(stk::rtrim_string(" \n word\nword \n "), " \n word\nword");
  EXPECT_EQ(stk::rtrim_string( " \r word word \r "),  " \r word word");
  EXPECT_EQ(stk::rtrim_string(" \r word\rword \r "), " \r word\rword");
}

TEST(TrimString, oneWordNoWhitespace)
{
  EXPECT_EQ(stk::trim_string("word"), "word");
}

TEST(TrimString, oneWordMixedWhitespace)
{
  EXPECT_EQ(stk::trim_string(      " word "), "word");
  EXPECT_EQ(stk::trim_string(    "  word  "), "word");
  EXPECT_EQ(stk::trim_string(    "\tword\t"), "word");
  EXPECT_EQ(stk::trim_string("\t  word\t  "), "word");
  EXPECT_EQ(stk::trim_string("  \tword  \t"), "word");
  EXPECT_EQ(stk::trim_string(" \t word \t "), "word");
  EXPECT_EQ(stk::trim_string(    "\nword\n"), "word");
  EXPECT_EQ(stk::trim_string("\n  word\n  "), "word");
  EXPECT_EQ(stk::trim_string("  \nword  \n"), "word");
  EXPECT_EQ(stk::trim_string(" \n word \n "), "word");
  EXPECT_EQ(stk::trim_string(    "\rword\r"), "word");
  EXPECT_EQ(stk::trim_string("\r  word\r  "), "word");
  EXPECT_EQ(stk::trim_string("  \rword  \r"), "word");
  EXPECT_EQ(stk::trim_string(" \r word \r "), "word");
}

TEST(TrimString, twoWordsMixedWhitespace)
{
  EXPECT_EQ(stk::trim_string(         "word word"), "word word");
  EXPECT_EQ(stk::trim_string(       " word word "), "word word");
  EXPECT_EQ(stk::trim_string(    "\tword\tword\t"), "word\tword");
  EXPECT_EQ(stk::trim_string( " \t word word \t "), "word word");
  EXPECT_EQ(stk::trim_string(" \t word\tword \t "), "word\tword");
  EXPECT_EQ(stk::trim_string( " \n word word \n "), "word word");
  EXPECT_EQ(stk::trim_string(" \n word\nword \n "), "word\nword");
  EXPECT_EQ(stk::trim_string( " \r word word \r "), "word word");
  EXPECT_EQ(stk::trim_string(" \r word\rword \r "), "word\rword");
}

TEST(SplitString, empty)
{
  EXPECT_EQ(stk::split_string("", ':'), std::vector<std::string>{});
}

TEST(SplitString, whitespace)
{
  EXPECT_EQ(stk::split_string("   ", ':'), std::vector<std::string>{"   "});
}

TEST(SplitString, emptyFields)
{
  EXPECT_EQ(stk::split_string(":", ':'),    (std::vector<std::string>{"", ""}            ));
  EXPECT_EQ(stk::split_string("::", ':'),   (std::vector<std::string>{"", "", ""}        ));
  EXPECT_EQ(stk::split_string(":::", ':'),  (std::vector<std::string>{"", "", "", ""}    ));
  EXPECT_EQ(stk::split_string("::::", ':'), (std::vector<std::string>{"", "", "", "", ""}));
}

TEST(SplitString, multipleNonemptyFields)
{
  EXPECT_EQ(stk::split_string("a", ':'),           (std::vector<std::string>{"a"}                ));
  EXPECT_EQ(stk::split_string("a:b", ':'),         (std::vector<std::string>{"a", "b"}           ));
  EXPECT_EQ(stk::split_string("a:b:c", ':'),       (std::vector<std::string>{"a", "b", "c"}      ));
  EXPECT_EQ(stk::split_string("a:b:c:d", ':'),     (std::vector<std::string>{"a", "b", "c", "d"} ));
  EXPECT_EQ(stk::split_string("a,b:c=d:e f", ':'), (std::vector<std::string>{"a,b", "c=d", "e f"}));
}

TEST(SplitString, mixedEmptyNonemptyFields)
{
  EXPECT_EQ(stk::split_string("a:", ':'),   (std::vector<std::string>{"a", ""}     ));
  EXPECT_EQ(stk::split_string(":a", ':'),   (std::vector<std::string>{"", "a"}     ));
  EXPECT_EQ(stk::split_string("a::b", ':'), (std::vector<std::string>{"a", "", "b"}));
  EXPECT_EQ(stk::split_string("a::", ':'),  (std::vector<std::string>{"a", "", ""} ));
  EXPECT_EQ(stk::split_string("::a", ':'),  (std::vector<std::string>{"", "", "a" }));
}

TEST(SplitString, mixedEmptyNonemptyWhitespaceFields)
{
  EXPECT_EQ(stk::split_string(" a : ", ':'),     (std::vector<std::string>{" a ", " "}       ));
  EXPECT_EQ(stk::split_string(" : a\n", ':'),    (std::vector<std::string>{" ", " a\n"}      ));
  EXPECT_EQ(stk::split_string(":\t:", ':'),      (std::vector<std::string>{"", "\t", ""}     ));
  EXPECT_EQ(stk::split_string("a: :", ':'),      (std::vector<std::string>{"a", " ", ""}     ));
  EXPECT_EQ(stk::split_string("\t : \n:a", ':'), (std::vector<std::string>{"\t ", " \n", "a"}));
}

TEST(SplitCsvString, empty)
{
  EXPECT_EQ(stk::split_csv_string(""), std::vector<std::string>{});
}

TEST(SplitCsvString, whitespace)
{
  EXPECT_EQ(stk::split_csv_string("   "), std::vector<std::string>{"   "});
}

TEST(SplitCsvString, emptyFields)
{
  EXPECT_EQ(stk::split_csv_string(","),    (std::vector<std::string>{"", ""}            ));
  EXPECT_EQ(stk::split_csv_string(",,"),   (std::vector<std::string>{"", "", ""}        ));
  EXPECT_EQ(stk::split_csv_string(",,,"),  (std::vector<std::string>{"", "", "", ""}    ));
  EXPECT_EQ(stk::split_csv_string(",,,,"), (std::vector<std::string>{"", "", "", "", ""}));
}

TEST(SplitCsvString, multipleNonemptyFields)
{
  EXPECT_EQ(stk::split_csv_string("a"),           (std::vector<std::string>{"a"}                ));
  EXPECT_EQ(stk::split_csv_string("a,b"),         (std::vector<std::string>{"a", "b"}           ));
  EXPECT_EQ(stk::split_csv_string("a,b,c"),       (std::vector<std::string>{"a", "b", "c"}      ));
  EXPECT_EQ(stk::split_csv_string("a,b,c,d"),     (std::vector<std::string>{"a", "b", "c", "d"} ));
  EXPECT_EQ(stk::split_csv_string("a:b,c=d,e f"), (std::vector<std::string>{"a:b", "c=d", "e f"}));
}

TEST(SplitCsvString, mixedEmptyNonemptyFields)
{
  EXPECT_EQ(stk::split_csv_string("a,"),   (std::vector<std::string>{"a", ""}     ));
  EXPECT_EQ(stk::split_csv_string(",a"),   (std::vector<std::string>{"", "a"}     ));
  EXPECT_EQ(stk::split_csv_string("a,,b"), (std::vector<std::string>{"a", "", "b"}));
  EXPECT_EQ(stk::split_csv_string("a,,"),  (std::vector<std::string>{"a", "", ""} ));
  EXPECT_EQ(stk::split_csv_string(",,a"),  (std::vector<std::string>{"", "", "a"} ));
}

TEST(SplitCsvString, mixedEmptyNonemptyWhitespaceFields)
{
  EXPECT_EQ(stk::split_csv_string(" a , "),     (std::vector<std::string>{" a ", " "}       ));
  EXPECT_EQ(stk::split_csv_string(" , a\n"),    (std::vector<std::string>{" ", " a\n"}      ));
  EXPECT_EQ(stk::split_csv_string(",\t,"),      (std::vector<std::string>{"", "\t", ""}     ));
  EXPECT_EQ(stk::split_csv_string("a, ,"),      (std::vector<std::string>{"a", " ", ""}     ));
  EXPECT_EQ(stk::split_csv_string("\t , \n,a"), (std::vector<std::string>{"\t ", " \n", "a"}));
}
