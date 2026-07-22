// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "parse_table.h"
#include "parser.h"

#include <ctype.h>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

namespace PAMGEN_NEVADA {

Parse_Table::Parse_Table(unsigned N, const Keyword *list)
{
  assert(list != 0);  // (kind of moot by this point)
  assert(std::find_if(list, list+N, Is_Invalid_Keyword)==list+N);
  merge( list, N );
}

Parse_Table::Parse_Table( const Parse_Table& x )
 : keywords(x.keywords)
{}

Parse_Table& Parse_Table::operator=( const Parse_Table& x )
{
  keywords = x.keywords;
  return *this;
}

/*****************************************************************************/
void Parse_Table::merge( const Keyword & key )
/*****************************************************************************/
{
  merge( &key, 1 );
}

/*****************************************************************************/
void Parse_Table::merge( const Parse_Table & table )
/*****************************************************************************/
{
  merge( table.begin(), table.size() );
}

/*****************************************************************************/
void Parse_Table::merge( const Parse_Table * table )
/*****************************************************************************/
{
  merge( table->begin(), table->size() );
}

class Parse_Table_kwd_less {
public:
  bool operator()(const Keyword& x, const Keyword& y) {
    return Token::Token_Match( x.name, y.name ) < 0;
  }
};

/*****************************************************************************/
void Parse_Table::merge( const Keyword * tab, unsigned N )
/*****************************************************************************/
{
  assert( Check_Keyword_Table() == 0 );

  unsigned addsize = 0;

  for ( unsigned k = 0; k < N; ++k)
    if ( ! binary_search( begin(), end(), tab[k] ) )
      ++addsize;

  keywords.reserve( size() + addsize );

  for ( unsigned k = 0; k < N; ++k)
  {
    assert( ! Is_Invalid_Keyword( tab[k] ) );

    // the lower bound function does a binary search to find the first
    // keyword such that all previous keywords are less than tab[k]
    const Keyword * itr = lower_bound( begin(), end(), tab[k] );
    assert( itr >= begin() && itr <= end() );

    if ( itr == end() )
    {
      keywords.push_back( tab[k] );
    }
    else if ( strcmp( (*itr).name, tab[k].name ) != 0 )
    {
      assert( itr >= begin() );

      // even if the keyword is not the same, the abreviation rules for
      // token matching can cause ambiguities
      if ( binary_search(begin(), end(), tab[k], Parse_Table_kwd_less()) )
      {
        std::cout << "ERROR Parse_Table::merge() ambiguous keywords detected: "
                  << (*itr).name << ", "
                  <<  tab[k].name
                  << std::endl;
        exit(0);
      }

      // keyword not found in my list; open up a slot at the right location

      unsigned i = ( itr - begin() );
      assert( i <= size() );

      keywords.push_back( tab[k] );  // just to increase the size by one
      for ( unsigned j = size()-1; j > i; --j )
        keywords[j] = keywords[j-1];

      keywords[i] = tab[k];
    }
  }

  assert( Check_Keyword_Table() == 0 );
}

/*****************************************************************************/
string Parse_Table::Concatenate_Legal_Commands() const
/*****************************************************************************/
{
  return PAMGEN_NEVADA::Concatenate_Legal_Commands(begin(), size());
}

/*****************************************************************************/
bool Parse_Table::Check_Keyword_Table() const
/*****************************************************************************/
{
  return PAMGEN_NEVADA::Check_Keyword_Table(begin(), size());
}


/*****************************************************************************/
bool Check_Keyword_Table(const Keyword *table, int N)
/*****************************************************************************/
{
  // This function tests for valid ordering of the Keyword table

  int error_flag = 0;

  for (int i=0; i<N-1; i++) {
    string entry1 = string(table[i].name);

  for (int j=i+1; j<N; j++) {
    string entry2 = string(table[j].name);

// check one way
    int test_result = Token::Token_Match(entry1.c_str(),entry2.c_str());
    if (test_result==0) {
      error_flag += 1;
      std::cout << "Check_Keyword_Table(...): "
                <<  entry1
                << " preceding "
                << entry2
                << " is ambiguous\n ";
      exit(1);
    }
    if (test_result>0) {
      error_flag += 1;
      std::cout << "Check_Keyword_Table(...): "
                <<  entry1
                << " preceding "
                << entry2
                << " is out of place in table\n ";
      exit(1);
    }

// check the other way
    test_result = Token::Token_Match(entry2.c_str(),entry1.c_str());
    if (test_result==0) {
      error_flag += 1;
      std::cout << "Check_Keyword_Table(...): "
                << entry1
                << " following "
                << entry2
                << " would be ambiguous\n";
      exit(1);
    }
    if (test_result<0) {
      error_flag += 1;
      std::cout << "Check_Keyword_Table(...): "
                << entry2
                << " following "
                << entry1
                << " is out of place in table\n";
      exit(1);
    }
  }
  }

  return error_flag;
}

/*****************************************************************************/
bool Parse_Table::Check_Keyword(const char *candidate)
/*****************************************************************************/
{
  char c = *candidate;
  if (isalpha(c) || c=='_'){
    while (1){
      while (isalnum(c) || c=='_'){
        c = *candidate++;
      }
      if (c==0) return true;
      if (c!=' ') return false;
      c = *candidate++;
    }
  }

  return false;
}


/*****************************************************************************/
bool Parse_Table::Is_Invalid_Keyword(const Keyword &key)
/*****************************************************************************/
{
  if (!key.name) return true;
  return !Check_Keyword(key.name) || !key.func;
}

}//end namespace
