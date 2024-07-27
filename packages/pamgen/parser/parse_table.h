// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef parse_tableH
#define parse_tableH

#include "keyword.h"
#include <vector>
namespace PAMGEN_NEVADA {

/*****************************************************************************/
//! Represents a table of keywords with their associated parsing functions.  
class Parse_Table
/*****************************************************************************/
{
  public:
    
    //! Create an empty Parse_Table to which keywords can subsequently
    //! be added via the push_back function.
    Parse_Table() {}

    //! Create a parse table with an initial set of keywords.
    Parse_Table(unsigned N, const Keyword *list) ;
    
    Parse_Table( const Parse_Table& );
    Parse_Table& operator=( const Parse_Table& );
    
    //! Inserts each unique keyword from the given table into this table.
    //! Complexity is linear in size of the given table.
    //@{
    void merge( const Keyword & key );
    void merge( const Parse_Table & table );
    void merge( const Parse_Table * table );
    void merge( const Keyword * list, unsigned N );
    //@}
    
    unsigned size() const { return keywords.size(); }

    const Keyword * begin() const { return size() > 0 ? &(keywords[0]) : 0; }
    const Keyword * end() const { return begin() + size(); }
    
    //! Produce a std::string, suitable for output, listing all the
    //! keywords in the table.
    std::string Concatenate_Legal_Commands() const;

    //! Check the table for ordering and uniqueness of keywords
    bool Check_Keyword_Table() const;

    //! Check that a keyword name is a legitimate ALEGRA keyword name.
    //! The keyword must consist of one or more words, separated
    //! by single spaces but with no preceeding or tailing spaces.
    //! Each word must consist of uppercase letters, digits, and
    //! underscores, but cannot begin with a digit. 
    static bool Check_Keyword(const char*);

    //! Check the entire keyword structure, returning true if it
    //! is no good.
    static bool Is_Invalid_Keyword(const Keyword &);

  private:
    std::vector<Keyword> keywords;
};
}//end namespace PAMGEN_NEVADA
#endif
