// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef InputBlock_hh
#define InputBlock_hh

#include <string>
#include <list>
#include <iostream>

namespace PAMGEN_NEVADA{


class InputBlock {
public:
  
  typedef std::pair< std::string, std::string > Attribute;
  
  //! The iterators are for the strings/tokens contained in the InputBlock.
  //@{
  typedef std::list< std::string >::iterator       iterator;
  typedef std::list< std::string >::const_iterator const_iterator;
  //@}
  
  typedef std::list< Attribute   >::const_iterator Attribute_iterator;
  typedef std::list< InputBlock* >::const_iterator InputBlock_iterator;
  
  InputBlock();
 ~InputBlock();
  
  const std::string& getName() const;
  
  //! Returns true if no attributes, content strings, or sub-blocks.
  bool empty() const;
  
  //! These refer to the content strings contained in the InputBlock.
  //@{
  unsigned       size()  const;
  const_iterator begin() const;
  const_iterator end()   const;
  //@}
  
  //! Static convenience method to look ahead to the next token.  If there is
  //! no next token, an empty string is returned.
  static const std::string& lookAhead( const_iterator current,
                                       const_iterator endtoks );
  
  //! The first three letters of each whitespace separate string are compared.
  //! The master string is the larger, unabbreviated string.  The last
  //! argument is the number of characters to force equality with on the
  //! abbreviation.
  static bool abbreviation( const std::string & str,
                            const std::string & master,
                            unsigned num_chars = 3 );
  
  unsigned           numAttributes()  const;
  Attribute_iterator beginAttributes() const;
  Attribute_iterator endAttributes()   const;
  
  //! Returns endAttributes() if not found.
  Attribute_iterator findAttribute( const char* attr_name ) const;
  
  unsigned            numInputBlocks()  const;
  InputBlock_iterator beginInputBlocks() const;
  InputBlock_iterator endInputBlocks()   const;
  
  //! Returns endInputBlocks() if not found.
  InputBlock_iterator findInputBlock( const char* block_name ) const;
  
  //! Returns NULL if no parent.
  InputBlock* getParent() { return parent; }
  
  /*****************************************************************/
  
  // construction methods
  //@{
  
  void set( const char* block_name, int line_number );
  
  void add( const char* attr_name, const char* attr_value );
  void add( const char* token );
  //! Ownership of the block memory is transferred to this object.
  void add( InputBlock* block );
  
  std::list< InputBlock* >::iterator beginInputBlocks();
  std::list< InputBlock* >::iterator endInputBlocks();
  
  void Display( std::ostream&, int indent_level=0 ) const;
  
  //@}
  
private:
  
  std::string name;
  int         lineno;
  
  std::list< Attribute   > attrs;
  std::list< std::string > tokens;
  std::list< InputBlock* > blocks;
  
  InputBlock* parent;
};
} // end namespace 
#endif
