// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "InputBlock.h"

#include <vector>

using namespace std;

namespace PAMGEN_NEVADA{


InputBlock::InputBlock() : lineno(0), parent(0) {}

InputBlock::~InputBlock()
{
  InputBlock_iterator itr = blocks.begin();
  for (; itr != blocks.end(); ++itr)
    delete (*itr);
}

const string& InputBlock::getName() const
{
  return name;
}

bool InputBlock::empty() const
{
  return attrs.empty() && tokens.empty() && blocks.empty();
}

unsigned InputBlock::numAttributes() const
{
  return attrs.size();
}

InputBlock::Attribute_iterator InputBlock::beginAttributes() const
{
  return attrs.begin();
}

InputBlock::Attribute_iterator InputBlock::endAttributes() const
{
  return attrs.end();
}

InputBlock::Attribute_iterator InputBlock::findAttribute(const char* n) const
{
  string sn;
  if (n != 0) sn = n;
  
  Attribute_iterator itr = attrs.begin();
  for (; itr != attrs.end(); ++itr)
    if ( (*itr).first == sn )
      return itr;
  
  return endAttributes();
}

unsigned InputBlock::size() const
{
  return tokens.size();
}

InputBlock::const_iterator InputBlock::begin() const
{
  return tokens.begin();
}

InputBlock::const_iterator InputBlock::end() const
{
  return tokens.end();
}

const string& InputBlock::lookAhead( InputBlock::const_iterator current,
                                     InputBlock::const_iterator endtoks )
{
  static string none;
  
  if (current == endtoks) return none;

  ++current;
  if (current == endtoks) return none;

  return *current;
}

static void IB__split( const string & s, vector<string> & L )
{
  L.resize(0);
  
  const char* sp = s.c_str();
  unsigned len = 0;
  
  for (; *(sp+len) != 0;)
  {
//    cout << "+" << (*(sp+len));
    if ( *(sp+len) == ' ' )
    {
      if (len > 0) {
//        cout << "push" << string(sp,len);
        L.push_back( string(sp,len) );
      }
      sp += (len+1);
      len = 0;
    }
    else ++len;
  }
  
  if (len > 0) {
//    cout << "push" << string(sp,len);
    L.push_back( string(sp,len) );
  }
  
//  cout << endl;
}

bool InputBlock::abbreviation( const string & s, const string & m,
                               unsigned num_chars )
{
//  cout << "*** " << s << ", " << m << endl;
  
  static vector<string> mL;  // not thread safe but nicer on memory
  IB__split( m, mL );
  
  static vector<string> sL;
  IB__split( s, sL );
  
//   cout << "len mL " << mL.size() << " len sL " << sL.size() << endl;
//   for (unsigned i = 0; i < mL.size(); ++i)
//     cout << "    " << mL[i];
//   cout << "  --  ";
//   for (unsigned i = 0; i < sL.size(); ++i)
//     cout << "    " << sL[i];
//   cout << endl;
  
  if (mL.size() != sL.size()) return false;
  
  for (unsigned i = 0; i < mL.size(); ++i)
  {
    const string& mS = mL[i];
    const string& sS = sL[i];
    
    if (mS.size() <= 3) {
      if (mS != sS) return false;
    }
    else {
      if (mS.compare(0, num_chars, sS, 0, num_chars) != 0) return false;
    }
  }
  
  return true;
}

unsigned InputBlock::numInputBlocks() const
{
  return blocks.size();
}

InputBlock::InputBlock_iterator InputBlock::beginInputBlocks() const
{
  return blocks.begin();
}

InputBlock::InputBlock_iterator InputBlock::endInputBlocks() const
{
  return blocks.end();
}

InputBlock::InputBlock_iterator InputBlock::findInputBlock(const char* n) const
{
  string sn;
  if (n != 0) sn = n;
  
  InputBlock_iterator itr = blocks.begin();
  for (; itr != blocks.end(); ++itr)
    if ( (*itr)->name == sn )
      return itr;
  
  return endInputBlocks();
}

void InputBlock::set( const char* n, int ln )
{
  if (n != 0) name = n;
  else name = "";
  lineno = ln;
}

void InputBlock::add( const char* n, const char* v )
{
  if (n == 0) n = "";
  if (v == 0) v = "";
  
  attrs.push_back( Attribute(n,v) );
}

void InputBlock::add( const char* t )
{
  if (t == 0) t = "";
  
  tokens.push_back(t);
}

void InputBlock::add( InputBlock* block )
{
  if (block != 0)
  {
    blocks.push_back( block );
    block->parent = this;
  }
}

std::list< InputBlock* >::iterator InputBlock::beginInputBlocks()
{
  return blocks.begin();
}

std::list< InputBlock* >::iterator InputBlock::endInputBlocks()
{
  return blocks.end();
}

void InputBlock::Display( ostream& os, int indent ) const
{
  string pad;
  for (int i = 0; i < indent; ++i)
    pad += "  ";
  
  os << pad << "Block: " << name << ", line " << lineno << "\n";
  if (size() > 0) {
    os << pad;
    for (const_iterator itr = begin(); itr != end(); ++itr) {
      if (itr == begin()) os << "    " << (*itr);
      else                os << ", " << (*itr);
    }
    os << endl;
  }
  
  for (InputBlock_iterator itr = beginInputBlocks();
       itr != endInputBlocks(); ++itr)
    (*itr)->Display( os, indent+1 );
}

}// end namespace
