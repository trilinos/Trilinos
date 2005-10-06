// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_OSTAB_HPP
#define TEUCHOS_OSTAB_HPP

#include "Teuchos_TestForException.hpp"

namespace Teuchos {

/** \brief Tabbing class for helping to create formated output.
 */
class OSTab
{
public:

  OSTab( const int tag = 0 )
    :tag_(tag)
    {
      TEST_FOR_EXCEPT( !( 0 <= tag && tag < static_cast<int>(tabIndent_.size()) ) );
      ++tabIndent_[tag_];
    }

  OSTab( const OSTab &osTab )
    :tag_(osTab.tag_)
    {
      ++tabIndent_[tag_];
    }

  ~OSTab()
    {
      --tabIndent_[tag_];
    }

  OSTab& operator=( const OSTab &osTab )
    {
      tag_ = osTab.tag_;
      ++tabIndent_[tag_];
      return *this;
    }

  static int getCurrTopTag()
    {
      return currTopTag_;
    }

  static int getNextTag()
    {
      ++currTopTag_;
      tabIndent_.push_back(0);
      return currTopTag_;
    }
  
  static int getCurrIndent( const int tag = 0 )
    {
      TEST_FOR_EXCEPT( !( 0 <= tag && tag < static_cast<int>(tabIndent_.size()) ) );
      return tabIndent_[tag];
    }

private:

  // //////////////////////
  // Private types

  typedef std::vector<int> tabIndent_t;

  // /////////////////////
  // Private data members

  static tabIndent_t    tabIndent_;
  static int            currTopTag_;

  int tag_;
  
};

} // namespace Teuchos

#endif // TEUCHOS_OSTAB_HPP
