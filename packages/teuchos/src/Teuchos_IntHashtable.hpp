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

#ifndef INTHASHTABLE_H
#define INTHASHTABLE_H

#include "Teuchos_ConfigDefs.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_HashUtils.hpp"


namespace Teuchos
{

  class IntPair
    {
    public:
      IntPair() : key_(0), value_(0) {;}
      IntPair(int key, int value) : key_(key), value_(value) {;}

      int key_;
      int value_;
    };

  /**
   * \ingroup Containers
   * Hashtable hardwired for integers
   */

  class IntHashtable
    {
    public:
      IntHashtable(int capacity=101);

      bool containsKey(int key) const ;
      int get(int key) const ;

      void put(int key, int value) ;

      int size() const {return count_;}


      string toString() const ;
      const Array<Array<IntPair> >& data() const {return data_;}
    private:
      void rehash();
      int nextPrime(int newCap) const ;

      Array<Array<IntPair> > data_;
      int count_;
      int capacity_;
      mutable int mostRecentValue_;
      mutable int mostRecentKey_;
    };

  inline string toString(const IntHashtable& h)
    {
      return h.toString();
    }


  inline ostream& operator<<(ostream& os, const IntHashtable& h)
    {
      return os << h.toString();
    }

}



#endif
