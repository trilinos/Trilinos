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
