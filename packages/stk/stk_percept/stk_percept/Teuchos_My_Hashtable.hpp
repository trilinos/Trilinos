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

#ifndef TEUCHOS_MY_HASHTABLE_H
#define TEUCHOS_MY_HASHTABLE_H

/// NOTE: this is copied from Teuchos_Hashtable and modified to enable some stl-like interfaces, such as iterators and operator[]

/*! \file Teuchos_Hashtable.hpp
    \brief Templated hashtable class
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_HashUtils.hpp"
//#include "Teuchos_Hashtable.hpp"

namespace Teuchos
{
  using std::string;

  /** \ingroup Containers
   * \brief Helper class for Teuchos::Hashtable, representing a single <key, value> pair.
   */
  template<class Key, class Value> class HashPair
    {
    public:
      //! Empty constructor 
      inline HashPair() : key_(), value_() {;}
      //! Basic <key, value> constructor
      inline HashPair(const Key& key, const Value& value)
        : key_(key), value_(value) {;}

      //! Templated key variable
      Key key_;
      //! Templated value variable
      Value value_;

      Key& first() { return key_; }
      Value& second() { return value_; }
    };

  /**
     \ingroup Containers
     \brief Templated hashtable class.
     @author Kevin Long
  */
  template<class Key, class Value> class Teuchos_Hashtable
    {
    public:

      //! Create an empty Teuchos_Hashtable
      inline Teuchos_Hashtable(int capacity=101, double rehashDensity = 0.8);

      //! Check for the presence of a key
      inline bool containsKey(const Key& key, bool debug=false) const ;

      inline Value& operator[](const Key& key) 
      { 
        int hashCodeOfKey = hashCode(key);
        int index = hashCodeOfKey % capacity_;

        //if (index < data_.size())
        {
          //TEST_FOR_EXCEPTION(index >= data_.size(), std::runtime_error, "100 hashCodeOfKey, capacity_, index, sz=" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_) << " " << Teuchos::toString(index) << " " << Teuchos::toString(data_.size()) );
          //VERIFY_OP_ON(index, <,  data_.size(), std::runtime_error, "100" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_));
          TEST_FOR_EXCEPTION(index >= data_.capacity(), std::runtime_error, "100 hashCodeOfKey, capacity_, index, sz=" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_) << " " << Teuchos::toString(index) << " " << Teuchos::toString(data_.size()) );
          Array<HashPair<Key, Value> >& candidates = data_[index];
          //TEST_FOR_EXCEPTION(index >= data_.size() && candidates.size(), std::runtime_error, "100 hashCodeOfKey, capacity_, index, sz=" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_) << " " << Teuchos::toString(index) << " " << Teuchos::toString(data_.size()) );      
          for (int i=0; i<candidates.length(); i++)
            {
              HashPair<Key, Value>& c = candidates[i];
              if (c.key_ == key)
                {
                  return c.value_;
                }
            }
        }

        //  put(key, Value());

        {
          // no duplicate key, so increment element count by one.
          count_++;
      
          // check for need to resize.
          if ((double) count_ > rehashDensity_ * (double) capacity_)
            {
              capacity_ = HashUtils::nextPrime(capacity_+1);
              rehash();
              // recaluate index
              index = hashCodeOfKey % capacity_;
            }
      
          //TEST_FOR_EXCEPTION(index >= data_.size(), std::runtime_error, "100 hashCodeOfKey, capacity_, index, sz=" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_) << " " << Teuchos::toString(index) << " " << Teuchos::toString(data_.size()) );
          TEST_FOR_EXCEPTION(index >= data_.capacity(), std::runtime_error, "100 hashCodeOfKey, capacity_, index, sz=" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_) << " " << Teuchos::toString(index) << " " << Teuchos::toString(data_.size()) );

          Array<HashPair<Key, Value> > & data_ref = data_[index];
          //TEST_FOR_EXCEPTION(index >= data_.size() && data_ref.size(), std::runtime_error, "100 hashCodeOfKey, capacity_, index, sz=" << Teuchos::toString(hashCodeOfKey) << " " << Teuchos::toString(capacity_) << " " << Teuchos::toString(index) << " " << Teuchos::toString(data_.size()) );

          int length = data_ref.length();
          data_ref.append(HashPair<Key, Value>(key, default_value_));
          accumulateAvgFill(length+1);

          TEST_FOR_EXCEPTION(length >= data_ref.size(), std::runtime_error, "3 " << Teuchos::toString(length) << "  " << Teuchos::toString(data_ref.size()));
          return data_ref[length].value_;
        }
      }

      inline const Value& operator[](const Key& key) const { return get(key); }

      //! Get the value indexed by key
      inline const Value& get(const Key& key) const ;

      //! Put a new (key, value) pair in the table.
      inline void put(const Key& key, const Value& value, bool debug = false);

      //! Remove from the table the element given by key.
      inline void remove(const Key& key);

      //! Get the number of elements in the table
      inline int size() const {return count_;}

      //! Get lists of keys and values in Array form
      inline void arrayify(Array<Key>& keys, Array<Value>& values) const ;

      //! Return the average degeneracy (average number of entries per hash code).
      inline double avgDegeneracy() const {return avgDegeneracy_;}

      //! Return the density of the hashtable (num entries / capacity)
      inline double density() const {return ((double)count_)/((double) capacity_);}

      //! Set the density at which to do a rehash
      inline void setRehashDensity(double rehashDensity);

      //! Write to a std::string
      inline std::string toString() const ;

      //! direct access to the data array (could be avoided with a friend statement)
      inline Array<Array<HashPair<Key, Value> > >& getData() { return data_; }

      inline void clear() { data_.clear() ; }

    private:

      inline void rehash();
      inline int nextPrime(int newCap) const ;
      inline void accumulateAvgFill(int n) const ;


      Array<Array<HashPair<Key, Value> > > data_;
      int count_;
      int capacity_;
      mutable Value mostRecentValue_;
      mutable Key mostRecentKey_;
      mutable Value default_value_;

      mutable size_t nHits_;
      mutable double avgDegeneracy_;
      double rehashDensity_;

      bool debug_;

    };

  template<class Key, class Value>
  std::string toString(const Teuchos_Hashtable<Key, Value>& h);

  /** \relates Teuchos_Hashtable 
      \brief Write Teuchos_Hashtable to a stream
  */
  template<class Key, class Value>
  std::ostream& operator<<(std::ostream& os, const Teuchos_Hashtable<Key, Value>& h);

  template<class Key, class Value> inline
    Teuchos_Hashtable<Key, Value>::Teuchos_Hashtable(int capacity, double rehashDensity)
    : data_(), count_(0), capacity_(HashUtils::nextPrime(capacity)),
      nHits_(0), avgDegeneracy_(0), rehashDensity_(rehashDensity),  debug_(false)

    {
      data_.resize(capacity_);
    }

  template<class Key, class Value> inline
  bool Teuchos_Hashtable<Key, Value>::containsKey(const Key& key, bool debug) const
    {
      int hashCodeOfKey = hashCode(key);
      int index = hashCodeOfKey % capacity_;
      if (debug)
        {
          std::cout << "Teuchos_Hashtable::containsKey index= " << index << " capacity_= " << capacity_ 
                    << " hashCodeOfKey= " << hashCodeOfKey
                    << std::endl;
        }
      const Array<HashPair<Key, Value> >& candidates
        = data_[index];

      if (debug)
        {
          std::cout << "Teuchos_Hashtable::containsKey index= " << index << " capacity_= " << capacity_ 
                    << " candidates.length()= " << candidates.length()
                    << " hashCodeOfKey= " << hashCodeOfKey
                    << std::endl;
        } 
      
      for (int i=0; i<candidates.length(); i++)
        {
          const HashPair<Key, Value>& c = candidates[i];
          if (c.key_ == key)
            {
              //          (Key&) mostRecentKey_ = key;
              //(Value&) mostRecentValue_ = c.value_;
              return true;
            }
        }
      return false;
    }

  template<class Key, class Value> inline
  void Teuchos_Hashtable<Key, Value>::put(const Key& key, const Value& value, bool debug)
    {
      int index = hashCode(key) % capacity_;
      
      Array<HashPair<Key, Value> >& local = data_[index];
      
      // check for duplicate key
      for (int i=0; i<local.length(); i++)
        {
          if (local[i].key_ == key)
            {
              if (debug)
                std::cout << "debug tmp Teuchos_Hashtable::put found key=true" << std::endl;
              local[i].value_ = value;
              return;
            }
        }
      
      // no duplicate key, so increment element count by one.
      count_++;
              if (debug)
                std::cout << "debug tmp Teuchos_Hashtable::put  count_ = " << count_ << std::endl;
      
      // check for need to resize.
      if ((double) count_ > rehashDensity_ * (double) capacity_)
        {
              if (debug)
                std::cout << "debug tmp Teuchos_Hashtable::put  rehash = " << count_ << std::endl;
          capacity_ = HashUtils::nextPrime(capacity_+1);
          rehash();
          // recaluate index
          index = hashCode(key) % capacity_;
        }
      
      data_[index].append(HashPair<Key, Value>(key, value));

              if (debug)
                std::cout << "debug tmp Teuchos_Hashtable::put  data_.size() = " << data_.size() << std::endl;
    }



  template<class Key, class Value> inline
    void Teuchos_Hashtable<Key, Value>::rehash()
    {

      //      std::cout << "rehashing" << std::endl;
      Array<Array<HashPair<Key, Value> > > tmp(capacity_);

      for (int i=0; i<data_.length(); i++)
        {
          for (int j=0; j<data_[i].length(); j++)
            {
              int newIndex = hashCode(data_[i][j].key_) % capacity_;
              tmp[newIndex].append(data_[i][j]);
            }
        }
      
      data_ = tmp;
    }


  template<class Key, class Value> inline
    void Teuchos_Hashtable<Key, Value>::arrayify(Array<Key>& keys, Array<Value>& values) const
    {
      keys.reserve(size());
      values.reserve(size());

      for (int i=0; i<data_.length(); i++)
        {
          for (int j=0; j<data_[i].length(); j++)
            {
              keys.append(data_[i][j].key_);
              values.append(data_[i][j].value_);
            }
        }
    }

  template<class Key, class Value>  inline
  std::string Teuchos_Hashtable<Key, Value>::toString() const 
  {
    Array<Key> keys;
    Array<Value> values;
    arrayify(keys, values);
    
    std::string rtn = "[";
    for (int i=0; i<keys.length(); i++)
      {
        rtn += "{" + Teuchos::toString(keys[i]) + ", " + Teuchos::toString(values[i])
          + "}";
        if (i < keys.length()-1) rtn += ", ";
      }
    rtn += "]";
    
    return rtn;
  }

  template<class Key, class Value>  inline
    std::string toString(const Teuchos_Hashtable<Key, Value>& h)
    {
      Array<Key> keys;
      Array<Value> values;
      h.arrayify(keys, values);

      std::string rtn = "keys.length = " + Teuchos::toString((int)keys.length()) + " [";
      for (int i=0; i<keys.length(); i++)
        {
          rtn += "{" + Teuchos::toString(keys[i]) + ", " + Teuchos::toString(values[i])
            + "}";
          if (i < keys.length()-1) rtn += ", ";
        }
      rtn += "]";

      return rtn;
    }

  template<class Key, class Value> inline
    const Value& Teuchos_Hashtable<Key, Value>::get(const Key& key) const
    {
      TEST_FOR_EXCEPTION(!containsKey(key),
                         std::runtime_error,
                         "Teuchos_Hashtable<Key, Value>::get: key " 
                         << Teuchos::toString(key) 
                         << " not found in Teuchos_Hashtable"
                         << toString());
      
      const Array<HashPair<Key, Value> >& candidates
        = data_[hashCode(key) % capacity_];

      accumulateAvgFill(candidates.length());
      
      for (int i=0; i<candidates.length(); i++)
        {
          const HashPair<Key, Value>& c = candidates[i];
          if (c.key_ == key)
            {
              return c.value_;
            }
        }
      return mostRecentValue_;
    }


  template<class Key, class Value> inline
    void Teuchos_Hashtable<Key, Value>::remove(const Key& key)
    {
      TEST_FOR_EXCEPTION(!containsKey(key),
                         std::runtime_error,
                         "Teuchos_Hashtable<Key, Value>::remove: key " 
                         << Teuchos::toString(key) 
                         << " not found in Teuchos_Hashtable"
                         << toString());

      count_--;
      int h = hashCode(key) % capacity_;
      const Array<HashPair<Key, Value> >& candidates = data_[h];

      for (int i=0; i<candidates.length(); i++)
        {
          const HashPair<Key, Value>& c = candidates[i];
          if (c.key_ == key)
            {
              data_[h].remove(i);
              break;
            }
        }
    }

  template<class Key, class Value> inline
  void Teuchos_Hashtable<Key, Value>::accumulateAvgFill(int n) const
  {
    avgDegeneracy_ = ((double) nHits_)/(nHits_ + 1.0) * avgDegeneracy_ + ((double) n)/(nHits_ + 1.0);
    nHits_++;
    }

  template<class Key, class Value>  inline
    std::ostream& operator<<(std::ostream& os, const Teuchos_Hashtable<Key, Value>& h)
    {
      return os << toString(h);
    }


}

#endif // TEUCHOS_HASHTABLE_H
