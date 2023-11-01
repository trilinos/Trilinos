// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_SubDimCell_hpp
#define adapt_SubDimCell_hpp

#include <set>
#include <iostream>

#include <percept/NoMallocArray.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

  namespace percept {

    template<class T>
    struct SubDimCellCompare
    {
      bool operator() (T i, T j) { return (i < j) ; }
    };

    template<class T, std::size_t N=4>
    class SDCHashCode
    {
    public:
      template <T,N>
      using base_type = percept::NoMallocArray<T,N>;

      int operator()(base_type<T,N>& sdc);
    };

    /// We assume we don't have any sub-dimensional entities with more than 4 nodes
    template<class T, std::size_t N=4, class CompareClass = SubDimCellCompare<T>, class HC = SDCHashCode<T,N>  >
    class SubDimCell : public percept::NoMallocArray<T,N>
    {
    protected:
      std::size_t m_hash;
    public:
      HC m_HashCode;
      CompareClass m_CompareClass;

      template<class T, std::size_t N=4>
      using based_type = percept::NoMallocArray<T,N>;

      using size_type = std::size_t; 

      template<class T, std::size_t N=4, class CompareClass = SubDimCellCompare<T>, class HC = SDCHashCode<T,N>  >
      using VAL = SubDimCell<T,N,CompareClass,HC>;

      //repo always init to 0 size: SubDimCell(unsigned n=4) : base_type<T,N,CompareClass,HC>(n), m_hash(0u) {}
      SubDimCell() : base_type<T,N,CompareClass,HC>(), m_hash(0u), m_HashCode(), m_CompareClass() {}
      SubDimCell(unsigned n) : base_type<T,N,CompareClass,HC>(), m_hash(0u), m_HashCode(), m_CompareClass() {}
      // behaves like std::set
      void insert(T val)
      {
        bool found = false;
        for (size_type i = 0; i < base_type<T,N,CompareClass,HC>::size(); i++)
          {
            if (val == (*this)[i])
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            base_type<T,N,CompareClass,HC>::insert(val);
            sort();
          }
        updateHashCode();
      }

      void sort()
      {
        std::sort( base_type<T,N,CompareClass,HC>::begin(), base_type<T,N,CompareClass,HC>::end(), m_CompareClass );
      }

      void updateHashCode()
      {
        m_hash = hashCode();
      }
      int hashCode()
      {
        return m_HashCode(*this);
      }

      inline unsigned getHash() const
      {
        return m_hash;
      }
      inline void setHash(std::size_t hash)
      {
        m_hash = hash;
      }
      void clear()
      {
        m_hash = 0u;
        base_type<T,N,CompareClass,HC>::clear();
      }

      bool operator==(const VAL& rhs) const;
      bool operator!=(const VAL& rhs) const
      {
        return !operator==(rhs);
      }

      bool operator<(const VAL& rhs) const;

    };

    template<class T, std::size_t N>
    inline int SDCHashCode<T,N>::operator()(SDCHashCode<T,N>::base_type<T,N>& sdc)
    {
      std::size_t sum = 0;

      for (typename base_type<T,N>::iterator i = sdc.begin(); i != sdc.end(); i++)
        {
          sum += (size_t)(*i);
        }
      return sum;
    }

    template<class T, std::size_t N, class CompareClass, class HC >
    inline bool SubDimCell<T,N,CompareClass,HC>::
    operator==(const VAL& rhs) const
    {
      if (base_type<T,N,CompareClass,HC>::size() != rhs.size())
        return false;
      for (size_type i = 0; i < base_type<T,N,CompareClass,HC>::size(); i++)
        {
          if ((*this)[i] != rhs[i])
            return false;
        }
      return true;
    }

    template<class T, std::size_t N, class CompareClass, class HC >
    inline bool SubDimCell<T,N,CompareClass, HC>::
    operator<(const VAL& rhs) const
    {
      if (base_type<T,N,CompareClass,HC>::size() < rhs.size())
        return true;
      else if (base_type<T,N,CompareClass,HC>::size() > rhs.size())
        return false;
      else
        {
          for (size_type i = 0; i < base_type<T,N,CompareClass,HC>::size(); i++)
            {
              if ((*this)[i] < rhs[i])
                return true;
            }
        }
      return false;
    }

#define GET(x,i) x[i]

    template<class T, std::size_t N>
    struct my_hash
    {
      using _Tp = SubDimCell<T,N>;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        if (x.getHash())
          {
            return x.getHash();
          }
        std::size_t sum = 0;
        typename  _Tp::const_iterator i = x.begin();

        for (i = x.begin(); i != x.end(); i++)
          {
            sum += static_cast<std::size_t>(*i);
          }

        (const_cast<_Tp *>(&x))->setHash(sum);

        return sum;
      }

    };

    template<class T, std::size_t N>
    struct my_fast_hash
    {
      using _TP = SubDimCell<T,N>;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        return x.getHash();
      }

    };

    template<class T, std::size_t N>
    struct my_equal_to
    {
      using _Tp = SubDimCell<T,N>;

      bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.size() != y.size()) return false;
        //for (unsigned i = 0; i < x.size(); i++)
        typename  _Tp::const_iterator ix = x.begin();
        typename  _Tp::const_iterator iy = y.begin();
        //for (unsigned i = 0; i < x.size(); i++)
        for (ix = x.begin(), iy = y.begin(); ix != x.end(); ix++, iy++)
          {
            //if (GET(x,i) != GET(y,i)) return false;
            if (*ix != *iy) return false;
          }
        return true;
      }
    };

    template<class T, std::size_t N>
    struct my_fast_equal_to
    {
      using _Tp = SubDimCell<T,N>;

      inline bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.getHash() != y.getHash()) return false;
        if (x.size() != y.size()) return false;
        typename  _Tp::const_iterator ix = x.begin();
        typename  _Tp::const_iterator iy = y.begin();
        //for (ix = x.begin(), iy = y.begin(); ix != x.end(); ix++, iy++)
        for (; ix != x.end(); ix++, iy++)
          {
            if (*ix != *iy) return false;
          }
        return true;
      }
    };

    template<class T, std::size_t N>
    std::ostream& operator<<(std::ostream& out, const SubDimCell<T,N>& c)
    {
      out << "SubDimCell size= " << c.size() << " vals= {";
      //for (unsigned i = 0; i < c.size(); i++)
      typename  SubDimCell<T,N>::const_iterator i = c.begin();
        //for (unsigned i = 0; i < x.size(); i++)
        for (i = c.begin(); i != c.end(); i++)
          {

            out << *i << " ";
            //out << c[i] << " ";
        }
        out << "}";
      return out;
    }

#undef GET

  }

#endif
