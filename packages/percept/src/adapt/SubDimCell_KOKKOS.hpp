// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_SubDimCell_KOKKOS_hpp
#define adapt_SubDimCell_KOKKOS_hpp

#include <set>
#include <iostream>


#define STK_ADAPT_SUBDIMCELL_KOKKOS_USES_STL_SET 0
#define STK_ADAPT_SUBDIMCELL_KOKKOS_USES_STL_VECTOR 0
#define STK_ADAPT_SUBDIMCELL_KOKKOS_USES_NO_MALLOC_ARRAY 1

#if STK_ADAPT_SUBDIMCELL_KOKKOS_USES_NO_MALLOC_ARRAY
#include <percept/NoMallocArray.hpp>
#endif

#if STK_ADAPT_SUBDIMCELL_KOKKOS_USES_BOOST_ARRAY
#include <percept/PerceptBoostArray.hpp>
#endif

#include <vector>
#include <algorithm>
#include <iostream>

#ifdef STK_HAVE_TBB
#include <tbb/scalable_allocator.h>
#endif


  namespace percept {


    template<class T>
    struct SubDimCell_KOKKOSCompare
    {
      bool operator() (T i, T j) { return (i < j) ; }
    };

    // only set one of these

#if STK_ADAPT_SUBDIMCELL_KOKKOS_USES_STL_SET

#  ifdef STK_HAVE_TBB
  typedef std::set<unsigned, std::less<unsigned>, tbb::scalable_allocator<unsigned> > SubDimCell_KOKKOSBaseClass;
#  else
  typedef std::set<unsigned> SubDimCell_KOKKOSBaseClass;
#  endif

    /// this class represents the set of node id's defining a sub-dimensional entity of an element (like a face or edge)
    template<typename Ids=unsigned, std::size_t N=4, class CompareClass = SubDimCell_KOKKOSCompare<T> >
    class SubDimCell_KOKKOS : public SubDimCell_KOKKOSBaseClass
    {
      std::size_t m_hash;
    public:
      //set<Ids> m_ids;
      SubDimCell_KOKKOS() : SubDimCell_KOKKOSBaseClass(), m_hash(0u) {}
      SubDimCell_KOKKOS(unsigned num_ids, Ids *ids) : SubDimCell_KOKKOSBaseClass(ids, ids+num_ids), m_hash(0u)
      {
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
        SubDimCell_KOKKOSBaseClass::clear();
      }
    };
#endif

#if STK_ADAPT_SUBDIMCELL_KOKKOS_USES_STL_VECTOR

    //typedef array<int, 3> SubDimCell_KOKKOS;

    /// We assume we don't have any sub-dimensional entities with more than 4 nodes
    template<class T, std::size_t N=4, class CompareClass = SubDimCell_KOKKOSCompare<T> >
    class SubDimCell_KOKKOS : public std::vector<T> //: public my_array<T,4> //: public boost::array<T,4>
    {
    public:
      //typedef boost::array<T,4> base_type;
      typedef std::vector<T> base_type;

      SubDimCell_KOKKOS(unsigned n=4) : base_type() { base_type::reserve(n); }
      void insert(T val)
      {
        bool found = false;
        for (unsigned i = 0; i < base_type::size(); i++)
          {
            if (val == (*this)[i])
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            base_type::push_back(val);
            std::sort( base_type::begin(), base_type::end(), CompareClass() );

          }
      }
      void sort()
      {
        std::sort( base_type::begin(), base_type::end(), CompareClass() );
      }
    };
#endif

#if STK_ADAPT_SUBDIMCELL_KOKKOS_USES_NO_MALLOC_ARRAY

    //typedef array<int, 3> SubDimCell_KOKKOS;



    template<class T, std::size_t N=4>
    class SDCHashCode_KOKKOS
    {
    public:
      typedef percept::NoMallocArray<T,N> base_type;
      int operator()(base_type& sdc);
    };

    /// We assume we don't have any sub-dimensional entities with more than 4 nodes
    template<class T, std::size_t N=4, class CompareClass = SubDimCell_KOKKOSCompare<T>, class HC = SDCHashCode_KOKKOS<T,N>  >
    class SubDimCell_KOKKOS : public percept::NoMallocArray<T,N>
    {
    protected:
      std::size_t m_hash;
    public:
      HC m_HashCode;
      CompareClass m_CompareClass;
      typedef percept::NoMallocArray<T,N> base_type;
      typedef std::size_t    size_type;

      typedef SubDimCell_KOKKOS<T,N,CompareClass,HC> VAL;

      //repo always init to 0 size: SubDimCell_KOKKOS(unsigned n=4) : base_type(n), m_hash(0u) {}
      SubDimCell_KOKKOS() : base_type(), m_hash(0u), m_HashCode(), m_CompareClass() {}
      SubDimCell_KOKKOS(unsigned n) : base_type(), m_hash(0u), m_HashCode(), m_CompareClass() {}
#if 0
      SubDimCell_KOKKOS(const SubDimCell_KOKKOS& sdc) : base_type(sdc), m_hash(sdc.m_hash) {}
#endif
      // behaves like std::set
      void insert(T val)
      {
        bool found = false;
        for (size_type i = 0; i < base_type::size(); i++)
          {
            if (val == (*this)[i])
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            //if (size() > max_size() ) throw std::runtime_error("SubDimCell_KOKKOS out of range");
            base_type::insert(val);
            sort();
          }
        //m_hash = hashCode();
        updateHashCode();
      }

      void sort()
      {
#ifdef __GNUC__
#if (__GNUC__ == 4) && (__GNUC_MINOR__ == 4)
            std::stable_sort( base_type::begin(), base_type::end(), m_CompareClass );
#else
            std::sort( base_type::begin(), base_type::end(), m_CompareClass );
#endif
#else
            std::sort( base_type::begin(), base_type::end(), m_CompareClass );
#endif
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
        base_type::clear();
      }

      bool operator==(const VAL& rhs) const;
      bool operator!=(const VAL& rhs) const
      {
        return !operator==(rhs);
      }

      bool operator<(const VAL& rhs) const;

    };

    template<class T, std::size_t N>
    inline int SDCHashCode_KOKKOS<T,N>::operator()(SDCHashCode_KOKKOS<T,N>::base_type& sdc)
    {
      std::size_t sum = 0;

      for (typename base_type::iterator i = sdc.begin(); i != sdc.end(); i++)
        {
          //sum += static_cast<std::size_t>(const_cast<T>(*i));
          //sum += static_cast<std::size_t>((*i)->identifier());
          sum += (size_t)(*i);
        }
      return sum;
    }


    template<class T, std::size_t N, class CompareClass, class HC >
    inline bool SubDimCell_KOKKOS<T,N,CompareClass,HC>::
    operator==(const VAL& rhs) const
    {
      if (base_type::size() != rhs.size())
        return false;
      for (size_type i = 0; i < base_type::size(); i++)
        {
          if ((*this)[i] != rhs[i])
            return false;
        }
      return true;
    }

    template<class T, std::size_t N, class CompareClass, class HC >
    inline bool SubDimCell_KOKKOS<T,N,CompareClass, HC>::
    operator<(const VAL& rhs) const
    {
      if (base_type::size() < rhs.size())
        return true;
      else if (base_type::size() > rhs.size())
        return false;
      else
        {
          for (size_type i = 0; i < base_type::size(); i++)
            {
              if ((*this)[i] < rhs[i])
                return true;
            }
        }
      return false;
    }


#endif

#if 1

#define DATA
#define GET(x,i) x[i]

    template<class T, std::size_t N=4>
    struct SubDimCell_KOKKOS_compare
    {
      bool operator() (const SubDimCell_KOKKOS<T,N>& lhs, const SubDimCell_KOKKOS<T,N>& rhs) const
      {
        if (lhs.size() < rhs.size()) return true;
        if (lhs.size() > rhs.size()) return false;
        for (unsigned i = 0; i < lhs.size(); i++)
          {
            if (GET(lhs,i) < GET(rhs,i)) return true;
            if (GET(lhs,i) > GET(rhs,i)) return false;
          }
        return false;
      }
    };

    template<class T, std::size_t N>
    struct my_hash_KOKKOS : public std::unary_function< SubDimCell_KOKKOS<T,N>, std::size_t>
    {
      typedef SubDimCell_KOKKOS<T,N> _Tp ;

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
    struct my_fast_hash_KOKKOS : public std::unary_function< SubDimCell_KOKKOS<T,N>, std::size_t>
    {
      typedef SubDimCell_KOKKOS<T,N> _Tp ;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        return x.getHash();
      }

    };

    template<class T, std::size_t N>
    struct my_equal_to_KOKKOS :  public std::binary_function<SubDimCell_KOKKOS<T,N>,
                                                      SubDimCell_KOKKOS<T,N>, bool>
    {
      typedef SubDimCell_KOKKOS<T,N> _Tp ;
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
    struct my_fast_equal_to_KOKKOS :  public std::binary_function<SubDimCell_KOKKOS<T,N>,
                                                           SubDimCell_KOKKOS<T,N>, bool>
    {
      typedef SubDimCell_KOKKOS<T,N> _Tp ;
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
    struct my_hash_KOKKOS_old : public std::unary_function< SubDimCell_KOKKOS<T,N> , std::size_t>
    {
      typedef SubDimCell_KOKKOS<T,N> _Tp ;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        std::size_t sum = 0;
        for (unsigned i = 0; i < x.size(); i++)
          {
            sum += static_cast<std::size_t>(GET(x,i));
          }
        return sum;
      }
    };

    template<class T, std::size_t N>
    struct my_equal_to_KOKKOS_old :  public std::binary_function< SubDimCell_KOKKOS<T,N>,  SubDimCell_KOKKOS<T,N>, bool>
    {
      typedef SubDimCell_KOKKOS<T,N> _Tp ;
      bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.size() != y.size()) return false;
        for (unsigned i = 0; i < x.size(); i++)
          {
            if (GET(x,i) != GET(y,i)) return false;
          }
        return true;
      }
    };


    template<class T, std::size_t N>
    std::ostream& operator<<(std::ostream& out, const SubDimCell_KOKKOS<T,N>& c)
    {
      out << "SubDimCell_KOKKOS size= " << c.size() << " vals= {";
      //for (unsigned i = 0; i < c.size(); i++)
      typename  SubDimCell_KOKKOS<T,N>::const_iterator i = c.begin();
        //for (unsigned i = 0; i < x.size(); i++)
        for (i = c.begin(); i != c.end(); i++)
          {

            out << *i << " ";
            //out << c[i] << " ";
        }
        out << "}";
      return out;
    }

//     template<class T, std::size_t N>
//     std::ostringstream& operator<< (std::ostringstream& out, SubDimCell_KOKKOS<T,N>& c)
//     {
//       out << t;
//       return out;
//     }

#undef DATA
#undef GET

#endif


  }

#endif
