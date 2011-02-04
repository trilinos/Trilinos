#ifndef stk_adapt_SubDimCell_hpp
#define stk_adapt_SubDimCell_hpp

#include <set>
#include <iostream>


#define STK_ADAPT_SUBDIMCELL_USES_STL_SET 0
#define STK_ADAPT_SUBDIMCELL_USES_STL_VECTOR 0
#define STK_ADAPT_SUBDIMCELL_USES_NO_MALLOC_ARRAY 1

#if STK_ADAPT_SUBDIMCELL_USES_NO_MALLOC_ARRAY
#include <stk_percept/NoMallocArray.hpp>
#endif

#if STK_ADAPT_SUBDIMCELL_USES_BOOST_ARRAY
#include <boost/array.hpp>
#endif

#include <vector>
#include <algorithm>
#include <iostream>

#ifdef STK_HAVE_TBB
#include <tbb/scalable_allocator.h>
#endif


namespace stk {
  namespace adapt {

    // only set one of these

#if STK_ADAPT_SUBDIMCELL_USES_STL_SET

#  ifdef STK_HAVE_TBB
  typedef std::set<unsigned, std::less<unsigned>, tbb::scalable_allocator<unsigned> > SubDimCellBaseClass;
#  else
  typedef std::set<unsigned> SubDimCellBaseClass;
#  endif

    /// this class represents the set of node id's defining a sub-dimensional entity of an element (like a face or edge)
    template<typename Ids=unsigned, std::size_t N=4>
    class SubDimCell : public SubDimCellBaseClass
    {
      std::size_t m_hash;
    public:
      //set<Ids> m_ids;
      SubDimCell() : SubDimCellBaseClass(), m_hash(0u) {}
      SubDimCell(unsigned num_ids, Ids *ids) : SubDimCellBaseClass(ids, ids+num_ids), m_hash(0u)
      {
      }
      unsigned getHash() const
      {
        return m_hash;
      }
      void setHash(std::size_t hash)
      {
        m_hash = hash;
      }
      void clear() 
      {
        m_hash = 0u;
        SubDimCellBaseClass::clear();
      }
    };
#endif

#if STK_ADAPT_SUBDIMCELL_USES_STL_VECTOR

    //typedef array<int, 3> SubDimCell;

    /// We assume we don't have any sub-dimensional entities with more than 4 nodes
    template<class T, std::size_t N=4>
    class SubDimCell : public std::vector<T> //: public my_array<T,4> //: public boost::array<T,4>
    {
    public:
      //typedef boost::array<T,4> base_type;
      typedef std::vector<T> base_type;
      
      SubDimCell(unsigned n=4) : base_type() { base_type::reserve(n); }
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
            std::sort( base_type::begin(), base_type::end() );
          }
      }

    };
#endif

#if STK_ADAPT_SUBDIMCELL_USES_NO_MALLOC_ARRAY

    //typedef array<int, 3> SubDimCell;

    /// We assume we don't have any sub-dimensional entities with more than 4 nodes
    template<class T, std::size_t N=4>
    class SubDimCell : public stk::percept::NoMallocArray<T,N>
    {
      std::size_t m_hash;
      
    public:
      typedef stk::percept::NoMallocArray<T,N> base_type;
      typedef std::size_t    size_type;

      typedef SubDimCell<T,N> VAL;
      
      //repo always init to 0 size: SubDimCell(unsigned n=4) : base_type(n), m_hash(0u) {}
      SubDimCell() : base_type(), m_hash(0u) {}
      SubDimCell(unsigned n) : base_type(), m_hash(0u) {}

      
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
            //if (size() > max_size() ) throw std::runtime_error("SubDimCell out of range");
            base_type::insert(val);
            std::sort( base_type::begin(), base_type::end() );
          }
        m_hash = hashCode();
      }

      int hashCode()
      {
        std::size_t sum = 0;

        for (typename base_type::const_iterator i = this->begin(); i != this->end(); i++)
          {
            sum += static_cast<std::size_t>(*i);
          }
        return sum;
      }

      unsigned getHash() const
      {
        return m_hash;
      }
      void setHash(std::size_t hash)
      {
        m_hash = hash;
      }
      void clear() 
      {
        m_hash = 0u;
        base_type::clear();
      }

      bool operator==(const VAL& rhs) const
      {
        if (base_type::size() != rhs.size())
          return false;
        //return true;
        for (size_type i = 0; i < base_type::size(); i++)
          {
            if ((*this)[i] != rhs[i])
              return false;
          }
        return true;
      }

      bool operator<(const VAL& rhs) const
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

    };
#endif

#if 1

#define DATA
#define GET(x,i) x[i]

    template<class T, std::size_t N=4>
    struct SubDimCell_compare 
    {
      bool operator() (const SubDimCell<T,N>& lhs, const SubDimCell<T,N>& rhs) const
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
    struct my_hash : public std::unary_function< SubDimCell<T,N>, std::size_t>
    {
      typedef SubDimCell<T,N> _Tp ;

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
    struct my_equal_to :  public std::binary_function<SubDimCell<T,N>,  
                                                      SubDimCell<T,N>, bool>
    {
      typedef SubDimCell<T,N> _Tp ;
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
    struct my_hash_old : public std::unary_function< SubDimCell<T,N> , std::size_t>
    {
      typedef SubDimCell<T,N> _Tp ;

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
    struct my_equal_to_old :  public std::binary_function< SubDimCell<T,N>,  SubDimCell<T,N>, bool>
    {
      typedef SubDimCell<T,N> _Tp ;
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

//     template<class T, std::size_t N>
//     std::ostringstream& operator<< (std::ostringstream& out, SubDimCell<T,N>& c)
//     {
//       out << t;
//       return out;
//     }

#undef DATA
#undef GET

#endif


  }
}

#endif
