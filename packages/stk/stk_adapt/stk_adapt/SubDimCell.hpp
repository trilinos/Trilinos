#ifndef stk_adapt_SubDimCell_hpp
#define stk_adapt_SubDimCell_hpp

#include <set>
#include <iostream>


#include <boost/array.hpp>

#include <vector>
#include <algorithm>
#include <iostream>



namespace stk {
  namespace adapt {

#define SUBDIMCELL_USES_STL_SET 1
#if SUBDIMCELL_USES_STL_SET
    /// this class represents the set of node id's defining a sub-dimensional entity of an element (like a face or edge)
    template<typename Ids=unsigned, std::size_t N=4>
    class SubDimCell : public std::set<Ids>
    {
      std::size_t m_hash;
    public:
      typedef std::set<Ids> SetIds;
      //set<Ids> m_ids;
      SubDimCell() : SetIds(), m_hash(0u) {}
      SubDimCell(unsigned num_ids, Ids *ids) : SetIds(ids, ids+num_ids), m_hash(0u)
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

    };
#else
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

#if 1
//#define DATA .data()
#define DATA
#define GET(x,i) x[i]
//#define GET(x,i) x.data()[i]
//#define GET(x,i) *( ((int *)&x) + i)

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

#if 1
        if (x.getHash())
          {
            return x.getHash();
          }
#endif
        std::size_t sum = 0;
        typename  _Tp::const_iterator i = x.begin();
        //for (unsigned i = 0; i < x.size(); i++)
        for (i = x.begin(); i != x.end(); i++)
          {
            //sum += static_cast<std::size_t>(GET(x,i));
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
    std::ostream& operator<<(std::ostream& out, SubDimCell<T,N>& c)
    {
      out << "SubDimCell size= " << c.size() << " vals= ";
      //for (unsigned i = 0; i < c.size(); i++)
      typename  SubDimCell<T,N>::const_iterator i = c.begin();
        //for (unsigned i = 0; i < x.size(); i++)
        for (i = c.begin(); i != c.end(); i++)
          {

            out << *i << " ";
            //out << c[i] << " ";
        }
      return out;
    }

#undef DATA
#undef GET

#endif


  }
}

#endif
