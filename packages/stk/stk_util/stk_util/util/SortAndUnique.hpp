#ifndef stk_util_util_SortAndUnique_hpp
#define stk_util_util_SortAndUnique_hpp

namespace stk
{
namespace util
{

template<typename VECTOR, typename COMPARE>
void sort_and_unique(VECTOR &vector, COMPARE compare )
{
    std::sort(vector.begin(), vector.end(), compare);
    auto endIter = std::unique(vector.begin(), vector.end());
    vector.resize(endIter - vector.begin());
}

template<typename VECTOR>
void sort_and_unique(VECTOR &vec)
{
    sort_and_unique(vec,std::less<typename VECTOR::value_type>());
}

template<class VECTOR, typename COMPARE>
bool is_sorted_and_unique(const VECTOR& vec, COMPARE compare)
{
    bool sorted_and_unique = true;
    for(size_t i=1; i<vec.size(); ++i) {
        if (!compare(vec[i-1],vec[i])) {
            sorted_and_unique = false;
        }
    }
    return sorted_and_unique;
}

template<class VECTOR>
bool insert_keep_sorted_and_unique(typename VECTOR::value_type p, VECTOR& procs)
{
  typename VECTOR::iterator iter = std::lower_bound(procs.begin(), procs.end(), p);
  if (iter == procs.end() || *iter != p) {
    procs.insert(iter, p);
    return true;
  }
  return false;
}

} //namespace util
} //namespace stk

#endif
