#ifndef SIERRA_SIERRA_UTIL_ALGORITHMS_HPP
#define SIERRA_SIERRA_UTIL_ALGORITHMS_HPP

namespace sierra {
namespace util {

template <class ForwardIterator >
bool sorted( ForwardIterator first, ForwardIterator last )
{
  bool is_sorted = true;
  ForwardIterator prev(first);
  ++first;
  for( ; first != last && is_sorted; ++first, ++prev) {
    is_sorted = *prev < *first || !(*first < *prev);
  }
  return is_sorted;
}

template <class ForwardIterator, class Compare >
bool sorted( ForwardIterator first, ForwardIterator last, const Compare & comp = Compare() )
{
  bool is_sorted = true;
  ForwardIterator prev(first);
  ++first;
  for( ; first != last && is_sorted; ++first, ++prev) {
    is_sorted = comp(*prev,*first) || !comp(*first,*prev);
  }
  return is_sorted;
}


} // namespace util
} // namespace sierra


#endif // SIERRA_SIERRA_UTIL_ALGORITHMS_HPP
