#ifndef STK_UTIL_RADIX_SORT_2_HPP
#define STK_UTIL_RADIX_SORT_2_HPP

#include <stk_util/util/ForceInline.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/array.hpp>

#include <inttypes.h>

namespace stk { namespace details {

template <typename Key>
FORCEINLINE
typename boost::enable_if_c< boost::is_unsigned<Key>::value, bool >::type
is_negative(const Key &)
{ return false; }

template <typename Key>
FORCEINLINE
typename boost::enable_if_c< !boost::is_unsigned<Key>::value, bool  >::type
is_negative(const Key & k)
{ return (k < (Key)0); }

inline
bool is_little_endian()
{
  const uint32_t i=0x01020304u;
  unsigned char const * const c = reinterpret_cast<unsigned char const * const>(&i);
  return *c == static_cast<unsigned char>(4u);
}


template <bool IsLittleEndian> struct GetRadix;

template <> struct GetRadix<true>
{
  template <typename Key>
  FORCEINLINE
  unsigned operator()( Key const& key, unsigned byte ) const
  {
    return reinterpret_cast<unsigned char const * const>(&key)[byte];
  }

};

template <> struct GetRadix<false>
{
  template <typename Key>
  unsigned operator()( Key const& key, unsigned byte ) const
  {
    return reinterpret_cast<unsigned char const * const>(&key)[sizeof(Key) - byte -1];
  }

};

template < typename KeyArray , bool IsLittleEndian >
struct RadixSort
{
  typedef KeyArray key_array;
  typedef typename KeyArray::value_type key_type;

  BOOST_STATIC_ASSERT(( boost::is_integral<key_type>::value || boost::is_enum<key_type>::value || boost::is_floating_point<key_type>::value ));


  static const unsigned key_size = sizeof(key_type);
  static const unsigned radix_bits = 8;
  static const unsigned num_bytes = ((key_size*8) + (radix_bits - 1))/radix_bits;
  static const unsigned radix_size = 1u << radix_bits;
  static const unsigned radix_mask = radix_size - 1u;

  typedef boost::array< boost::array<size_t, radix_size>, num_bytes> histogram_type;
  typedef histogram_type offset_type;
  typedef boost::array<bool,num_bytes> skip_type;

  struct histogram_result
  {
    size_t unsorted;
    size_t num_negatives;
  };

  // compute the histogram
  // returns true if the array is already sorted
  histogram_result compute_histogram( KeyArray const & keys
                                     ,histogram_type & histogram
                                    ) const
  {
    GetRadix<IsLittleEndian> get_radix;

    histogram_result results = {0,0};

    for (size_t i=0, ie=keys.size(); i<ie; ++i) {
      results.unsorted += (i>0) ? (keys[i] <= keys[i-i]) : 0;
      results.num_negatives += is_negative(keys[i]);
      for (unsigned byte=0; byte<num_bytes; ++byte) {
        unsigned radix = get_radix(keys[i],byte);
        ++histogram[byte][radix];
      }
    }

    return results;
  }

  void compute_offsets_unsigned( const histogram_type & histogram
                                ,const size_t size
                                ,offset_type & offsets
                                ,skip_type & skip_byte
                               ) const
  {
    for (unsigned byte=0; byte<num_bytes; ++byte) {
      for (unsigned radix=1; radix<radix_size; ++radix) {
        // if all the keys are in the same radix
        // we can skip the byte
        const size_t hist_count = histogram[byte][radix-1];
        skip_byte[byte] |= (hist_count == size);

        // starting offset for radix => prefix sum
        offsets[byte][radix] = offsets[byte][radix-1] + hist_count;
      }
      skip_byte[byte] |= histogram[byte][radix_size-1] == size;
    }
  }

  void compute_offsets_signed( const histogram_type & histogram
                              ,const size_t size
                              ,offset_type & offsets
                              ,skip_type & skip_byte
                             ) const
  {
    for (unsigned byte=0; byte<num_bytes-1; ++byte) {
      for (unsigned radix=1; radix<radix_size; ++radix) {
        // if all the keys are in the same radix
        // we can skip the byte
        const size_t hist_count = histogram[byte][radix-1];
        skip_byte[byte] |= (hist_count == size);

        // starting offset for radix => prefix sum
        offsets[byte][radix] = offsets[byte][radix-1] + hist_count;
      }
      skip_byte[byte] |= histogram[byte][radix_size-1] == size;
    }
    {
      const unsigned byte = num_bytes-1;
      // move negative integers before positive integers
      for (unsigned radix=1; radix<radix_size; ++radix) {
        // if all the keys are in the same radix
        // we can skip the byte
        const int prev = ((radix_size>>1)+radix-1u) & radix_mask;
        const int curr = ((radix_size>>1)+radix) & radix_mask;
        const size_t hist_count = histogram[byte][prev];
        skip_byte[byte] |= (hist_count == size);

        // starting offset for radix => prefix sum
        offsets[byte][curr] = offsets[byte][prev] + hist_count;
      }
      skip_byte[byte] |= histogram[byte][(radix_size>>1)] == size;
    }
  }

  // this will only work with floats that conform to the IEEE 734 standard
  void compute_offsets_signed_float( const histogram_type & histogram
                                    ,const size_t size
                                    ,const size_t num_negatives
                                    ,offset_type & offsets
                                    ,skip_type & skip_byte
                                   ) const
  {
    for (unsigned byte=0; byte<num_bytes-1; ++byte) {
      for (unsigned radix=1; radix<radix_size; ++radix) {
        // if all the keys are in the same radix
        // we can skip the byte
        const size_t hist_count = histogram[byte][radix-1];
        skip_byte[byte] |= (hist_count == size);

        // starting offset for radix => prefix sum
        offsets[byte][radix] = offsets[byte][radix-1] + hist_count;
      }
      skip_byte[byte] |= histogram[byte][radix_size-1] == size;
    }
    {
      const unsigned byte = num_bytes-1;
      // move negative integers before positive integers
      offsets[byte][0] = num_negatives;
      for (unsigned radix=1; radix<=(radix_size>>1); ++radix) {
        // if all the keys are in the same radix
        // we can skip the byte
        const size_t hist_count = histogram[byte][radix-1];
        skip_byte[byte] |= (hist_count == size);

        // starting offset for radix => prefix sum
        offsets[byte][radix] = offsets[byte][radix-1] + hist_count;
      }
      offsets[byte][radix_size-1] = 0;

      for (unsigned inv_radix=radix_size-2; inv_radix>(radix_size>>1); --inv_radix) {
        // if all the keys are in the same radix
        // we can skip the byte
        const size_t hist_count = histogram[byte][inv_radix+1];
        skip_byte[byte] |= (hist_count == size);

        // starting offset for radix => prefix sum
        offsets[byte][inv_radix] = offsets[byte][inv_radix+1] + hist_count;
      }
      skip_byte[byte] |= histogram[byte][(radix_size>>1)] == size;
      skip_byte[byte] |= histogram[byte][(radix_size-1)]  == size;
    }
  }

  void sort_keys(key_array & keys, offset_type & offsets, const skip_type & skip_byte) const
  {
    GetRadix<IsLittleEndian> get_radix;

    const size_t size = keys.size();
    key_array scratch(size);

    for (unsigned byte=0; byte<num_bytes; ++byte) {
      if (!skip_byte[byte]) {
        for (size_t i=0; i<size; ++i) {
          unsigned radix = get_radix(keys[i],byte);
          scratch[offsets[byte][radix]++] = keys[i];
        }
        // swap buffers
        std::swap(keys,scratch);
      }
    }
  }

  void sort_keys_signed_floats(key_array & keys, offset_type & offsets, const skip_type & skip_byte) const
  {
    GetRadix<IsLittleEndian> get_radix;

    const size_t size = keys.size();
    key_array scratch(size);

    for (unsigned byte=0; byte<num_bytes-1; ++byte) {
      if (!skip_byte[byte]) {
        for (size_t i=0; i<size; ++i) {
          unsigned radix = get_radix(keys[i],byte);
          scratch[offsets[byte][radix]++] = keys[i];
        }
        // swap buffers
        std::swap(keys,scratch);
      }
    }
    {
      unsigned byte = num_bytes-1;
      if (!skip_byte[byte]) {
        for (size_t i=0; i<size; ++i) {
          unsigned radix = get_radix(keys[i],byte);
          size_t index = radix <= (radix_size>>1) ? offsets[byte][radix]++ : --offsets[byte][radix-1];
          scratch[index] = keys[i];
        }
        // swap buffers
        std::swap(keys,scratch);
      }
    }
  }

  void operator()(key_array & keys) const
  {
    const size_t size = keys.size();

    offset_type offsets = {};
    skip_type skip_byte = {};

    bool has_negative_float = false;

    // compute histogram and offsets
    {
      histogram_type histogram = {};
      const histogram_result r = compute_histogram(keys,histogram);

      if (r.unsorted == 0u) {
        // already sorted
        return;
      }
      else if (r.unsorted == size) {
        // reversed sorted
        for (size_t i = 0, j= size-1; i < j; ++i, --j) {
          std::swap(keys[i],keys[j]);
        }
        return;
      }

      if (!r.num_negatives) {
        // no negative values present
        compute_offsets_unsigned(histogram, size, offsets, skip_byte);
      }
      else if(boost::is_integral<key_type>::value || boost::is_enum<key_type>::value) {
        // negative integers
        compute_offsets_signed(histogram, size, offsets, skip_byte);
      }
      else {
        // negative floats
        has_negative_float = true;
        compute_offsets_signed_float(histogram, size, r.num_negatives, offsets, skip_byte);
      }
    }

    if (!has_negative_float) {
      sort_keys(keys, offsets, skip_byte);
    }else {
      sort_keys_signed_floats(keys, offsets, skip_byte);
    }
  }

};

}} //namespace stk::details

namespace stk {

template <typename KeyArray>
typename boost::enable_if_c<(   boost::is_integral<typename KeyArray::value_type>::value
                             || boost::is_enum<typename KeyArray::value_type>::value
                             || boost::is_floating_point<typename KeyArray::value_type>::value
                            ),void>::type
radix_sort( KeyArray & keys)
{
  if (details::is_little_endian()) {
    details::RadixSort<KeyArray, true> rsort;
    rsort(keys);
  }
  else {
    details::RadixSort<KeyArray, false> rsort;
    rsort(keys);
  }
}

} // namespace stk


#endif // STK_UTIL_RADIX_SORT_2_HPP
