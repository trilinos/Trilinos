#ifndef STK_UTIL_UTIL_FArrayPrint_h
#define STK_UTIL_UTIL_FArrayPrint_h

#include <ostream>
#include <iomanip>

#include <stk_util/diag/FArray.hpp>


namespace sierra {

namespace {

template<unsigned N> struct ArrayVerbosePrint ;

template<>
struct ArrayVerbosePrint<0>
{
  template<typename T>
  static std::ostream &dump(std::ostream &os, const unsigned * const, const T *, const unsigned * const)
  {
    return os;
  }
};

template<>
struct ArrayVerbosePrint<1>
{
  template<typename T>
  static std::ostream &dump(std::ostream &os, const unsigned * const dim,
                            const T *a_ptr, const unsigned * const a_inc)
  {
    const unsigned inc = *a_inc ;
    const T * ptr = a_ptr;
    const T * const end = ptr + *dim * inc ;
    os << "(";
    while ( ptr != end ) {
      if (ptr != a_ptr)
        os << " ";
      os << *ptr;
      ptr += inc ;
    }
    os << ")";

    return os;
  }
};

template<>
struct ArrayVerbosePrint<2>
{
public:
  template<typename T>
  static std::ostream &dump(std::ostream &os, const unsigned * const dim,
                            const T *a_ptr, const unsigned * const a_inc)
  {
    if (dim[0] < 8) { // less than 8 columns wide
      const unsigned inc = a_inc[1] ;
      const T * r_ptr = a_ptr;
      const T * const end = r_ptr + dim[1] * inc ;
      os << "+-\n";
      while ( r_ptr != end ) {
        {
          const unsigned inner_inc = *a_inc ;
          const T * c_ptr = r_ptr;
          const T * const inner_end = c_ptr + dim[0] * inner_inc ;
          os << "| ";
          while ( c_ptr != inner_end ) {
            if (c_ptr != r_ptr)
              os << " ";
            os << *c_ptr;
            c_ptr += inner_inc ;
          }
          os << "\n";
        }
          
        r_ptr += inc ;
      }
      os << "+-\n";
    }
    else {
      const unsigned inc = a_inc[1] ;
      const T * ptr = a_ptr;
      const T * const end = ptr + dim[1] * inc ;
      while ( ptr != end ) {
        ArrayVerbosePrint<1>::dump(os, dim, ptr, a_inc );
        os << std::endl;
        ptr += inc ;
      }
    }

    return os;
  }
};

template<>
struct ArrayVerbosePrint<3>
{
public:
  template<typename T>
  static std::ostream &dump(std::ostream &os, const unsigned * const dim,
                            const T *a_ptr, const unsigned * const a_inc)
  {
//       if (a_inc[3] < 50) {
// 	const unsigned ia = a_inc[2] ;
// 	unsigned index = 0;
// 	const T * ptr = a_ptr;
// 	const T * const end = ptr + a_inc[3];
// 	while ( ptr != end ) {
// 	  if (ptr != a_ptr)
// 	    os << ",";
// 	  os << "(";
// 	  ArrayVerbosePrint<2>::dump(os, dim, ptr, a_inc );
// 	  os << ")";
// 	  ptr += ia ;
// 	  ++index;
// 	}
//       }
//       else {
    const unsigned ia = a_inc[2] ;
    unsigned index = 0;
    const T * const a_end = a_ptr + a_inc[3];
    while ( a_end != a_ptr ) {
      os << "(";
      for (unsigned i = 0; i < 2; ++i)
        os << "0:" << dim[i] - 1 << ", ";
      os << index << ")" << std::endl;
      ArrayVerbosePrint<2>::dump(os, dim, a_ptr, a_inc );
      os << std::endl << std::endl;
      a_ptr += ia ;
      ++index;
    }
//       }

    return os;
  }
};

template<unsigned N>
struct ArrayVerbosePrint
{
public:
  template<typename T>
  static std::ostream &dump(std::ostream &os, const unsigned * const dim,
                            const T *a_ptr, const unsigned * const a_inc)
  {
    const unsigned ia = a_inc[N - 1] ;
    unsigned index = 0;
    const T * const a_end = a_ptr + a_inc[N];
    while ( a_end != a_ptr ) {
      os << "(";
      for (unsigned i = 0; i < N - 1; ++i)
        os << "0:" << dim[i] - 1 << ", ";
      os << index << ")" << std::endl;
      ArrayVerbosePrint<N - 1>::dump(os, dim, a_ptr, a_inc );
      os << std::endl << std::endl;
      a_ptr += ia ;
      ++index;
    }

    return os;
  }
};

} // namespace <unnamed>

template< class ElementType,
	  class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7 >
std::ostream &
operator<<(
  std::ostream &        os,
  const sierra::Array<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > &array)
{
  typedef sierra::Array<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > X;

  ArrayVerbosePrint<X::NumDim>::dump(os, array.dimension(), array.ptr(), array.stride());
  os << std::endl;

  return os;
}


template< class ElementType,
	  class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7 >
std::ostream &
operator<<(
  std::ostream &        os,
  const sierra::ArrayContainer<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > &array)
{
  typedef sierra::ArrayContainer<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > X;

  ArrayVerbosePrint<X::NumDim>::dump(os, array.dimension(), array.ptr(), array.stride());
  os << std::endl;

  return os;
}


template< class ElementType, int Dimension>
std::ostream &
operator<<(
  std::ostream &        os,
  const sierra::FArray<ElementType, Dimension> &array)
{
  typedef sierra::FArray<ElementType, Dimension> X;

  ArrayVerbosePrint<X::NumDim>::dump(os, array.dimension(), array.ptr(), array.stride());
  os << std::endl;

  return os;
}


template< class ElementType, int Dimension>
std::ostream &
operator<<(
  std::ostream &        os,
  const sierra::FArrayContainer<ElementType, Dimension> &array)
{
  typedef sierra::FArrayContainer<ElementType, Dimension> X;

  ArrayVerbosePrint<X::NumDim>::dump(os, array.dimension(), array.ptr(), array.stride());
  os << std::endl;

  return os;
}

} // namespace sierra

#endif // STK_UTIL_UTIL_FArrayPrint_h
