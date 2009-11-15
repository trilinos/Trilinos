INCLUDE(CheckCXXSourceCompiles)
INCLUDE(MultilineSet)

FUNCTION(CHECK_FOR_BOOST_IS_POLYMORPHIC  VARNAME)

  SET(SOURCE
  "
  #include <boost/type_traits/is_polymorphic.hpp>
  
  int main()
  {
     // If this compiles, I am happy
     return (boost::is_polymorphic<int>::value == true);
  }
  "
  )
  
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()
