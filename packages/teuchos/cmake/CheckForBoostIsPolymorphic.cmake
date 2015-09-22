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
  
  SET(CMAKE_REQUIRED_INCLUDES ${TPL_Boost_INCLUDE_DIRS})
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()
