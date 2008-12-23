INCLUDE(CheckCXXSourceRuns)
INCLUDE(MultilineSet)

FUNCTION(CHECK_FOR_ABI_DEMANGLE VARNAME)


  SET(SOURCE
  "
  #include <cxxabi.h>
  #include <string>
  namespace MyNamespace {
    class MyClass {};
  }
  
  int main()
  {
    const std::string
      mangledName = typeid(MyNamespace::MyClass).name();
    int status;
    char
      *_demangledName = abi::__cxa_demangle(mangledName.c_str(),0,0,&status);
    const std::string demangledName(_demangledName);
    free(_demangledName);
    return ( demangledName == \"MyNamespace::MyClass\" ? 0 : 1 );
  }
  "
  )
  
  CHECK_CXX_SOURCE_RUNS("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()