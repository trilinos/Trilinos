
/*

Base file created and maintained by hand to handle all of the common options
by all other XXX_config.h files.

2008/07/03: rabartl: I created this file as follows:

1) I copied the defines for F77_FUNC(...), F77_FUNC_(...) from
/sierra/Release/Trilinos/8.0.6/include/Teuchos_config.h.

I then created the rest of defines as follows:

1) First, I did a cat of all the config files in
/sierra/Release/Trilinos/8.0.6/include and put them in one file as:

  [rabartl@sais503 include]$ pwd
  /sierra/Release/Trilinos/8.0.6/include
  [rabartl@sais503 include]$ cat *_Config.h *_config.h > /var/scratch2/rabartl/PROJECTS/Sierra/Aria_Trilinos/Nbtools/Trilinos/dev/config_headers/glob.out

2) I removed the comments, sorted, and removed duplicates by doing:

  [rabartl@sais503 config_headers]$ pwd
  /var/scratch2/rabartl/PROJECTS/Sierra/Aria_Trilinos/Nbtools/Trilinos/dev/config_headers
  [rabartl@sais503 config_headers]$ cat glob.out | grep '#define' | sort > glob.clean.sorted.out 
  [rabartl@sais503 config_headers]$ cp glob.clean.sorted.out glob.clean.sorted.tmp
  [rabartl@sais503 config_headers]$ cat glob.clean.sorted.out | uniq -c > glob.clean.uniqu.out

3) I then went through and manually removed all of the package/configuration
specific items to form TrilinosPlatform_config.h

 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#ifndef F77_FUNC
#  ifndef _AIX
#    define F77_FUNC(name,NAME) name ## _
#  else
#    define F77_FUNC(name,NAME) name
#  endif
#endif

/* As F77_FUNC, but for C identifiers containing underscores. */
#ifndef F77_FUNC_
#  ifndef _AIX
#    define F77_FUNC_(name,NAME) name ## _
#  else
#    define F77_FUNC_(name,NAME) name
#  endif
#endif

/* Define the Fortran name mangling to be used for the BLAS */
#define F77_BLAS_MANGLE(name,NAME) F77_FUNC(name,NAME)

#define HAVE_ALGORITHM 1
#define HAVE_ALGO_H 1
#define HAVE_BOOL 
#define HAVE_CASSERT 1
#define HAVE_CCTYPE
#define HAVE_CERRNO 1
#define HAVE_CFLOAT 1
#define HAVE_CFUNC 
#define HAVE_CLIMITS 1
#define HAVE_CMATH 1
#define HAVE_COMPLEX 1
#define HAVE_CSTDARG 1
#define HAVE_CSTDIO 1
#define HAVE_CSTDLIB 1
#define HAVE_CSTRING 1
#define HAVE_CTIME 1
#define HAVE_EXPORT_MAKEFILES 
#define HAVE_FPU_CONTROL_H 1
#define HAVE_FSTREAM 1
#define HAVE_GCC_ABI_DEMANGLE 
#define HAVE_GNUMAKE 

#if ! defined(__sun) && ! defined(__sgi)
#  define HAVE_INF_SUPPORT 
#endif

#define HAVE_INOUT 
#define HAVE_INTTYPES_H 1
#define HAVE_IOMANIP 1
#define HAVE_IOSTREAM 1
#define HAVE_ITERATOR 1
#define HAVE_LIST 1
#define HAVE_MALLOC_H 1
#define HAVE_MAP 1
#define HAVE_MATH_H 1
#define HAVE_MEMORY 1
#define HAVE_MEMORY_H 1
#define HAVE_MPI 
#define HAVE_MUTABLE 
#define HAVE_NAMESPACES 
#define HAVE_NEW_FOR_SCOPING 
#define HAVE_NUMERIC 1
#define HAVE_NUMERIC_LIMITS 
#define HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS 
#define HAVE_SET 1
#define HAVE_SSTREAM 1
#define HAVE_STDEXCEPT 1
#define HAVE_STDINT_H 1
#define HAVE_STDIO_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STD_IOS_BASE_FMTFLAGS 
#define HAVE_STD_NEW_COUNT_SYNTAX 
#define HAVE_STD_SPRINTF 
#define HAVE_STL 
#define HAVE_STRING 1
#define HAVE_STRINGS_H 1
#define HAVE_STRING_H 1
#define HAVE_SYS_RESOURCE_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_SYS_TIME_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_SYS_UTSNAME_H 1
#define HAVE_TIME_H 1
#define HAVE_TYPEINFO 1
#define HAVE_UNISTD_H 1
#define HAVE_VALGRIND_VALGRIND_H 1
#define HAVE_VECTOR 1
#define INVALID_TEMPLATE_QUALIFIER 
#define STDC_HEADERS 1

#if defined(_AIX)
#  define TEMPLATE_FRIENDS_NOT_SUPPORTED
#  define TEUCHOS_PRIVIATE_DELETE_NOT_SUPPORTED
#endif

#if defined(__PGI)
#  define THYRA_DEFAULT_PRODUCT_VECTOR_SPACE_EXPLICIT_INSTANTIATION
#  define THYRA_DEFAULT_PRODUCT_VECTOR_EXPLICIT_INSTANTIATION
#  define THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_EXPLICIT_INSTANTIATION
#endif

#if defined(__sun) 
#  define THYRA_INJECT_USING_DECLARATIONS
#endif
