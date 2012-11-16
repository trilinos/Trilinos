#ifndef STK_UTIL_UTIL_Fortran_h
#define STK_UTIL_UTIL_Fortran_h

#if ! defined(SIERRA_FORTRAN) && ! defined(SIERRAFORTRAN)

#if defined(FORTRAN_NO_UNDERSCORE)
# define SIERRA_FORTRAN(subname) subname
# define SIERRAFORTRAN(subname) subname
# define SIERRA_FORTRAN_SUFFIX ""
#elif defined(FORTRAN_ONE_UNDERSCORE)
# define SIERRA_FORTRAN(subname) subname##_
# define SIERRAFORTRAN(subname) subname##_
# define SIERRA_FORTRAN_SUFFIX "_"
#elif defined(FORTRAN_TWO_UNDERSCORES)
# define SIERRA_FORTRAN(subname) subname##__
# define SIERRAFORTRAN(subname) subname##__
# define SIERRA_FORTRAN_SUFFIX "__"
#endif

#endif // SIERRA_FORTRAN

#endif // STK_UTIL_UTIL_Fortran_h
