#ifndef STK_UTIL_UTIL_FortranMod_Hpp
#define STK_UTIL_UTIL_FortranMod_Hpp

#if ! defined(FORTRAN_MODULE)
// defined in Jamfile #define LINUX_INTEL_F90_MODULES
#ifdef LINUX_GCC_F90_MODULES
  #define FORTRAN_MODULE(module_name,subroutinename) __ ## module_name ## _MOD_ ## subroutinename
#elif defined(LINUX_INTEL_F90_MODULES)
  #define FORTRAN_MODULE(module_name,subroutinename)  module_name ## _mp_ ## subroutinename ## _
#else
  #define FORTRAN_MODULE(module_name,subroutinename) no_module_defs_for_this_compiler ## module_name
#endif

#endif // FORTRAN_MODULE

#endif  // end STK_UTIL_UTIL_FortranMod_Hpp
