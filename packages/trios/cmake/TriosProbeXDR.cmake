
########## CHECK FOR HEADER FILES ############

INCLUDE(CheckIncludeFiles)

# Probe for header files
CHECK_INCLUDE_FILES("mach/mach_time.h" HAVE_TRIOS_MACH_TIME_H)
CHECK_INCLUDE_FILES("time.h" HAVE_TRIOS_TIME_H)
CHECK_INCLUDE_FILES("sys/time.h" HAVE_TRIOS_SYS_TIME_H)
CHECK_INCLUDE_FILES("papi.h" HAVE_TRIOS_PAPI_H)

########## CHECK FOR FUNCTIONS ############

INCLUDE(CheckLibraryExists)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXSourceCompiles)

CHECK_FUNCTION_EXISTS(xdr_u_int16_t HAVE_TRIOS_XDR_U_INT16_T)
CHECK_FUNCTION_EXISTS(xdr_u_int32_t HAVE_TRIOS_XDR_U_INT32_T)
CHECK_FUNCTION_EXISTS(xdr_u_int64_t HAVE_TRIOS_XDR_U_INT64_T)

# XDR_SIZEOF
check_c_source_compiles(
    "#include <rpc/xdr.h>\nint main(){xdr_sizeof(NULL,NULL);return 0;}"
    HAVE_TRIOS_XDR_SIZEOF
)
